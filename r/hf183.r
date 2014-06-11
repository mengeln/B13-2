source("r/helpers.r")

hf <- function (samplefile, sketafile, format, org) {

  
  # Read in data
  sampledata <- readAndClean(samplefile, format)
  sketaData <- readAndClean(sketafile$datapath, format)
  metadata <- getMeta(samplefile, format)
  
  # Subset by target
  HFData <- sampledata[sampledata$Target == "hf183", ]
  IACdata <- sampledata[sampledata$Target == "iac", ]
  
  
  # Standard Curve

  HFStandard <- HFData[grepl("Std", HFData$Content), ]
  HFStandard$Log10CopyPeruL <- log10(as.numeric(HFStandard$CopyPeruL)) 
  
  HF.model <- lm(data = HFStandard,  Cq ~ Log10CopyPeruL) # create standard curve model
  HF.yint <- coef(HF.model)[[1]]
  HF.Slope <- coef(HF.model)[[2]]
  HF.r2 <- summary(HF.model)$r.squared
  HF.Efficiency <- 10^(-1/HF.Slope)
  
  
  # Make report 
  stdReport <- standardQC(10^(-1/coef(HF.model)[[2]]),
                             summary(HF.model)$r.squared,
                             "HF183")
  
  # Controls
  controlFrame <- function (data, assay) {
    data$Sample <- toupper(data$Sample)
    cData <- data[grepl("NTC|NEC", data$Sample),] 
    cData <- ddply(cData, .(Sample), function(df){
      df$Replicate <- paste0("Ct$_{Rep", 1:nrow(df), "}$")
      df
    })
    controlDF <- dcast(cData, Sample ~ Replicate, value.var="Cq")
    controlDF$"PASS?" <- apply(controlDF[, -1], 1, function(x)ifelse(all(x == m), "PASS", "FAIL"))
    controlDF[controlDF ==m] <- "N/A"
    controlDF$Assay <- assay
    controlDF
  }
  controlDF <- controlFrame(HFData, "HF183")
  controlSk <- controlFrame(sketaData, "Sketa22")
  
  controlsDF <- rbind(controlDF, controlSk[controlSk$Sample == "NTC",])
  numCols <- ncol(controlsDF)
  controlsDF <- controlsDF[, c(numCols, 1:(numCols - 1))]
  
  dbCtrl <- melt(controlsDF[, names(controlsDF) %nin% "PASS?"], id.vars=c("Assay", "Sample"))
  dbCtrl$variable <- as.numeric(gsub(".*?Rep(\\d).*", "\\1", dbCtrl$variable))
  names(dbCtrl)[names(dbCtrl) %in% c("Assay", "Sample", "variable", "value")] <- c("Target", "Type", "Rep", "Ct")
  # Sketa Inhibition QC
  
  sketaQC <- function(data=sketaData, threshold=thres){
    sk.unkn <- data[grepl("Unkn", data$Content), ]
    calibrators <- sketaData$Cq[grepl("NEC", toupper(sketaData$Sample))]
    Ct.sk.calibrator <- mean(calibrators)
    Ct.sk.sd <<- sd(calibrators)
    
    sk.calibrator <<- Ct.sk.calibrator
    sk.unkn$sk.dct <- sk.unkn$Cq - Ct.sk.calibrator
    sk.unkn$Inhibition <- ifelse(sk.unkn$sk.dct > threshold | sk.unkn$sk.dct < (-threshold),  # need to update variable name to "SketaQC"
                                 "FAIL", "PASS")
    names(sk.unkn)[names(sk.unkn)=="Cq"] <- "sk.Ct"
    sk.unkn
  }
  
  sketaData <- sketaQC(sketaData)
  
  sketaDataTrim <- sketaData[, c("Sample", "sk.Ct", "sk.dct", "Inhibition")]
  
  sketaDataTrim <- ddply(sketaDataTrim, .(Sample), function(df){
    data.frame(Sample = unique(df$Sample),
               sk.Ct = mean(df$sk.Ct, na.rm=TRUE),
               sk.dct = mean(df$sk.dct, na.rm=TRUE),
               Inhibtion = ifelse(all(df$Inhibition == "PASS"), "PASS", "FAIL")
    )
  })
  names(sketaDataTrim) <- c("Sample", "sketa Ct$_{mean}$", "$\\Delta$Ct$_{mean}$", "Pass?")
  
  # IAC Inhibition QC
  
  IACNEC <- IACdata[grepl("NEC", IACdata$Sample), c("Sample", "Cq")]
  names(IACNEC) <- c("IAC NEC", "Ct")
  
  IACstd <- IACdata[IACdata$Content == "Std", ]
  IACstd <- IACstd[order(IACstd$CopyPeruL, decreasing =FALSE), ]
  Istd <- ddply(IACstd, .(CopyPeruL), summarize, Cq = mean(Cq, na.rm=TRUE))
  Istd$compare <- Istd$Cq - Istd$Cq[Istd$CopyPeruL == min(Istd$CopyPeruL)]
  ROQ <- Istd$Cq[Istd$compare < 0.75 & Istd$compare > - 0.75]
  Max.comp <- max(Istd$CopyPeruL[Istd$compare < 0.75 & Istd$compare > - 0.75])
  
  IACcompetition <- predict(HF.model, data.frame(Log10CopyPeruL = log10(Max.comp)))
  
  iacsd <- sd(ROQ, na.rm=TRUE)
  IACinterference <- mean(ROQ) + 4*ifelse(is.na(iacsd), 0, iacsd)
  
  IACinhib <- ddply(IACdata[IACdata$Content == "Unkn", ], .(Sample), function(df){ 
    mean <- mean(df$Cq, na.rm=TRUE)
    inhibited <- mean >= IACinterference 
    data.frame(Sample = unique(df$Sample),
               Ctmean = mean,
               "PASS?" = ifelse(!inhibited, "PASS", "FAIL")
    )
  })
  IACinhib$Ctmean[IACinhib$Ctmean == m] <- "ND"
  names(IACinhib) <- c("Sample", "IAC Ct$_{mean}$", "Pass?")
  
  # CCE interpolation
  directCT <- function(data, ulPerRxn=2, mlFiltered=100, ulCE=600, ulCErecovered=380, ulPE=100){
    data$r2 <- HF.r2
    data$Eff <- HF.Efficiency
    data$Slope <- HF.Slope
    data$yint <- HF.yint
    data$HF.predict <- (data$Cq - HF.yint) / HF.Slope  # direct interpolation from standard curve
    data$copiesPerRxn <- 10^(data$HF.predict)
    data$copiesPer100ml <- data$copiesPerRxn/ulPerRxn * ulPE * (ulCE/ulCErecovered) * 100/mlFiltered
    data$log10copiesPer100ml <- round(log10(data$copiesPer100ml), digits=3)
    
    data[data$Cq == m, c("log10copiesPer100ml")] <- "ND"
    
    data[!grepl("Std", data$Content), ]
    
  }
  
  HFData2 <- directCT(HFData)
  if(nrow(HFData2) > 0) {
    HFData2$Date <- metadata["Run Started"]
  }
  
  
  
  ### Integrate results ###
  result <- Reduce(function(x,y)merge(x,y, by="Sample"), list(HFData2, sketaDataTrim, IACinhib))
  
  result$Competition <- result$"Pass?.y" == "FAIL" & result$Cq < IACcompetition  # need to modify in the future for accidental overdose of iac
  IACinhib$Competition <- ifelse(result$Competition[match(IACinhib$Sample, result$Sample)], "Yes", NA)
  
  resultsTrim <- rbind.fill(lapply(split(result, result$Sample), function(df){
    inhibition <- !all(rbind(df$"Pass?.x" == "PASS", df$"Pass?.y" == "PASS")) & !df$Competition
    
    res <- df[, c("Sample", "Target", "Cq", "log10copiesPer100ml", "copiesPer100ml")]
    res$Mean <- NA
    res$Mean[1] <- ifelse(any(inhibition), "inhibited", round(mean(log10(res$copiesPer100ml)), digits=2))
    res$Cq[res$Cq == m] <- "N/A"
    res$Replicate <- 1:nrow(res)
    if(any(inhibition))res$log10copiesPer100ml <- "inhibited"
    subset(res, select=-c(copiesPer100ml))
  }))
  resultsTrim <- arrange(resultsTrim, Sample, Mean)
  
  names(resultsTrim)[3] <- c("Ct")
  
  rmelt <- melt(resultsTrim, id.vars=c("Sample", "Target", "Replicate"))
  resultsTrim2 <- dcast(rmelt, Sample + Target  ~ Replicate + variable, value.var="value")
  
  resultsTrim2 <- resultsTrim2[, c("Sample", "Target", "1_Ct",
                                   "2_Ct", "3_Ct", "1_log10copiesPer100ml",
                                   "2_log10copiesPer100ml",
                                   "3_log10copiesPer100ml",
                                   "1_Mean")]
  resultsTrim2[is.na(resultsTrim2$"1_Ct") & is.na(resultsTrim2$"2_Ct"), c("1_Mean")] <- "ND"
  names(resultsTrim2) <- c("Sample", "Target",
                           "Ct$_{Rep 1}$",
                           "Ct$_{Rep 2}$",
                           "Ct$_{Rep 3}$",
                           "$\\log_{10}$ copies/100 \\si{\\milli\\litre}$_{Rep1}$",
                           "$\\log_{10}$ copies/100 \\si{\\milli\\litre}$_{Rep2}$",
                           "$\\log_{10}$ copies/100 \\si{\\milli\\litre}$_{Rep3}$",
                           "Mean $\\log_{10}$ copies/100 \\si{\\milli\\litre}")
  print(stdReport)
  list(metadata=metadata,
       org=org,
       r2.min=r2.min,
       eff.min=eff.min,
       eff.max=eff.max,
       stdTable = stdReport,
       controlsDF=controlsDF,
       sk.calibrator=sk.calibrator,
       Ct.sk.sd=Ct.sk.sd,
       thres=thres,
       sketaDataTrim=sketaDataTrim,
       ROQ=ROQ,
       IACcompetition=IACcompetition,
       IACinhib=IACinhib,
       IACNEC=IACNEC,
       HF.Slope=HF.Slope,
       m=m,
       resultsTrim2=resultsTrim2)
  
}