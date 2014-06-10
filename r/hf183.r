source("r/helpers.r")

hf183 <- function (samplefile, sketafile) {
  
  eff.max <- 2.1
  eff.min <- 1.87
  r2.min <- 0.98
  m <- 45  # Ct to assign to unamplified wells
  thres <- 3.0 # sketa sample processing control failure threshold
  
  # Read in data
  sampledata <- readAndClean(samplefile)
  sketadata <- readAndClean(sketafile)
  
  # Subset by target
  hf183Data <- sampledata[sampledata$Target == "hf183", ]
  IACdata <- sampledata[sampledata$Target == "iac", ]
  
  
  ### Standard Curve QC ###
  
  # Build linear model
  hf183Std <- hf183Data[grepl("Std", hf183Data$Content), ]
  hf183Std$Log10CopyPeruL <- log10(as.numeric(hf183Std$CopyPeruL)) 
  hf183Model <- lm(data = hf183Std,  Cq ~ Log10CopyPeruL)

  # Make report 
  stdReport <- standardQC(10^(-1/coef(hf183Model )[[2]]),
                         summary(HF.model)$r.squared,
                         "HF183")
  
  ### Negative Controls QC ###

  controlDF <- controlFrame(hf183Data, "HF183")
  controlSk <- controlFrame(sketaData, "Sketa22")
  controlsReport <- rbind(controlDF, controlSk[controlSk$Sample == "NTC",])

  ### Sketa Inhibition QC ###
  
  sketaQReport <- sketaQC(sketaData)

  ### IAC Inhibition QC ###
  
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
}