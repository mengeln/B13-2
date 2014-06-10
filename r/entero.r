source("r/helpers.r")

ent <- function (samplefile, format, org) {

  
  # Read in data
  sampledata <- readAndClean(samplefile, format)
  metadata <- getMeta(samplefile, format)
  
  # Subset by target
  enteroData <- sampledata[sampledata$Target == "ent", ]
  sketaData <- sampledata[grepl("sketa", sampledata$Target), ]
  
  # Standard Curve
  
  # Build linear model
  enteroStd <- enteroData[grepl("Std", enteroData$Content), ]
  enteroStd$Log10CopyPeruL <- log10(as.numeric(enteroStd$CopyPeruL)) 
  enteroModel <- lm(data = enteroStd,  Cq ~ Log10CopyPeruL)
  
  sketaStd <- sketaData[grepl("Std", sketaData$Content), ]
  sketaStd$Log10CopyPeruL <- log10(as.numeric(sketaStd$CopyPeruL)) 
  sketaModel <- lm(data = sketaStd,  Cq ~ Log10CopyPeruL)
  
  # Make report 
  entstdReport <- standardQC(10^(-1/coef(enteroModel)[[2]]),
                          summary(enteroModel)$r.squared,
                          "Entero1A")
  skstdReport <- standardQC(10^(-1/coef(sketaModel)[[2]]),
                             summary(sketaModel)$r.squared,
                            "Sketa")
  stdReport <- rbind(entstdReport, skstdReport)
  
  
  # Negative Control QC
  
  controlDF <- controlFrame(enteroData, "Entero1A")
  controlSk <- controlFrame(sketaData, "Sketa")
  
  controlsReport <- rbind(controlDF, controlSk[controlSk$Sample == "NTC",])

  
  ### Sketa Inhibition QC ###
  NECmean <- mean(sketaData$Cq[grepl("NEC", sketaData$Sample)], na.rm=TRUE)
  calibratorQC <- data.frame(CalibratorCt = sketaData$Cq[grepl("calibrator", sketaData$Sample)])
  calibratorQC$delta <- calibratorQC$CalibratorCt - NECmean
  calibratorQC$PASS <- ifelse(calibratorQC$delta > thres | calibratorQC$delta < (-thres), "FAIL", "PASS")
  names(calibratorQC) <- c("Calibrator Ct", "$\\Delta$ Ct", "QC")

  sketaQCReport <- sketaQC(sketaData)

  # Ct to copy number (dct quantification model)
  
  dct <- function(data, mlFiltered=100, cal=1e5){
    Ct.ent.calibrator <- mean(data$Cq[data$Sample == "calibrator"])

    data$ent.dct <- data$Cq - Ct.ent.calibrator
    data$log10cellPerFilter <- data$ent.dct/coef(enteroModel)[[2]]  + log10(cal)
    data$log10cellPer100ml <- round((data$log10cellPerFilter + log10(100/mlFiltered)), digits=3)
    data$cellPer100ml <- 10^data$log10cellPer100ml
    data[!grepl("Std", data$Content), ]
    
  }

  entData <- dct(enteroData)
  result  <- ddply(entData[entData$Content == "Unkn", ], .(Sample), function(df){
    Cqs <- df$cellPer100ml[df$Cq != m]
    df$Mean <- round(log10(mean(Cqs)), digits=3)
    df
  })
  
  # Inhibited flag
  result$log10cellPer100ml[result$Inhibition == "FAIL"] <- "inhibited"
  
  resultsTrim <- subset(result, select = c(Sample, Target, Cq, log10cellPer100ml, Mean))
  names(resultsTrim)[3:4] <- c("Ct", "log10")
  resultsTrim[resultsTrim$Ct == m, c("Ct", "Mean", "log10")] <- "N/A"
  
  resultsTrim <- ddply(resultsTrim, .(Sample), function(df){
    df$Replicate <- 1:nrow(df)
    if(any(df$log10 == "inhibited"))
      df$Mean[!is.na(df$Mean)] <- "inhibited" 
    df
  })
  resultsTrim <- dcast(melt(resultsTrim,
                            id.vars=names(resultsTrim[names(resultsTrim) %in% c("Sample", "Target", "Replicate")])),
                       Sample + Target ~ Replicate + variable, value.var="value")
  resultsTrim <- resultsTrim[, c("Sample", "Target", "1_Ct", "2_Ct", "1_log10", "2_log10", "1_Mean")]
  resultsTrim$Target <- "Entero1A"
  resultsTrim$"1_Mean"[resultsTrim$Sample %in% sketaQCReport$Sample[sketaQCReport$QC == "FAIL"]] <- "Inhibited"
  names(resultsTrim) <- c("Sample", "Target", "Ct$_{Rep 1}$", "Ct$_{Rep 2}$", "$\\log_{10}$ cells/100 \\si{\\milli\\litre}$_{Rep1}$",
                          "$\\log_{10}$ cells/100 \\si{\\milli\\litre}$_{Rep2}$", "Mean $\\log_{10}$ cells/100 \\si{\\milli\\litre}")
  
  
  list(metadata = metadata,
       org = org,
       
       r2.min = r2.min,
       eff.min = eff.min,
       eff.max = eff.max,
       m = m,
       thres= thres,
       
       ent.r2 = summary(enteroModel)$r.squared,
       ent.Efficiency =  10^(-1/coef(enteroModel)[[2]]),
       ent.Slope = coef(enteroModel)[[2]],
       
       stdTable = stdReport,
       sketa.r2 =  summary(sketaModel)$r.squared,
       sketa.Efficiency = 10^(-1/coef(sketaModel)[[2]]),
       
       controlsDF = controlsReport,
       NECmean = NECmean,
       calibratorQC = calibratorQC,
       calSD = sd(calibratorQC$CalibratorCt),
       
       sketaDataTrim = sketaQCReport,
       sk.calibrator =  mean(sketaData$Cq[sketaData$Sample == "calibrator"]),

       resultsTrim2 = resultsTrim)
}
