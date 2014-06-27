### helpers

library(plyr)
library(reshape2)
library(lubridate)

# Globals
eff.max <- 2.1
eff.min <- 1.87
r2.min <- 0.98
m <- 45  # Ct to assign to unamplified wells
thres <- 1.7 # inhibition threshold
sketaCal <- 10

`%nin%` <- Negate(`%in%`)

errorBlock <- function(message, block){
  test <- try(block)
  if(class(test) == "try-error")
    stop(message, call. = FALSE)
  else
    test
}

abiToCfx <- function (abiFile) {
  data <- read.csv(abiFile , skip =8)
  cols <- c("Well", "Sample.Name", "Target.Name", "Task", "Reporter", "Quantity")
  data_subset <- data[, cols]
  names(data_subset) <- c("Well", "Sample", "Target", "Content", "Fluor", "Starting Quantity (SQ)")
  data_subset$Cq <- data[, 7]
  data_subset$Cq[data_subset$Cq == "Undetermined"] <- "N/A"
  data_subset$Content[data_subset$Content == "UNKNOWN"] <- "Unkn"
  data_subset$Content[data_subset$Content == "STANDARD"] <- "Std"
  
  data_subset
}

getMeta <- function(file, format){
  if(format == "abi"){
    metadata <- read.csv(abiFile, nrow=5)[, 1]
    file <- tail(strsplit(metadata[2], "\\\\")[[1]], 1)
    endtime <- strsplit(metadata[3], "=")[[1]][2]
    c("File Name" = file,
      "Run Started" = endtime,
      "Protocol File Name" = NULL,
      "Sample Vol" = NULL)
  } 
  else if(format == "cfx"){
    meta <- read.csv(file, nrows = 13, header = FALSE)[, 1:2]
    metadata <- meta[, 2]
    names(metadata) <- meta[, 1]
    metadata <- gsub("_", " ", metadata)
    metadata
  }
}

# read in CFX formatted data
readAndClean <-  function(file, format){
  options(stringsAsFactors=FALSE)
  
  if(format == "cfx") {
    dat <- read.csv(file,
                    skip=19,
                    stringsAsFactors=FALSE)
  } else {
    dat <- abiToCfx(file)
  }
  
  
  # data Clean Up 
  
  names(dat)[names(dat) == "Starting.Quantity..SQ."] <- "CopyPeruL"
  dat$CopyPeruL <- as.numeric(dat$CopyPeru)
  dat$Cq[dat$Cq == "N/A"] <- m
  dat$Cq <- as.numeric(dat$Cq)
  dat$Target <- tolower(dat$Target)
  dat$Content[grepl("Std", dat$Content)] <- "Std"
  dat$Content[grepl("Unkn", dat$Content)] <- "Unkn"
  
  dat
} 

# PASS/FAIL standard curve
standardQC <- function(eff, r2, target, ef.max=eff.max, ef.min=eff.min, r2.m=r2.min) {
  effQC <- ifelse((eff > ef.min & eff < ef.max), "PASS", "FAIL")
  r2QC <- ifelse(r2 > r2.m, "PASS", "FAIL")
  data.frame(Target = target,
             Parameter = c("r-squared", "Amplification Factor"),
             Value = c(r2, eff),
             QC = c(r2QC, effQC))
}


controlFrame <- function (data, assay) {
  data$Parameter <- toupper(data$Sample)
  cData <- data[grepl("NTC|NEC", data$Sample),] 
  cData <- ddply(cData, .(Sample), function(df){
    df$Replicate <- paste0("Ct$_{Rep", 1:nrow(df), "}$")
    df
  })
  controlDF <- dcast(cData, Sample ~ Replicate, value.var="Cq")
  controlDF$QC <- apply(controlDF[, -1], 1, function(x)ifelse(all(x == m), "PASS", "FAIL"))
  controlDF[controlDF == m] <- "N/A"
  controlDF$Target <- assay
  controlDF[, c(ncol(controlDF), 1:(ncol(controlDF) - 1))]
}


sketaQC <- function(data, threshold=thres){
  sk.unkn <- data[grepl("Unkn", data$Content), ]
  calibrators <- data$Cq[grepl("NEC", toupper(data$Sample))]
  Ct.sk.calibrator <- mean(calibrators)
  
  sk.unkn$sk.dct <- sk.unkn$Cq - Ct.sk.calibrator
  sk.unkn$Inhibition <- ifelse(sk.unkn$sk.dct > threshold | sk.unkn$sk.dct < (-threshold),
                               "FAIL", "PASS")
  names(sk.unkn)[names(sk.unkn)=="Cq"] <- "sk.Ct"
  
  res <- ddply(sk.unkn, .(Sample), function(df){
    data.frame(Sample = unique(df$Sample),
               sk.Ct = mean(df$sk.Ct, na.rm=TRUE),
               sk.dct = mean(df$sk.dct, na.rm=TRUE),
               Inhibtion = ifelse(all(df$Inhibition == "PASS"), "PASS", "FAIL")
               
    )
  })
  names(res) <- c("Sample", "sketaCt$_{mean}$", "$\\Delta$Ct$_{mean}$", "QC")
  res
}