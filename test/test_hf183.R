library(knitr)
library(xtable)

source("r/hf183.r")
sketasource <- list(datapath = "www/hfsketafile.csv")
hfresults <- hf("www/hfdatafile.csv", sketasource, "cfx", "test")
results <- function()hfresults

knit2pdf(input="templates/hf/report.Rtex", output = "templates/hf/report.tex",
         compiler="xelatex")
