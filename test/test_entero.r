library(knitr)
library(xtable)

source("r/entero.r")

entresults <- ent("www/ent_example.csv", NULL, "cfx", "test")
results <- function()entresults

knit2pdf(input="templates/ent/report.Rtex", output = "templates/ent/report.tex",
         compiler="xelatex")


