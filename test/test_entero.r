library(knitr)
library(xtable)

source("r/entero.r")

results <- entero("data/enttest.csv", "cfx", "test")

knit2pdf(input="templates/ent/report.Rtex", output = "templates/ent/report.tex",
         compiler="xelatex")


