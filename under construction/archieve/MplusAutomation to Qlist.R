

library(lpa.mi.src)
library("MplusAutomation")
setwd("C:/Users/marcu/Dropbox/Dissertation/lpa-mi-src/under construction")
params_df = readModels(filefilter = "k-3naive dfcom p1 z1 rep1", what = "parameters")$parameters$unstandardized

Mplus2Qlist(params_df)