rm(list = ls())



library(lpa.mi.src)

setwd("D:/Dropbox/Dissertation/lpa-mi-src/under construction/approx-bayes-proper")
load("environment.RData")
load(paste0("list-observed rep",91, " z",1," pm", 1, ".RData"))
load(paste0("list-complete rep",91, " z",1, ".RData"))


t00<-proc.time()
obj_call = EM_with_sampling_FMM(obsdf = list_observed$list_obsdf[[1]], M = 100, z = 1, data_conditions = data_conditions, 
                                tempdir_mplus = c("C:/Users/marcu/Documents"), rep = 91)
proc.time()-t00
