rm(list = ls())
library(plyr)
library(tidyverse)
library(lpa.mi.src)

Processors = 2
Replications = 15

dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
pool_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-pool", sep = "")

setwd(pool_wd)
load("image lpa-mi-pool.RData")

z=2; p=1; rep = 1;
out_comp = get_complete_data(z = z,data_conditions = data_conditions, save_it = FALSE)  
out_obs = get_obs_data(z = z, df = out_comp$dfcom, pctmiss_vec = pctmiss_vec, data_conditions = data_conditions, save_it = FALSE)
out_imp = get_imputed_data(z=z,list_get_obs = out_obs, list_get_complete = out_comp, methods_list = methods_list, data_conditions = data_conditions, pctmiss_vec = pctmiss_vec,  save_it = FALSE)

