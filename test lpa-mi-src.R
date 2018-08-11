rm(list = ls())
library(plyr)
library(tidyverse)
library(lpa.mi.src)
library(MplusAutomation)


dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
pool_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-pool", sep = "")

setwd(pool_wd)
load("image lpa-mi-pool.RData")

z=1; p=1; rep = 1;


out_comp = get_complete_data(z = z,data_conditions = data_conditions, rep = rep, p = p, save_it = FALSE, temp_wd_p = temp_wd_p_vec)  




out_obs = get_obs_data(z = z, df = out_comp$dfcom, pctmiss_vec = pctmiss_vec, data_conditions = data_conditions, save_it = TRUE, temp_wd_p_vec = temp_wd_p_vec, rep = rep, p = p)

out_imp = get_imputed_data(z=z,list_get_obs = out_obs, list_get_complete = out_comp, methods_list = methods_list, data_conditions = data_conditions, pctmiss_vec = pctmiss_vec,  save_it = TRUE, temp_wd_p_vec = temp_wd_p_vec, rep = rep, p = p)


out_get_FMM = get_FMM_params(z, data_conditions)

inp_Mplus_comp <- create_naiveMplus_inpfile(out_get_FMM = out_get_FMM, dffolderfiles = out_comp$dffolderfiles, temp_wd_p_vec = temp_wd_p_vec, type_imputation = FALSE, starts_txt = "0;", symbol_txt = "*")
inp_Mplus_imp <- create_naiveMplus_inpfile(out_get_FMM = out_get_FMM, dffolderfiles = out_imp$dffolderfiles, temp_wd_p_vec = temp_wd_p_vec, type_imputation = TRUE, starts_txt = "0;", symbol_txt = "*")
