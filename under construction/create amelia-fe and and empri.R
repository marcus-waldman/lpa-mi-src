rm(list = ls())

library(plyr)
library(tidyverse)
library(lpa.mi.src)


path_wd = "C:/Users/marcu/Dropbox/Dissertation/saved/data-gen/test-mgAmelia/test-mgAmelia-BigB Jan 03 2019 13 44"

setwd(path_wd)
load(file = "environment.RData")
stage0_dir = "C:/Users/marcu/Dropbox/Dissertation/data-manipulation/lists-stage0"
dir_unzip = "C:/Users/marcu/Documents/zip-files/super-tmp" 

z = 37
head(data_conditions[z, ])

file_unzip = paste0("list-stage0 rep", rep, " z",z,".RData")
unzip(zipfile = paste(stage0_dir, stage0_zip, sep = "/"), 
      files = file_unzip, 
      exdir = dir_unzip)
load(paste0(dir_unzip,"/",file_unzip))


methods_list = list(procedure = list(NULL), name = list(NULL), args = list(NULL))
pva_vec = c(11,12,21,22)
methods_list$procedure[[1]] = "amelia-FE";           methods_list$name[[1]] = "Amelia-FE";    methods_list$args[[1]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4)
methods_list$procedure[[2]] = "amelia-FE";           methods_list$name[[2]] = "Amelia-FE";    methods_list$args[[2]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4, empri = 0.01)
methods_list$procedure[[3]] = "stratamelia";         methods_list$name[[3]] = "MG-Amelia";    methods_list$args[[3]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4)
methods_list$procedure[[4]] = "stratamelia";         methods_list$name[[4]] = "MG-Amelia";    methods_list$args[[4]] = list(m = M, p2s = as.numeric(VerboseImpute), tolerance = 1E-4, empri = 0.01)


methods_list$pva_vec = pva_vec


list_imputed = lpa.mi.src::get_imputed_data(z = z, 
                                            list_get_obs = list_observed,
                                            list_get_complete = list_complete,
                                            methods_list = methods_list,
                                            data_conditions = data_conditions,
                                            pctmiss_vec = pctmiss_vec,
                                            save_it = FALSE)

list_get_obs = list_observed
list_get_complete = list_complete
save_it = FALSE

