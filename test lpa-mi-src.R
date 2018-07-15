
library(plyr)
library(tidyverse)
library(lpa.mi.src)

Processors = 2
Replications = 15

dropbox_wd =  "C:/Users/marcu/Dropbox"  #  Dropbox folder on local machine
pool_wd = paste(dropbox_wd, "/Dissertation/lpa-mi-pool", sep = "")

setwd(pool_wd)
load("image lpa-mi-pool.RData")

