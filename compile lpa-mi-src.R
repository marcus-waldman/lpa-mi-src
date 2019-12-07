rm(list = ls())
remove.packages("lpa.mi.src")
.rs.restartR()


library(roxygen2)
library(devtools)

# #Acer-Laptop
  #dropbox = "C:/Users/marcu/Dropbox"
# Big Bertha
  #dropbox="D:/Dropbox/Dropbox"
# Gaby/QMR/Mar-Cobra
  dropbox = "D:/Dropbox"
# IQSS-RCE NoMachine
  #dropbox = "/nfs/home/M/mwaldman/diss-nomx/Dropbox"
  

src_dir= paste(dropbox,"/Dissertation/lpa-mi-src", sep = "")
#create("lpa.mi.src")
setwd(src_dir)
setwd("lpa.mi.src")
document()

setwd(src_dir)
install("lpa.mi.src", upgrade = FALSE)

rm(list = ls())


#install_github("marcus-waldman/lpa-mi-src/lpa.mi.src")
#library("lpa.mi.src")
