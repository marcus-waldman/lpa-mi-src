rm(list = ls())


remove.packages("lpa.mi.src")
.rs.restartR()

library(roxygen2)
library(devtools)

dropbox = "C:/Users/marcu/Dropbox"
#dropbox="D:/Dropbox/Dropbox"
#dropbox = "D:/Dropbox"
#dropbox = "Z:/Dropbox"

src_dir= paste(dropbox,"/Dissertation/lpa-mi-src", sep = "")
#create("lpa.mi.src")
setwd(src_dir)
setwd("lpa.mi.src")
document()

setwd(src_dir)
install("lpa.mi.src")

rm(list = ls())


#install_github("marcus-waldman/lpa-mi-src/lpa.mi.src")
#library("lpa.mi.src")
