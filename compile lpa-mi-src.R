rm(list = ls())

library(roxygen2)
library(devtools)

dropbox = "C:/Users/marcu/Dropbox"

src_dir= paste(dropbox,"/Dissertation/lpa-mi-src", sep = "")

setwd(src_dir)

create("lpa.mi.src")
document(setwd)
