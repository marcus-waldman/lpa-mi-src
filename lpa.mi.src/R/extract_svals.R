#' Extract starting values from Mplus .out file.
#'
#' Extract starting values from Mplus .out file.
#' @param file (character) The filename of the Mplus .out file where the starting values are located
#' @param path (character) The path where the file is located
#' @param return_summaries (logical) If request a list with svals and summary statistics is parsed.
#' @return (character) vector with Model sytax at starting values from the .out file.
#' @export
#' @examples
#' extract_svalues(file = "example.out")

extract_svals <- function(file, path=getwd(), return_summaries = FALSE){

  # Input:
  #   A) file - The filename of the Mplus .out file where the starting values are located
  #   B) path - The path where the file is located
  # Output:
  #   hi - Character vector with starting values.

  hi = readLines(con = paste(path,"/",file, sep = ""))

  i_0 = which(hi == "MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES")
  if (length(i_0)==0){
      stop(paste("svalues in ", file,  " not found.", sep = ""))
    }
  flag = 0
  i = 0
  while(flag <2){
    i = i+1
    if ((hi[i_0+ i]) == ""){
      flag = flag+1
    } else {
      flag = 0
    }

  }
  i_f = i_0+i

  svals = hi[(i_0+1):i_f]

  if (return_summaries==FALSE){
    return(svals)
  } else{
    out_list = list(svals = svals)

    require(stringr)
    require(tidyverse)
    i_0 = which(str_detect(hi, "H0 Value"))
    LL = hi[i_0] %>% str_remove("H0 Value") %>% str_squish() %>% as.numeric()
    out_list$LL = LL

    return(out_list)
  }

}
