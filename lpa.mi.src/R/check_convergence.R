#' Check convergence
#'
#' Scans an .out file from an Mplus model and looks for signs of convergence.
#' @param file (integer)
#' @param folder (character)
#' @param starts_txt (character)
#' @param type_imputation (logical)
#' @return (logical)
#' @export
#' @examples


check_convergence<-function(file, folder_wd = getwd(), starts_txt = "0;", type_imputation = FALSE, Rcond_min = 1E-6){

  require(stringr)

#folder_wd = "C:/Users/marcu/Dropbox/Dissertation Proposal/Paper 1/Simulation Code/Impute LPA Sim - known K -  v1_0/Error diagnosis/p1/Complete data"
#file = "Naive dfcom p1 z1 rep15 v2.out"

    #setwd(folder_wd)

    hi = readLines(con = paste0(folder_wd, "/",file))
    i_0 = which(hi ==  "THE MODEL ESTIMATION TERMINATED NORMALLY" )
    i_1 = which(hi ==  "THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED.  RERUN WITH AT LEAST TWICE THE" )

    Terminate_Normal = FALSE

    if(type_imputation == FALSE){

        if ("Number of replications" %in% hi){
          warning("Type=imputation; appears to be specified according to .out file, but type_imputation = FALSE is set in check_convergence.")
        }

        if (length(i_0)>0){
          if (starts_txt=="0;"){Terminate_Normal = TRUE}
          if (starts_txt!="0;" & length(i_1)==1){Terminate_Normal = TRUE}
        }
    }

    if(type_imputation == TRUE){
      if (starts_txt == "0;"){ Terminate_Normal = TRUE}
      if (starts_txt != "0;" & length(i_1)==1){Terminate_Normal = TRUE}
    }

    Rcond = -9
    j_0 = which(hi == "QUALITY OF NUMERICAL RESULTS")
    tmp = NA
    if(length(j_0)>0){
      tmp = stringr::str_extract_all(hi[j_0+2],"[0-9]+")[[1]]
      Rcond = as.numeric(paste0(tmp[1],".",tmp[2],"E-",tmp[3]))
    }


    poor_Rcond = Rcond<=Rcond_min

    normal = Terminate_Normal


    return(data.frame(normal=normal,
                      Rcond = ifelse(Rcond==-9, NA, Rcond),
                      poor_Rcond = poor_Rcond))
}
