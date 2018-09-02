#' Create Mplus .inp file
#'
#' Creates an input file that fits a simple (i.e., naive) LPA model
#' @param z (integer) Condition number
#' @param out_get_FMM (list) Gaussian finite mixture model parameters in format exported by get_FMM_params
#' @param dffolderfiles  (data.frame) with the files and folders of the saved data
#' @param temp_wd_p (character). Specifies processor-specific temporary directory
#' @param type_imputation (logical) Add a "TYPE = IMPUTATION;" command in the Data command of Mplus. Defaults to FALSE (i.e. not TYPE statement in DATA command)
#' @param savedata (logical) Save the cpros or not (savedata = FALSE). Defaults to false.
#' @param estimates (logical) Save the parameter estimates (i.e. include a ESTIMATES = ... statement in the Savedata Mplus command). Defaults to FALSE.
#' @param results (logical) Save the results (i.e. include a RESULTS = ... statement in the Savedata Mplus command). Defaults to FALSE.
#' @param save_tech3 (logical) Save the asympototic variance covariance matrix. Defaults to FALSE.
#' @param starts_txt (character) Custom starting values in Mplus (e.g "20 10;" correpsonds to "STARTS = 20 10;"). Defaults to "0;"
#' @param symbol_txt ("character") Either "*" or the "(at)" symbol. The latter for fixing parameters.
#' @param Model_txt (character vector). Defaults to NULL.
#' @param weight (logical) Indicates if a sample weight is included as the last variable in the .dat files listed in dffolderfiles.
#' @return Saves the Mplus input files for Naive LPA analysis in the the folders specified of dffolderfiles input. Returns a character vector with the the Mplus code.
#' @export
#' @examples
#' create_naiveMplus_inpfile(z, data_conditions, dffolderfiles, temp_wd_p)
#' 
create_naiveMplus_inpfile<-function(z, out_get_FMM, dffolderfiles, temp_wd_p, savedata = FALSE, estimates = FALSE, results = FALSE, save_tech3 = FALSE, type_imputation = FALSE, starts_txt = "0;", symbol_txt = "*", Model_txt = NULL, weight = FALSE){

  #
  J = out_get_FMM$J
  J_Xcom_z = out_get_FMM$J_Xcom_z
  J_Xinc_z = out_get_FMM$J_Xinc_z
  J_Y_z = out_get_FMM$J_Y_z
  K_z = out_get_FMM$K_z
  if (is.null(out_get_FMM$mu_z)){
    full_mu_z = matrix(runif(J*K_z, min = 0, max = out_get_FMM$MD_z+1.96), nrow = J, nc = K_z)
  } else {
    full_mu_z = round(out_get_FMM$mu_z,5)
  }
  if (is.null(out_get_FMM$su_z)){
    full_S_z = array(data = 0, dim = c(J,J,K_z))
    for (kk in 1:K_z){
      full_S_z[,,kk] = diag(J)
    }
  } else {
    full_S_z = round(out_get_FMM$S_z,5)
  }
  if (is.null(out_get_FMM$pi_z)){
    pi_z = rep(1/K_z,K_z)
  } else {
    pi_z = as.numeric(round(out_get_FMM$pi_z,5))
  }

    # Output: Constructs the Mplus input files for Naive LPA analysis

#Error diagnostic stuff
# dffolderfiles = out_comp$dffolderfiles;
# type_imputation = FALSE;
# starts_txt = "0;"
# symbol_txt = "*";


    #### Pre-processing ####
    if(!endsWith(starts_txt, ";")){stop("starts_txt must end with a semicolon")}

    #### Title Command ####
        Title_txt = paste("TITLE: Code for Naive FMM model for data condition z=", z, sep = "")

    #### Variable Command ####

        # Names:
        names_Xinc = ifelse(J_Xinc_z>0,  paste("Xinc", seq(1,J_Xinc_z), sep = ""), "")
        names_Xcom = ifelse(J_Xcom_z>0,  paste("Xcom", seq(1,J_Xcom_z), sep = ""), "")
        names_vec = c(paste("Y",seq(1,J_Y_z), sep = ""),
                      names_Xinc,
                      names_Xcom, 
                      "subpop")
        names_txt = names_vec[1]
        for(t in seq(2,length(names_vec))){
          names_txt = paste(names_txt, names_vec[t], sep = " ")
        }
        if (weight){
          names_txt = paste(names_txt, "wgt", sep = " ")
        }
        names_txt = paste(names_txt, ";", sep = "")

        # Usevariables:
        usev_vec = names_vec[startsWith(names_vec,"Y")]
        usev_txt = usev_vec[1];
        for(t in seq(2, length(usev_vec))){
          usev_txt = paste(usev_txt, usev_vec[t], sep = " ")
        }
        usev_txt = paste(usev_txt, ";", sep = "")

        # Classes
        classes_txt = paste("c(", K_z,");", sep = "")

        # Put it all together
        Variables_txt = c("VARIABLE:",
                          paste("NAMES = ", names_txt, sep = ""),
                          paste("USEV = ", usev_txt, sep = ""),
                          "AUXILIARY = subpop;",
                          paste("CLASSES = ", classes_txt, sep = ""),
                          c("MISSING = .;"))
        if (weight) {Variables_txt = c(Variables_txt, "WEIGHT = wgt;")}


    #### Analysis Command ####
            Estimator_txt = ifelse(weight, "Estimator = MLR;", "Estimator = ML;")
            Analysis_txt = c("ANALYSIS: ",
                             "TYPE = MIXTURE;",
                             Estimator_txt,
                             paste("STARTS = ", starts_txt, sep = ""))


    #### Model Command ####
        if (is.null(Model_txt)){
            # write the overall Mplus command syntax
            marginals_txt = NULL
            if (K_z>1){
              marginals_txt = paste("[c#1*", "];", sep = "")
              if(K_z>2){
                for(k in seq(2,K_z-1)){
                  marginals_txt = c(marginals_txt,
                                   paste("[c#",k,symbol_txt, "];", sep = ""))
                }
              }
            }
            overall_txt = c("%OVERALL%",
                            marginals_txt,
                            usev_txt)

            # write Mplus syntax for the class-specific models
            classesmodels_txt = NULL
            for (k in c(1:K_z)){

              withinclass_txt = naive_withinclass_Mplus_syntax(mu_k = full_mu_z[1:J_Y_z,k],
                                                         S_k = full_S_z[1:J_Y_z,1:J_Y_z, k],
                                                         usev_vec = usev_vec,
                                                         Model_Covariance = ifelse(K_z==1,TRUE,FALSE))
              classesmodels_txt = c(classesmodels_txt,
                                    paste("\n %c#",k,"%",sep = ""),
                                    withinclass_txt)
            }

         # Put the overall and class-specific models together
         Model_txt = c("MODEL:",overall_txt,classesmodels_txt)
        } # end If model_txt

    #### Make and save the input files ####

        # For each simulated data file in the data condition z, construct and save an Mplus input file
        dffolderfiles_z = subset(dffolderfiles, data_condition==z)
        Htot = nrow(dffolderfiles_z)
        for(h in seq(1,Htot)){

          # Construct the Data and Savedata commands
          Data_txt = c("DATA:",
                         paste("FILE = ", dffolderfiles_z$files[h], ";",sep = ""))
          if(type_imputation){
            Data_txt = c(Data_txt, "TYPE = Imputation;")
          }

          Savedata_txt = c("SAVEDATA:")
          if(savedata==TRUE){
            Savedata_txt = c(Savedata_txt,
                               paste("FILE = cprob ", dffolderfiles_z$files[h], ";",sep = ""),
                               "SAVE = cprob;")
          }
          if (estimates==TRUE){
            Savedata_txt = c(Savedata_txt,
                             paste("ESTIMATES = est-", dffolderfiles_z$files[h], ";",sep = ""))
          }
          if (results==TRUE){
            Savedata_txt = c(Savedata_txt,
                             paste0("RESULTS = results-", dffolderfiles_z$files[h], ";"))
          }
          if (save_tech3==TRUE){
            Savedata_txt = c(Savedata_txt,
                             paste0("TECH3 = tech3-", dffolderfiles_z$files[h], ";"))
          }
          if (length(Savedata_txt)>1){Savedata_txt = c(Savedata_txt, "FORMAT = F10.6;") }

          # Create the output
          Output_txt = c("OUTPUT:",
                         "tech1;",
                         "svalues;")

          # Put all the Mplus commands together into a single text vector
          mplus_txt = c(Title_txt, "\n",
                        Data_txt, "\n",
                        Variables_txt, "\n",
                        Analysis_txt, "\n",
                        Model_txt, "\n",
                        Output_txt, "\n",
                        Savedata_txt)

          # Change directory and save the Mplus input file
          setwd(paste(temp_wd_p,"/", dffolderfiles_z$folders[h], sep = ""))
          fileConn<-file(paste("k-",K_z,"Naive ", sub(pattern = ".dat", ".inp",x = dffolderfiles_z$files[h]), sep = ""))
          writeLines(mplus_txt, fileConn)
          close(fileConn)

        }

        return(mplus_txt)

}

