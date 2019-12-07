#' Create Mplus .inp file from a tracker
#'
#' Creates an input file that fits a simple (i.e., naive) LPA model given information from a tracker data frame.
#' @param oneline_df (data.frame) One row from the tracker data frame for monitoring covernce, specifying dat file, specifying .inp file path, etc.
#' @param pop_params_nn (list) Optional. Gaussian finite mixture model parameters in format exported by get_FMM_params for the nn-th condition in the tracker. Will obtain if not specified.
#' @param data_conditions (data.frame) Optional, but must be provided if pop_params_nn not given.
#' @param savedata (logical) Save the cpros or not (savedata = FALSE). Defaults to false.
#' @param estimates (logical) Save the parameter estimates (i.e. include a ESTIMATES = ... statement in the Savedata Mplus command). Defaults to FALSE.
#' @param results (logical) Save the results (i.e. include a RESULTS = ... statement in the Savedata Mplus command). Defaults to FALSE.
#' @param save_tech3 (logical) Save the asympototic variance covariance matrix. Defaults to FALSE.
#' @param symbol_txt ("character") Either "*" or the "(at)" symbol. The latter for fixing parameters. Defaults to star.
#' @param Model_txt (character vector). Syntax for inside MODEL: command (exclude a "MODEL:" line when specifying). Defaults to NULL.
#' @param weight (logical) Indicates if a sample weight is included as the last variable in the .dat files listed in dffolderfiles. Defaults to FALSE.
#' @param output_txt (character) Optional. Text for Output block to be requested in addition to the default 'tech1;', 'svalues;', 'tech3;'
#' @param processor_cores (integer) Number of processor cores to use. Defaults to 1.
#' @param do_dk (logical) Tells me whether or not fit a model with smaller number of components
#' @return Saves the Mplus input files for Naive LPA analysis in the the folders specified of dffolderfiles input. Returns a character vector with the the Mplus code.
#' @export
#' @examples
#' tracker2inpfile_naive(oneline_df)

tracker2inpfile_naive<-function(oneline_df, pop_params_nn = NULL, data_conditions = NULL, savedata = FALSE, estimates = FALSE, results = FALSE, save_tech3 = FALSE,
                                symbol_txt = "*", Model_txt = NULL, weight = FALSE, output_txt = NULL, processor_cores = 1, do_dk = FALSE){

    if(is.null(pop_params_nn)){
      if(!is.null(data_conditions)){
        pop_params_nn = with(oneline_df, lpa.mi.src::get_FMM_params(z, data_conditions = data_conditions))
      } else {
        stop("data_conditions must be specified since pop_params_nn not specified")
      }
    }

    J = pop_params_nn$J
    J_Xcom_z = pop_params_nn$J_Xcom_z
    J_Xinc_z = pop_params_nn$J_Xinc_z
    J_Y_z = pop_params_nn$J_Y_z
    K_z = pop_params_nn$K_z
    if (do_dk==TRUE){K_z=oneline_df$kfit}
    if (is.null(pop_params_nn$mu_z)){
      full_mu_z = matrix(runif(J*K_z, min = 0, max = pop_params_nn$MD_z+1.96), nrow = J, nc = K_z)
    } else {
      full_mu_z = round(pop_params_nn$mu_z,5)
    }
    if (is.null(pop_params_nn$su_z)){
      full_S_z = array(data = 0, dim = c(J,J,K_z))
      for (kk in 1:K_z){
        full_S_z[,,kk] = diag(J)
      }
    } else {
      full_S_z = round(pop_params_nn$S_z,5)
    }
    if (is.null(pop_params_nn$pi_z)){
      pi_z = rep(1/K_z,K_z)
    } else {
      pi_z = as.numeric(round(pop_params_nn$pi_z,5))
    }


    #### Pre-processing ####
    starts0_nn = oneline_df$starts0
    if (starts0_nn == 0){
      starts_txt = "0;"
    } else {
      starts_txt = with(oneline_df,paste0(starts0," ",floor(starts0/4), ";"))
    }
    if(!endsWith(starts_txt, ";")){stop("starts_txt must end with a semicolon")}

    #### Title Command ####
    Title_txt = with(oneline_df,
                     paste0("TITLE: Code for Naive LPA model\n",
                            outfile)
    )

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
                      c("MISSING = .;"))
    if (K_z>1){Variables_txt = c(Variables_txt, paste("CLASSES = ", classes_txt, sep = ""))}
    if (weight) {Variables_txt = c(Variables_txt, "WEIGHT = wgt;")}


    #### Analysis Command ####
    Estimator_txt = ifelse(weight, "Estimator = MLR;", "Estimator = ML;")
    if(K_z>1){
      Analysis_txt = c("ANALYSIS: ",
                       "TYPE = MIXTURE;",
                       Estimator_txt,
                       "INFORMATION = obs;",
                       paste("STARTS = ", starts_txt, sep = ""),
                       paste0("PROCESSORS = ", processor_cores,";"))
    } else {
      Analysis_txt = c("ANALYSIS: ",
                       "TYPE = GENERAL;",
                       Estimator_txt,
                       "INFORMATION = obs;",
                       paste0("PROCESSORS = ", processor_cores,";"))
    }


    #### Model Command ####
    if (is.null(Model_txt) & K_z>1){
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

    if (is.null(Model_txt) & K_z==1){
      withinclass_txt = naive_withinclass_Mplus_syntax(mu_k = full_mu_z[1:J_Y_z,1],
                                                       S_k = full_S_z[1:J_Y_z,1:J_Y_z, 1],
                                                       usev_vec = usev_vec,
                                                       Model_Covariance = ifelse(K_z==1,TRUE,FALSE))
      Model_txt = c("MODEL:", withinclass_txt)
    }

    #### Make and save the input files ####


      # Construct the Data and Savedata commands
      datwd_nn = with(oneline_df, gsub("/","/\n",as.character(datwd)))
      datfolder_nn = with(oneline_df, gsub("/","/\n",as.character(datfolder)))
      datfile_nn = with(oneline_df, as.character(datfile))
      Data_txt = c("DATA:",
                   paste0("FILE = '", datwd_nn,"/\n",datfolder_nn,"/\n",datfile_nn,"';"))



      Savedata_txt = c("SAVEDATA:")
      if(savedata==TRUE & K_z>1){
        Savedata_txt = c(Savedata_txt,
                         paste("FILE = cprob ", datfile_nn, ";",sep = ""),
                         "SAVE = cprob;")
      }
      if (estimates==TRUE){
        Savedata_txt = c(Savedata_txt,
                         paste("ESTIMATES = est-", datfile_nn, ";",sep = ""))
      }
      if (results==TRUE){
        Savedata_txt = c(Savedata_txt,
                         paste0("RESULTS = results-", datfile_nn, ";"))
      }
      if (save_tech3==TRUE){
        Savedata_txt = c(Savedata_txt,
                         paste0("TECH3 = tech3-", datfile_nn, ";"))
      }
      if (length(Savedata_txt)>1){Savedata_txt = c(Savedata_txt, "FORMAT = F10.6;") }

      # Create the output
      Output_txt = c("OUTPUT:",
                     "tech1;",
                     "tech3;",
                     "svalues;",
                     output_txt)

      # Put all the Mplus commands together into a single text vector
      mplus_txt = c(Title_txt, "\n",
                    Data_txt, "\n",
                    Variables_txt, "\n",
                    Analysis_txt, "\n",
                    Model_txt, "\n",
                    Output_txt, "\n",
                    Savedata_txt)

      # Change directory and save the Mplus input file
      fileConn<-file(with(oneline_df,paste0(outwd,"/",outfolder,"/",outfile,".inp")))
      writeLines(mplus_txt, fileConn)
      close(fileConn)

      return(mplus_txt)

}
