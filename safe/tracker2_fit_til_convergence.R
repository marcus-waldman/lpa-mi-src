#' Fit to convergence using a tracker
#'
#' @param data_conditions (data.frame)
#' @param tracker_df (data.frame)
#' @param Processors (integer) Defaults to parallel::detectCores()-1
#' @param itermax (numeric) Default is 4^3 (or 64). Specifying the number of minutes a model should take to fit (on average) before terminating. Defaults to 1 minute.
#' @param multiplier (numeric) Multiplier across each iteration by which random starts should increase
#' @param Rcond_min (numeric) Threshold below which results are flagged for poor conditioning number.
#' @param cl cluster object from makeCluster
#'
#' @return tracker_df
#' @export
#'
#' @examples tracker2_fit_til_convergence(data_conditions, tracker_df)
#'
tracker2_fit_til_convergence<-function(data_conditions, tracker_df, itermax = 2, multiplier = 2^5, Rcond_min=1E-6, cl = NULL){

    require(parallel)
    require(pbapply)
    require(MplusAutomation)
    require(lpa.mi.src)

    # Get the number of cores on the machine
    numCores = parallel::detectCores()
    processor_cores = 1

    # Get the population parameters for each data condition
    list_pop_params = lapply(sort(unique(tracker_df$z)), FUN = function(x){get_FMM_params(z=x,data_condition=data_conditions)})
    names(list_pop_params) = paste0("z",sort(unique(tracker_df$z)))


    tracker_df = transform(tracker_df, normal = NA, Rcond = NA, poor_Rcond = NA)
    nn_vec = 1:nrow(tracker_df); flag = 0; iter = 1;
    while(flag == 0){

      if (iter==1){nn_0 = length(nn_vec)}
      if (is.null(cl) | nn_0 < numCores){
        processor_cores = parallel::detectCores()
      }


      print(" "); print(" ")
      print(paste0("Iteration ", iter, ": Fitting ", length(nn_vec)," Mplus Models"))

      # Create directories
      print("Writing .inp files...")
      tmp = pblapply(X = nn_vec, FUN = function(x){require(lpa.mi.src);
        idz = with(tracker_df[x, ], which(names(list_pop_params)==paste0("z",z)));
        pop_params_nn = list_pop_params[[idz]];
        tracker2inpfile_naive(nn=x,
                              tracker_df = tracker_df,
                              data_conditions = data_conditions,
                              processor_cores = processor_cores,
                              pop_params_nn = pop_params_nn)
        },
        cl = NULL
        )
      rm(tmp)


      #x = 1 ;
      #with(tracker_df, runModels(target = as.character(outfolder[x]), filefilter = as.character(outfile[x])))

      # Run the files
      print("Fitting LPA models in Mplus...")


      if (is.null(cl)){
          tmp = pblapply(X = nn_vec,
                        FUN = function(x){require(MplusAutomation); require(lpa.mi.src);
                                          with(tracker_df[x, ], runModels(target = paste0(outwd,"/",as.character(outfolder)),
                                               filefilter = as.character(paste0(outfile,".inp")),
                                               logFile = NULL));
                             },
                        cl = NULL)
      } else {
             print(cl)
             pb <- pbapply::timerProgressBar(max = length(nn_vec), style = 1, width = getOption("width")/4)
             progress <- function(x){setTimerProgressBar(pb, x)}
             opts <- list(progress = progress)
             tmp<-foreach(n = nn_vec, .packages= ("MplusAutomation"), .inorder = TRUE, .options.snow = opts) %dopar% {
                  with(tracker_df[n, ], runModels(target = paste0(outwd,"/",as.character(outfolder)),
                                                  filefilter = as.character(paste0(outfile,".inp")),
                                                  logFile = NULL));
                    return(NULL)
             }
             close(pb)
      }
      rm(tmp)


      # Check convergence
      print("Assessing convergence...")
      convergence_list = pblapply(X = nn_vec,
                                  FUN = function(x){require(lpa.mi.src);
                                    convergence_df = with(tracker_df[x, ], check_convergence(file = as.character(paste0(outfile,".out")),
                                                                                             folder_wd =  paste0(as.character(outwd),"/",as.character(outfolder)),
                                                                                             starts_txt = ifelse(starts0>0, paste0(starts0, " ", starts0/4, ";"), "0;"),
                                                                                             Rcond_min = Rcond_min));
                                    return(convergence_df)
                                  },
                                  cl = NULL)
      convergence_df = rbindlist(convergence_list)
      tracker_df$normal[nn_vec] = convergence_df$normal
      tracker_df$Rcond[nn_vec] = convergence_df$Rcond
      tracker_df$poor_Rcond[nn_vec] = convergence_df$poor_Rcond

      nn_vec = with(tracker_df, which(normal == FALSE | poor_Rcond == TRUE))
      nn_iter = length(nn_vec)

      iter = iter+1
      if ( nn_iter == 0 | (iter>itermax)){
        flag = 1
        print(" ")
        print(" ")
        print(paste0(">> Terminating with ",nn_iter, " (",100*round(nn_iter/nn_0,4),"%) model(s) having failed to converge."))
        print(paste0("...", sum(tracker_df$normal==FALSE), " continue to display abnormal termination."))
        print(paste0("...", sum(tracker_df$poor_Rcond==TRUE), " continue to display a condition number less than minimum of ", Rcond_min,"."))
      } else {
        print(paste0(nn_iter, " (",100*round(nn_iter/nn_0,4),"%) models have yet to converge."))
        print(paste0("...", sum(tracker_df$normal==FALSE), " displayed abnormal termination."))
        print(paste0("...", sum(tracker_df$poor_Rcond==TRUE), " displayed a condition number less than minimum of ", Rcond_min,"."))
        tracker_df$starts0[nn_vec] = multiplier*tracker_df$starts0[nn_vec]
      }



    } #end while

    return(tracker_df)
}
