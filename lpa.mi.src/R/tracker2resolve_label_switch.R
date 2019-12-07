#' Resolve label switching with population values known.
#'
#' @param nn (integer) Identifier in tracker_df
#' @param tracker_df (data.frame) The tracker
#' @param get_FMM_params (list) Optional. From get_FMM_params
#' @param args_KL_check_switching (list) Optional arguments pasted to KL_check_switching()
#' @param args_tracker2inpfile_naive (list) Optional arguments passed to tracker2inpfile_naive()
#' @return Refits and generates new .out file correcting for switching. Returns a data.frame.
#' @export
#' @examples tracker2resolve_label_switch<-function(nn,tracker_df, data_conditions,

tracker2resolve_label_switch<-function(nn,tracker_df, data_conditions,
                                       out_get_FMM=NULL,
                                       args_KL_check_switching = list(NULL),
                                       args_tracker2inpfile_naive = list(NULL)){

      require(MplusAutomation)
      if(is.null(out_get_FMM)){out_get_FMM<-get_FMM_params(tracker_df$z[nn],data_conditions)}


      out_Mplus_nn = with(tracker_df[nn, ], readModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                               filefilter = as.character(outfile),
                                               what = c("parameters"), quiet = TRUE));
      Qlist_nn = Mplus2Qlist(out_Mplus_nn$parameters$unstandardized)
      switched_df = with(Qlist_nn, KL_check_switching(nn = nn, mu_est_mat = mu, S_est_array = S, z = tracker_df$z[nn], data_conditions = data_conditions))
      if (switched_df$switched==TRUE){
          K_z = out_get_FMM$K_z
          map_to = as.numeric(switched_df[,startsWith(names(switched_df), "c")])

          # If class switching, refit the model
            # Extract model text
            svals_txt = extract_svals(file = as.character(paste0(tracker_df$outfile[nn],".out")),
                                      path = paste0(as.character(tracker_df$outwd[nn]),"/",as.character(tracker_df$outfolder[nn]))
                        )
            model_txt = c("MODEL: \n", svals_txt)

            # Create a properly switched model text
              list_model_txt = list(NULL)
              for(k in 1:K_z){
                i_start = which(endsWith(model_txt, paste0("%C#",k,"%")))+1
                i_end = ifelse(k==K_z,length(model_txt),which(endsWith(model_txt, paste0("%C#",k+1,"%"))))-1
                k_to = map_to[k]
                list_model_txt[[k_to]] = c(paste0("%C#",k_to,"%"), model_txt[i_start:i_end])
              }
              gamma_vec = with(Qlist_nn, log(pi[-K_z]/pi[K_z]))
              model_txt_to = c("MODEL:","\n",
                               "%OVERALL%",paste0("[C#",seq(1,K_z-1),"*",round(gamma_vec,5),"];"), "",
                               unlist(list_model_txt))

            # Fit the model specified by tne new model text
              tracker_nn =   tracker_df[nn, ] %>% transform(starts0 = 0)
              inp_txt = do.call(what = "tracker2inpfile_naive",
                                args = list(nn = 1, tracker_df = tracker_nn, pop_params_nn = out_get_FMM,
                                            data_conditions = data_conditions, Model_txt = model_txt_to,
                                            processor_cores = 1))
              with(tracker_nn[1, ], runModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                              filefilter = as.character(paste0(outfile,".inp")),
                                              logFile = NULL))

        }
      return(switched_df)

}
