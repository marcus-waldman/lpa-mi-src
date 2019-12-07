#' Resolve label switching with population values known for a given file (row) value in a tracker.
#'
#' @param oneline_df (data.frame) The tracker
#' @param Plist (list) List with population values to compare against
#' @param out_get_FMM Optional. Default to NULL. If not provided, calculated using Plist with the Plist2getFMM() function.
#' @param args_KL_check_switching (list) Optional arguments pasted to KL_check_switching() beyond nn, mu_est_mat, S_est_array, z, data_conditoins
#' @param args_tracker2inpfile_naive (list) Optional arguments passed to tracker2inpfile_naive() beyond nn, tracker_df, pop_params_nn, data_conditions, Model_txt, process_cores
#' @param showOutput (logical) Defaluts to FALSE. This sets the showOutput argument in MplusAutomation::runModels()
#' @return Refits and generates new .out file correcting for switching. Returns a data.frame.
#' @export
#' @examples tracker2_resolve_label_switch<-function(nn,tracker_df, data_conditions,

tracker2_resolve_label_switch<-function(oneline_df,
                                        Plist,
                                        out_get_FMM = NULL,
                                        args_KL_check_switching = list(NULL),
                                        args_tracker2inpfile_naive = list(NULL),
                                        showOutput = FALSE){


      if(nrow(oneline_df)!=1){stop("oneline_df must be a single row")}

      require(lpa.mi.src)
      require(MplusAutomation)



      out_Mplus_nn = with(oneline_df, readModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                               filefilter = as.character(outfile),
                                               what = c("summaries"), quiet = TRUE));
      out_extract = extract_svals(file = as.character(paste0(oneline_df$outfile,".out")),
                                path = paste0(as.character(oneline_df$outwd),"/",as.character(oneline_df$outfolder)),
                                return_summaries = TRUE)
      svals_txt = out_extract$svals;
      LL_i = out_extract$LL

      model_txt = c("MODEL: \n", svals_txt)

      Qlist_n = syntax2Qlist(syntax_x = svals_txt, Plist_x = Plist)

      switched_df = data.frame(nn = oneline_df$nn, switched = NA, singular = FALSE, t(rep(NA,dim(Qlist_n$S)[3])))
      names(switched_df)[-c(1:3)] = paste0("c",seq(1,length(Qlist_n$pi)))
      switched_df = transform(switched_df,
                              L1_pi = NA,
                              L1_mu = NA,
                              L1_S = NA,
                              L1_LL = NA)



        # Covariance matrix is non singular
          args = args_KL_check_switching
          switched_df = KL_check_switching(nn = oneline_df$nn,
                                           Qlist_x = Qlist_n, Plist_x = Plist,
                                           name = ifelse(is.null(args$name),"glpk", args$name),
                                           t_max = ifelse(is.null(args$t_max), 30, args$tmax),
                                           approximate = ifelse(is.null(args$approximate), 0, args$approximate),
                                           round_cplex = ifelse(is.null(args$round_cplex), 0, args$round_cplex),
                                           trace_cplex = ifelse(is.null(args$trace_cplex), 0, args$trace_cplex))
          switched_df = transform(switched_df,
                                  L1_pi = NA,
                                  L1_mu = NA,
                                  L1_S = NA,
                                  L1_LL = NA)

          if (switched_df$switched==TRUE){
#stop(paste0("bookholder for nn=",oneline_df$nn))
            # Covariance matrix and label swithcing as occured
              K_z = dim(Plist$S)[3]

              map_to = as.numeric(switched_df[,startsWith(names(switched_df), "c")])

              # If class switching, refit the model
                # Create a properly switched model text
                  list_model_txt = list(NULL)
                  for(k in 1:K_z){
                    i_start = which(endsWith(model_txt, paste0("%C#",k,"%")))+1
                    i_end = ifelse(k==K_z,length(model_txt),which(endsWith(model_txt, paste0("%C#",k+1,"%"))))-1
                    k_to = map_to[k]
                    list_model_txt[[k_to]] = c(paste0("%C#",k_to,"%"), model_txt[i_start:i_end])
                  }
                  id_ref = which(map_to == K_z)
                  id_oth = match(seq(1,K_z-1), map_to)
                  gamma_vec = with(Qlist_n, log(pi[id_oth]/pi[id_ref]))
                  model_txt_to = c("MODEL:","\n",
                                   "%OVERALL%",paste0("[C#",seq(1,K_z-1),"*",round(gamma_vec,5),"];"), "",
                                   unlist(list_model_txt))

                # Fit the model specified by tne new model text
                  oneline_df =   oneline_df %>% transform(starts0 = 0)
                  args = args_tracker2inpfile_naive;
                  if(is.null(out_get_FMM)){out_get_FMM = Plist2getFMM(Plist)}
                  inp_txt = tracker2inpfile_naive(oneline_df = oneline_df,
                                                  pop_params_nn = out_get_FMM,
                                                  data_conditions = data_conditions,
                                                  savedata = ifelse(is.null(args$savedata), FALSE, args$savedata),
                                                  estimates = ifelse(is.null(args$estimates), FALSE, args$estimates),
                                                  results = ifelse(is.null(args$results), FALSE, args$results),
                                                  save_tech3 = ifelse(is.null(args$save_tech3), FALSE, args$save_tech3),
                                                  symbol_txt = ifelse(is.null(args$symbol_txt), "*",args$symbol_txt),
                                                  Model_txt = model_txt_to,
                                                  weight = ifelse(is.null(args$weight), FALSE, args$weight),
                                                  output_txt = args$output_txt,
                                                  processor_cores = 1)


                  with(oneline_df, runModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                             filefilter = as.character(paste0(outfile,".inp")),
                                             logFile = NULL,
                                             showOutput = showOutput))

                  # Compare with discrepancy metrics to see if the classes have moved after resolving label switching
                  out_Mplus_switched = with(oneline_df, readModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                                                     filefilter = as.character(outfile),
                                                                     what = c("summaries"), quiet = TRUE));
                  out_extract2 = extract_svals(file = as.character(paste0(oneline_df$outfile,".out")),
                                               path = paste0(as.character(oneline_df$outwd),"/",as.character(oneline_df$outfolder)),
                                               return_summaries = TRUE
                                              )
                  Qlist_model_to = syntax2Qlist(syntax_x = model_txt_to, Plist_x = Plist)
                  Qlist_switched = syntax2Qlist(syntax_x = out_extract2$svals, Plist_x = Plist)

                  switched_df = transform(switched_df,
                                          L1_pi = max(abs(Qlist_model_to$pi-Qlist_switched$pi)),
                                          L1_mu = max(abs(Qlist_model_to$mu-Qlist_switched$mu)),
                                          L1_S = max(abs(Qlist_model_to$S - Qlist_switched$S)),
                                          L1_LL = abs(out_extract2$LL-LL_i))


      } #end if (switched)


      return(switched_df)

}
