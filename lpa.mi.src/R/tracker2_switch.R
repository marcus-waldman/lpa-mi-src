#' Address label switching from a tracker of .out files
#'
#' @param tracker_df (data.frame) The tracker
#' @param data_conditions (data.frame) data conditions
#' @param cl optional cluster object
#' @param args_KL_check_switching (list) Additional arguments to lpa.mi.src::KL_check_switching
#' @param args_tracker2inpfile_naive (list) additional argument to lpa.mi.src::tracker2inpfile_naive()
#'
#' @return
#' @export
#'
#' @examples
tracker2_switch<-function(tracker_df, data_conditions, cl = NULL,
                                   args_KL_check_switching = NULL, args_tracker2inpfile_naive = NULL){

  require(psych)
  require(parallel)
  require(MplusAutomation)
  require(lpa.mi.src)
  require(pbapply)

  if(!("converged" %in% names(tracker_df))){stop("converged must be a variable in tracker_df")}

  nn_vec = with(tracker_df, nn[converged == TRUE])
  # Get the population parameters for each data condition
  list_pop_params = lapply(sort(unique(tracker_df$z)), FUN = function(x){get_FMM_params(z=x,data_condition=data_conditions)})
  names(list_pop_params) = paste0("z",sort(unique(tracker_df$z)))


  # diagnose_list = NULL; iter = 0
  #  for (x in nn_vec){
  #    idz = with(tracker_df[x, ], which(names(list_pop_params)==paste0("z",z)));
  #    switched_df = lpa.mi.src::tracker2_resolve_label_switch(nn=x,
  #                                                            tracker_df=tracker_df,
  #                                                            data_conditions=data_conditions,
  #                                                            out_get_FMM = list_pop_params[[idz]])
  #    iter = iter+1
  #    diagnose_list[[iter]] = switched_df
  #  }
  #
  # diagnose_df = rbindlist(diagnose_list)

  print(" ")
  print("Identifying and remedying label switching...")
  if (is.null(cl)){
  list_switched = pblapply(X = 1:length(nn_vec),
                           FUN = function(x){
                             require(MplusAutomation); require(lpa.mi.src);
                             idz = with(tracker_df[x, ], which(names(list_pop_params)==paste0("z",z)));
                             switched_df = lpa.mi.src::tracker2_resolve_label_switch(oneline_df=tracker_df[x, ],
                                                                        data_conditions=data_conditions,
                                                                        out_get_FMM = list_pop_params[[idz]],
                                                                        args_KL_check_switching = args_KL_check_switching,
                                                                        args_tracker2inpfile_naive = args_tracker2inpfile_naive)
                             return(switched_df)
                           },
                           cl = NULL)
  } else {
    print(cl)
    pb <- pbapply::timerProgressBar(max = length(nn_vec), style = 1, width = getOption("width")/4)
    progress <- function(x){setTimerProgressBar(pb, x)}
    opts <- list(progress = progress)
    list_switched<-foreach(x =  1:length(nn_vec), .packages= c("MplusAutomation","lpa.mi.src"), .inorder = TRUE, .options.snow = opts) %dopar% {
                      idz = with(tracker_df[x, ], which(names(list_pop_params)==paste0("z",z)));
                      switched_df = lpa.mi.src::tracker2_resolve_label_switch(oneline_df = tracker_df[x, ],
                                                                              data_conditions=data_conditions,
                                                                              out_get_FMM = list_pop_params[[idz]],
                                                                              args_KL_check_switching = args_KL_check_switching,
                                                                              args_tracker2inpfile_naive = args_tracker2inpfile_naive)
                      return(switched_df)
    }
    close(pb)
  }
  switched_df = rbindlist(list_switched)
  tracker_df = merge(x = tracker_df, y = switched_df, by = "nn", all.x =  TRUE)
  tracker_df$converged[tracker_df$singular == TRUE] = FALSE

  N_switched = sum(as.numeric(tracker_df$switched), na.rm = TRUE)
  print(paste0("...",N_switched, " models displayed label switching."))

  N_singular = sum(as.numeric(tracker_df$singular), na.rm = TRUE)
  print(paste0("...",N_singular , " models that displayed convergence demonstrated at least one singular within-class covariance matrix. Convergence status changed appropriately."))

  return(tracker_df)
}
