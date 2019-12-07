#' tracker2fitfixed
#'
#' Creates an input file that fits an LPA model fixed at the pooled parameter estiamtes given information from a tracker data frame.
#' @param oneline_df (data.frame) One row from the tracker data frame for monitoring covernce, specifying dat file, specifying .inp file path, etc.
#' @param args_tracker2inpfile_naive (list) Optional arguments passed to tracker2inpfile_naive() beyond nn, tracker_df, pop_params_nn, data_conditions, Model_txt, process_cores
#' @param showOutput (logical) Defaluts to FALSE. This sets the showOutput argument in MplusAutomation::runModels()
#' @return Saves the Mplus input files for Naive LPA analysis in the the folders specified of dffolderfiles input. Returns a character vector with the the Mplus code.
#' @export
#' @examples
#' tracker2fitfixed(oneline_df)

tracker2fitfixed<-function(oneline_df,
                           args_tracker2inpfile_naive = list(NULL),
                           showOutput = FALSE){

  require(tidyverse)

  load(file = paste(oneline_df$poolwd,oneline_df$poolfolder,oneline_df$poolfile, sep = "/"))
  M_x = ifelse(oneline_df$M_condition>ncol(pooled_list$details$theta_mat), ncol(pooled_list$details$theta_mat), oneline_df$M_condition)

  t_bar_x = apply(pooled_list$details$theta_mat[,1:M_x],1,"mean")
  pooled_parameters_x = pooled_list$details$list_params$mm1[,c("paramHeader","param","LatentClass")] %>%
    transform(est = t_bar_x)
  Qlist_x = lpa.mi.src::Mplus2Qlist(params_df = pooled_parameters_x)
  svals_x =  lpa.mi.src::Qlist2syntax(Qlist = Qlist_x, symbol = "@")

  args = args_tracker2inpfile_naive
  inp_txt = tracker2inpfile_naive(oneline_df = oneline_df,
                                  pop_params_nn = Plist2getFMM(Qlist_x),
                                  savedata = ifelse(is.null(args$savedata), FALSE, args$savedata),
                                  estimates = ifelse(is.null(args$estimates), FALSE, args$estimates),
                                  results = ifelse(is.null(args$results), FALSE, args$results),
                                  save_tech3 = ifelse(is.null(args$save_tech3), FALSE, args$save_tech3),
                                  symbol_txt = ifelse(is.null(args$symbol_txt), "@",args$symbol_txt),
                                  Model_txt = c("MODEL:", svals_x),
                                  weight = ifelse(is.null(args$weight), FALSE, args$weight),
                                  output_txt = args$output_txt,
                                  processor_cores = 1)

  with(oneline_df, runModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                             filefilter = as.character(paste0(outfile,".inp")),
                             logFile = NULL,
                             showOutput = showOutput))

}
