#' Get ICs
#'
#' @param is_Complete  (logical) If true, then looks in ~/Complete data folder; otherwise, reads out of ~/Observed data
#' @param p  (integer) Processor number
#' @param temp_wd_p_vec
#' @param kk (integer) Fitted number of classes
#' @param z
#' @param rep
#' @param starts_txt
#' @param pm Defaults to NA
#' @param pva Defaults to NA
#' @return (data.frame) results_df
#' @export
#' @examples
#'
#'
#'

extract_ICs_naive<-function(is_Complete,p, temp_wd_p_vec, kk, z, rep, starts_txt, pm = NA, pva = NA){

  if (is_Complete == TRUE){DATATYPE = "Complete data"; TYPE = "Complete"}
  if (is_Complete == FALSE){DATATYPE = "Observed data"; TYPE = "Observed"}

  start_wd= getwd()
  ### Recored complete data results ####
  setwd(paste0(temp_wd_p_vec[p],"/", DATATYPE))
  if(!is.na(pm)){setwd(paste0("pm",pm))}
  if(!is.na(pva)){setwd(paste0("pva",pva))}
  target_wd = getwd()


  out_naive<-invisible(readModels(target = target_wd ,
                                     what = c("warn_err","summaries")))
  files_naive = list.files(target_wd, pattern = ".out", full.names = TRUE)

  # Get conditioning number
  hi = readLines(con = files_naive)
  i_0 = which(hi == "QUALITY OF NUMERICAL RESULTS")
  tmp = NA
  if(length(i_0)>0){
    tmp = str_extract_all(hi[i_0+2],"[0-9]+")[[1]]
    tmp = as.numeric(paste0(tmp[1],".",tmp[2],"E-",tmp[3]))
  }

  # record summary data
  sm = out_naive$summaries
  results_df = data.frame(Pool_Method = NA, Parameters = NA, LL = NA, AIC = NA, BIC = NA, aBIC = NA, Entropy = NA, AICC = NA, Rcond = NA, Converged = NA, Error = NA, m = 0)
  results_df$Parameters = ifelse(!is.null(sm$Parameters), sm$Parameters, NA)
  results_df$LL = ifelse(!is.null(sm$LL), sm$LL, NA)
  results_df$AIC = ifelse(!is.null(sm$AIC), sm$AIC, NA)
  results_df$BIC = ifelse(!is.null(sm$BIC), sm$BIC, NA)
  results_df$aBIC = ifelse(!is.null(sm$aBIC), sm$aBIC, NA)
  results_df$Entropy = ifelse(!is.null(sm$Entropy), sm$Entropy, NA)
  results_df$AICC = ifelse(!is.null(sm$AICC), sm$AICC, NA)
  results_df$Rcond = tmp
  if (kk == 1){
    results_df$Converged = TRUE
  } else {
    results_df$Converged = check_convergence(files_naive,
                                            folder_wd = paste0(temp_wd_p_vec[p],"/",DATATYPE),
                                            starts_txt = starts_txt)
  }
  results_df$Error = ifelse(length(out_naive$errors)>0,TRUE,FALSE)
  results_df$m = NA
  results_df$Pool_Method = as.factor(TYPE)

  results_df = transform(results_df, kk = kk, z = z, rep = rep, p = p)
  results_df = transform(results_df, pm = pm, pva = pva)

  setwd(start_wd)
  return(results_df)
}
