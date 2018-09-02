#' Get individual ICs for stacked type pooling method
#'
#' @param p  (integer) Processor number
#' @param temp_wd_p_vec
#' @param kk (integer) Fitted number of classes
#' @param z
#' @param rep
#' @param starts_txt
#' @param mids_pm_pva A mids object frmo get_imputed_ function
#' @param pm Defaults to NA
#' @param pva Defaults to NA
#' @return (data.frame) results_df
#' @export
#' @examples
#'
#'
#'
extract_ICs_stacked<-function(p, temp_wd_p_vec, kk, z, rep, starts_txt, mids_pm_pva, pm = NA, pva = NA){


  M = mids_pm_pva$m


  DATATYPE = "Stacked data"
  start_wd= getwd()
  ### Recored complete data results ####
  setwd(paste0(temp_wd_p_vec[p],"/", DATATYPE))
  if(!is.na(pm)){setwd(paste0("pm",pm))}
  if(!is.na(pva)){setwd(paste0("pva",pva))}
  target_wd = getwd()

  ### Record Stacked data ###

  out_stacked<-invisible(readModels(target = target_wd,
                                    what = c("warn_err","summaries"),
                                    quiet = TRUE))

  files_stacked = list.files(path = target_wd, pattern = ".out", full.names = TRUE)

  # Get conditioning number
  hi = readLines(con = files_stacked)
  i_0 = which(hi == "QUALITY OF NUMERICAL RESULTS")
  tmp = NA
  if(length(i_0)>0){
    tmp = str_extract_all(hi[i_0+2],"[0-9]+")[[1]]
    tmp = as.numeric(paste0(tmp[1],".",tmp[2],"E-",tmp[3]))
  }

  # record summary data
  sm = out_stacked$summaries
  results_df = data.frame(Pool_Method = NA, Parameters = NA, LL = NA, AIC = NA, BIC = NA, aBIC = NA, Entropy = NA, AICC = NA, Rcond = NA, Converged = NA, Error = NA, m = 0)

  LL_raw = ifelse(!is.null(sm$LL), sm$LL, NA);


  P = ifelse(!is.null(sm$Parameters), sm$Parameters, NA)
  LL = LL_raw/M
  N = sm$Observations/M
  Nstar = (N+2)/24;

  results_df$Parameters = P
  results_df$LL = LL
  results_df$AIC = 2*P - 2*LL
  results_df$BIC = log(N)*P - 2*LL
  results_df$AICC = NA  #Unknown if it is AIC or CAIC in MPlus
  results_df$aBIC = log(Nstar)*P - 2*LL
  results_df$Entropy = NA #Unknown if appropriate
  results_df$Rcond = tmp
  if (kk == 1){
    results_df$Converged = TRUE
  } else {
    results_df$Converged = check_convergence(files_stacked,
                                            folder_wd = paste0(temp_wd_p_vec[p],"/Stacked data/"),
                                            starts_txt = starts_txt,
                                            type_imputation = TRUE)
  }
  results_df$Error = ifelse(length(out_stacked$errors)>0,TRUE,FALSE)
  results_df$m = NA
  results_df$Pool_Method = as.factor("Stacked")

  results_df = transform(results_df, kk = kk, z = z, rep = rep, p = p)
  results_df = transform(results_df, pm = pm, pva = pva)


  return(results_df)

}

