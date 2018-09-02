#' Get ICs
#'
#' @param TYPE (character) Pooling type. Defaults to "NONE"
#' @param DATATYPE (charcter) One of (a) "Complete data", (b) "Observed data", (c) "Imputed data", or (d) "Stacked data"
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

extract_ICs_naive<-function(TYPE = "None", DATATYPE=NULL,p, temp_wd_p_vec, kk, z, rep, starts_txt, pm = NA, pva = NA, type_imputation = FALSE){

  start_wd= getwd()
  ### Recored complete data results ####
  setwd(paste0(temp_wd_p_vec[p],"/", DATATYPE))
  if(!is.na(pm)){setwd(paste0("pm",pm))}
  if(!is.na(pva)){setwd(paste0("pva",pva))}
  target_wd = getwd()

  files_naive = list.files(target_wd, pattern = ".out")
  results_out = NULL
  for (i in 1:length(files_naive)){

    out_i = gsub(".out","",files_naive[i])

    out_naive<-readModels(target = target_wd, what = c("warn_err","summaries"),
                          filefilter = out_i)

    # Get conditioning number
    hi = readLines(con = files_naive[i])
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
    if (type_imputation == FALSE){
      results_df$LL = ifelse(!is.null(sm$LL), sm$LL, NA)
      results_df$AIC = ifelse(!is.null(sm$AIC), sm$AIC, NA)
      results_df$BIC = ifelse(!is.null(sm$BIC), sm$BIC, NA)
      results_df$aBIC = ifelse(!is.null(sm$aBIC), sm$aBIC, NA)
      results_df$Entropy = ifelse(!is.null(sm$Entropy), sm$Entropy, NA)
      results_df$AICC = ifelse(!is.null(sm$AICC), sm$AICC, NA)
    }
    if(type_imputation == TRUE){
      results_df$LL = ifelse(!is.null(sm$LL_Mean), sm$LL_Mean, NA)
      results_df$AIC = ifelse(!is.null(sm$AIC_Mean), sm$AIC_Mean, NA)
      results_df$BIC = ifelse(!is.null(sm$BIC_Mean), sm$BIC_Mean, NA)
      results_df$aBIC = ifelse(!is.null(sm$aBIC_Mean), sm$aBIC_Mean, NA)
      results_df$Entropy = NA
      results_df$AICC = NA
    }
    results_df$Rcond = tmp
    if (kk == 1){
      results_df$Converged = TRUE
    } else {
      results_df$Converged = check_convergence(files_naive[i],
                                              folder_wd = target_wd,
                                              starts_txt = starts_txt,
                                              type_imputation = type_imputation)
    }
    results_df$Error = ifelse(length(out_naive$errors)>0,TRUE,FALSE)
    results_df$m = NA
    results_df$Pool_Method = as.factor(TYPE)

    results_df = transform(results_df, kk = kk, z = z, rep = rep, p = p)
    results_df = transform(results_df, pm = pm, pva = pva)

    results_out = rbind(results_out, results_df)

  }

  setwd(start_wd)
  return(results_out)
}
