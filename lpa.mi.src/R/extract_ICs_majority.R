#' Get individual ICs for majority vote
#'
#'
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
extract_ICs_majority<-function(p, temp_wd_p_vec, kk, z, rep, starts_txt, pm = NA, pva = NA){

      DATATYPE = "Imputed data"
      start_wd= getwd()
      ### Recored complete data results ####
      setwd(paste0(temp_wd_p_vec[p],"/", DATATYPE))
      if(!is.na(pm)){setwd(paste0("pm",pm))}
      if(!is.na(pva)){setwd(paste0("pva",pva))}
      target_wd = getwd()

     ### Record imputed data ###

    out_imputed<-invisible(readModels(target = target_wd,
                            what = c("warn_err","summaries")))

    imputed_files_m = list.files(path = target_wd, pattern = ".out")

    results_df =NULL
    for(m in 1:length(imputed_files_m)){
      results_m = data.frame(Pool_Method = NA, Parameters = NA, LL = NA, AIC = NA, BIC = NA, aBIC = NA, Entropy = NA, AICC = NA, Rcond = NA, Converged = NA, Error = NA, m = NA)

      file_m = imputed_files_m[which(endsWith(imputed_files_m,paste0("imp_",m,".out")))]

      hi = readLines(con = paste(target_wd,file_m, sep = "/"))
      i_0 = which(hi == "QUALITY OF NUMERICAL RESULTS")
      tmp = NA
      if(length(i_0)>0){
        tmp = str_extract_all(hi[i_0+2],"[0-9]+")[[1]]
        tmp = as.numeric(paste0(tmp[1],".",tmp[2],"E-",tmp[3]))
      }

      out_m = out_imputed[[which(endsWith(names(out_imputed),paste0("imp_",m,".out")))]]

      sm = out_m$summaries
      results_m$Parameters = ifelse(!is.null(sm$Parameters), sm$Parameters, NA)
      results_m$LL = ifelse(!is.null(sm$LL), sm$LL, NA)
      results_m$AIC = ifelse(!is.null(sm$AIC), sm$AIC, NA)
      results_m$BIC = ifelse(!is.null(sm$BIC), sm$BIC, NA)
      results_m$aBIC = ifelse(!is.null(sm$aBIC), sm$aBIC, NA)
      results_m$Entropy = ifelse(!is.null(sm$Entropy), sm$Entropy, NA)
      results_m$AICC = ifelse(!is.null(sm$AICC), sm$AICC, NA)
      results_m$Rcond = tmp
      results_m$Converged = check_convergence(file_m,
                                              folder_wd =  target_wd,
                                              starts_txt = starts_txt)
      results_m$Error = ifelse(length(out_m$errors)>0,TRUE,FALSE)
      results_m$m = m
      results_m$Pool_Method = "Majority"
      results_df = rbind(results_df, results_m)
    }

    results_df = transform(results_df, kk = kk, z = z, rep = rep, p = p)
    results_df$Pool_Method = as.factor(results_m$Pool_Method)

    results_df = transform(results_df, pm = pm, pva = pva)


    setwd(start_wd)
    return(results_df)

}

