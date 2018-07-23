#' Get results 
#'
#' 
#' @param p  (integer) Processor number
#' @param temp_wd_p_vec
#' @param kk (integer) Fitted number of classes
#' @param z 
#' @param rep
#' @param starts_txt
#' @return (data.frame) results_df
#' @export
#' @examples
#' 
#' 
#' 
results_pool<-function(p, temp_wd_p_vec, kk, z, rep, starts_txt){

    ### Recored complete data results ####
    
    
    out_complete<-invisible(readModels(target = paste0(temp_wd_p_vec[p],"/Complete data") , 
                             what = c("warn_err","summaries")))
    files_complete = list.files(paste0(temp_wd_p_vec[p],"/Complete data"), pattern = ".out")
    
    # Get conditioning number
    hi = readLines(con = paste(temp_wd_p_vec[p],"/Complete data/",files_complete, sep = ""))
    i_0 = which(hi == "QUALITY OF NUMERICAL RESULTS")
    tmp = NA
    if(length(i_0)>0){ 
      tmp = str_extract_all(hi[i_0+2],"[0-9]+")[[1]]
      tmp = as.numeric(paste0(tmp[1],".",tmp[2],"E-",tmp[3]))
    }
    
    # record summary data
    sm = out_complete$summaries
    results_m = data.frame(Parameters = NA, LL = NA, AIC = NA, BIC = NA, aBIC = NA, Entropy = NA, AICC = NA, Rcond = NA, Converged = NA, Error = NA, m = 0)
    results_m$Parameters = ifelse(!is.null(sm$Parameters), sm$Parameters, NA)
    results_m$LL = ifelse(!is.null(sm$LL), sm$LL, NA)
    results_m$AIC = ifelse(!is.null(sm$AIC), sm$AIC, NA)
    results_m$BIC = ifelse(!is.null(sm$BIC), sm$BIC, NA)
    results_m$aBIC = ifelse(!is.null(sm$aBIC), sm$aBIC, NA)
    results_m$Entropy = ifelse(!is.null(sm$Entropy), sm$Entropy, NA)
    results_m$AICC = ifelse(!is.null(sm$AICC), sm$AICC, NA)
    results_m$Rcond = tmp
    if (kk == 1){
      results_m$Converged = TRUE
    } else {
      results_m$Converged = check_convergence(files_complete, 
                                              folder_wd = paste0(temp_wd_p_vec[p],"/Complete data/"), 
                                              starts_txt = starts_txt) 
    }
    results_m$Error = ifelse(length(out_complete$errors)>0,TRUE,FALSE)
    results_m$m = NA
    
    ### Record imputed data ###
    
    out_imputed<-invisible(readModels(target = paste0(temp_wd_p_vec[p],"/Imputed data"), 
                            what = c("warn_err","summaries")))
    
    imputed_files_m = list.files(path = paste0(temp_wd_p_vec[p],"/Imputed data/"), pattern = ".out")
    
    results_df = results_m
    for(m in 1:length(imputed_files_m)){
      results_m = data.frame(Parameters = NA, LL = NA, AIC = NA, BIC = NA, aBIC = NA, Entropy = NA, AICC = NA, Rcond = NA, Converged = NA, Error = NA, m = NA)
      
      file_m = imputed_files_m[which(endsWith(imputed_files_m,paste0("imp_",m,".out")))]
      
      hi = readLines(con = paste(temp_wd_p_vec[p],"/Imputed data/",file_m, sep = ""))
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
                                              folder_wd =  paste0(temp_wd_p_vec[p],"/Imputed data/"), 
                                              starts_txt = starts_txt) 
      results_m$Error = ifelse(length(out_complete$errors)>0,TRUE,FALSE)
      results_m$m = m
      results_df = rbind(results_df, results_m)
    }
    
    results_df = transform(results_df, kk = kk, z = z, rep = rep, p = p)
    
    return(results_df)

}

