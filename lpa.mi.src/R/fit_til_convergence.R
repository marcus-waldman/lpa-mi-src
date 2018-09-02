#' Fit mixture models in Mplus until convergence
#'
#'
#' @param p Processor
#' @param z Data conditioning number
#' @param temp_wd_p_vec Processor-specific folder
#' @param target_wd Character of where data is located
#' @param dff_target List of files from the get_*_data functions
#' @param pop_params_kk Population parameters
#' @param starts0 Must be a multiple of 10. Defaults to 20.
#' @param mX_max Highest exponents to go to (ie.e max mX for 2^mx)
#' @param type_imputation (logical) passed to create_NaiveMplus_inpfile() or check_convergence()
#' @param ... Other arguments to send to create_naiveMplus_inpfile()
#' @return (data.frame) results_df
#' @export
#' @examples
#' fit_til_convergence<(z, target_wd, dff_target, pop_params_kk)



fit_til_convergence<-function(p, z, temp_wd_p_vec, target_wd, dff_target, pop_params_kk, starts0 = 20, mX_max = 7, type_imputation = FALSE, m = NA, ...){


    if(starts0%%10 != 0){error("starts0 needs to be a multiple of 10")}

    mX = -1;
    terminate = FALSE
    while(terminate==FALSE){
        mX = mX+1

        starts_txt = paste0(starts0*2^(mX)," ",starts0*(4/10)*2^(mX),";")

        pattern_inp = ifelse(is.na(m), ".inp", paste0(m,".inp"))
        pattern_out = ifelse(is.na(m), ".out", paste0(m,".out"))

        file.remove(list.files(path = target_wd  , pattern = pattern_inp, full.names = TRUE))
        file.remove(list.files(path = target_wd  , pattern = pattern_out, full.names = TRUE))

        inp_complete = do.call(what = "create_naiveMplus_inpfile",
                              args = list(z = z, out_get_FMM = pop_params_kk,
                                          dffolderfiles = dff_target, temp_wd_p = temp_wd_p_vec[p],
                                          starts_txt = starts_txt, type_imputation = type_imputation))

        runModels(target = target_wd, filefilter = list.files(path = target_wd, pattern = pattern_inp),logFile = NULL)
        converged = check_convergence(file = list.files(path = target_wd, pattern = pattern_out),
                          folder_wd = target_wd,
                          starts_txt = starts_txt,
                          type_imputation = type_imputation)

        if (converged == TRUE | mX >= mX_max){
          terminate = TRUE
        }

    }

    # Get conditioning number
    hi = readLines(con = list.files(path = target_wd, pattern = pattern_out, full.names = TRUE))
    i_0 = which(hi == "QUALITY OF NUMERICAL RESULTS")
    tmp = NA
    if(length(i_0)>0){
      tmp = str_extract_all(hi[i_0+2],"[0-9]+")[[1]]
      tmp = as.numeric(paste0(tmp[1],".",tmp[2],"E-",tmp[3]))
    }

    problem = (mX>=mX_max & converged == FALSE)

    out_list = list(Converged = converged, Rcond = tmp, Starts = starts0*2^(mX), problem = problem)
    return(out_list)
}
