#' Simulate complete LPA data for simulation study
#'
#' This function simulates complete data.
#' @param z (integer) condition identifier for the complete data (as specified by data_conditions)
#' @param list_imputed (list) from calling get_imputed_data() function
#' @param methods_list (list) List w/ pva elements. Designates the imputation methods to try.
#' @param pctmiss_vec  (vector) Vector w/ pm elements. Designates percent missing simulation conditions
#' @param rep (integer) Replication number (defaults to 1).
#' @param save_it (logical) if TRUE, then it saves a .dat data set to load into Mplus; rep, p, & temp_wd_p_vec need to be specified.
#' @param p (integer) Processor number (defaults to 1).
#' @param temp_wd_p_vec (character vector). Required if save_it==TRUE.  processor-specific temporary directory
#' @return out_list  (list) with the following elements
#'       (A) stack_df - (data.frame) stacked imputed data set
#'       (b) dffolderfiles (data.frame) with the files and folders of the saved data
#' @export
#' @examples
#'get_stacked_data(z,data_conditions,methods_list,pct_miss_vec,save_it = TRUE)


get_stacked_data<-function(z, list_imputed, methods_list, pctmiss_vec,
                           save_it = FALSE, rep = NA, p = NA, temp_wd_p_vec = NULL){


  temp_wd_p = temp_wd_p_vec[p]

  require(tidyverse)
  require(mice)
  if (save_it==TRUE){require(MplusAutomation)}
  if(save_it == TRUE & (is.na(rep) | is.null(temp_wd_p) | is.na(p))){
    stop("rep, p, & temp_wd_p arguments cannot be NULL with save_it = TRUE")
  }


    vv=0
    out_list = list(stacked_list=NULL, dffolderfiles = NULL)
    for (pm in 1:length(pctmiss_vec)){
      for (pva in 1:length(methods_list$procedure)){

          idx_mids = which(names(list_imputed$obj_call)==paste0("mids_pm",pm,"pva",pva))

          mids = list_imputed$obj_call[[idx_mids]]
          stack_df = mice::complete(mids, action = "long", include = FALSE) %>% transform(wgt = 1/mids$m)
          ns = names(stack_df); ns[1] = "m"; ns[2] = "id";
          names(stack_df) = ns;


          cols_discard = which(ns == "m" | ns == "id")

          if (save_it){
            invisible(
              prepareMplusData(stack_df, keepCols = setdiff(1:ncol(stack_df), cols_discard),
                               filename = paste0(temp_wd_p,"/Stacked data/pm",pm,"/pva",pva,"/stack p", p," z",z," rep",rep, " pm", pm, " pva", pva, ".dat"),
                               inpfile = FALSE,
                               overwrite = TRUE)
            )
          }

          vv = vv+1
          out_list$stacked_list[[vv]] = stack_df
          names(out_list$stacked_list)[vv] = paste0("stack_pm",pm,"pva",pva)
          out_list$dffolderfiles = rbind(out_list$dffolderfiles,
                                         data.frame(folders = paste0("Stacked data/pm",pm,"/pva",pva),
                                                    files = paste("stack p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat",sep = ""),
                                                    data_condition = z, m = NA) )

      }
    }

    return(out_list)

}


