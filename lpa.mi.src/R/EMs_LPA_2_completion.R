#' Imputation by EM with Sampling
#'
#' @param obsdf (data.frame) with the observed data
#' @param M (integer) Number of imputations
#' @param z (integer)
#' @param data_conditions (data.frame)
#' @param tmax (numeric). Defaults to a maximum of 30 minutes to achieve the imputations.
#' @param ... Optional arguments to EMs_LPA() function (e.g., starts0, bayes_boot, boot_2x, tempdir_mplus, rep, itermax, multiplier, cl)
#'
#' @return A mids object with the imputations
#' @export
#'
#' @examples
#' EMs_LPA_2_completion<-function(obsdf, M, z, data_conditions)
#'
EMs_LPA_2_completion<-function(obsdf, M, z, data_conditions, tmax = 30, ...){

    require("mice")
    require("lpa.mi.src")

    tic = proc.time()

    MM = M; flag = 0; iter = 0
    while(flag == 0){


      args_EMs = list(obsdf=obsdf, M=MM, z=z, data_conditions=data_conditions, ...)
      out_mids = do.call(what = "EMs_LPA", args = args_EMs)

      if(!is.null(out_mids)){
        iter = iter+1
        if (iter==1){
          EMs_mids = out_mids
        } else {
          EMs_mids = ibind(EMs_mids, out_mids)
        }
        MM = M - EMs_mids$m
      }

      toc = proc.time()-tic
      toc = toc[[3]]/60

      if (MM <= 0){
        flag = 1
      } else {
        if (toc>=tmax){warning("Time exceeded. Breaking EMs_LPA-2_completion without full imputations"); flag=1;}
        rm(out_mids)
      }

    }
    EMs_mids$method = gsub(pattern = "pmm", replacement = "EMs_LPA", EMs_mids$method)
    return(EMs_mids)

}
