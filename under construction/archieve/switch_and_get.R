#' Resolve label switching and extracts parameter estimates.
#'
#' @param out_ftc (list) An object from fit_til_convergence. Currently the options readModels = TRUE and savedata = TRUE must be specified in fit_til_convergence.
#' @param z (integer) Data conditioning number
#' @param data_conditions (data.frame) 
#' @return (data.frame) with parameter estimates
#' @export
#' @examples
#' param_complete = switch_and_get(out_ftc = out_ftc_complete, z = z, data_conditions = data_conditions)


switch_and_get<-function(out_ftc, z, data_conditions){
    agg_cprob = ddply(out_ftc$out_readModels$savedata[,c("SUBPOP","CPROB1","CPROB2","CPROB3")],
                .(SUBPOP), summarize, CPROB1 = mean(CPROB1), CPROB2 = mean(CPROB2), CPROB3 = mean(CPROB3))
    DMAT = 1-as.matrix(agg_cprob[,-1])
    out_switch = resolve_label_switch(mu_est_mat = NULL, S_est_array = NULL, z = z, data_conditions = data_conditions, 
                         parameters_df = out_ftc$out_readModels$parameters$unstandardized, 
                         DMAT = DMAT)
    return(out_switch)
}
