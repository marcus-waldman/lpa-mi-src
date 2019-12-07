#' Format population parameters to P list
#'
#'
#' @param z  (integer) condition identifier (number) for the complete data (as specified by data_conditions)
#' @param data_conditions  (data.frame) simulation conditions pertaining to the complete data
#' @param out_FMM  (list) Optional. Defaults to NULL. If NULL, then extracted from get_FMM_params
#' @return (list) with pi, mu, and S as arguments defining the mixture model
#' @export
#' @examples
#' get_Plist(z, data_conditions)

get_Plist<-function(z,data_conditions, out_FMM = NULL){

    if(is.null(out_FMM)){out_FMM = get_FMM_params(z, data_conditions)}

    if (!("J_Y_z" %in% names(out_FMM))){stop("J_Y_z not a variable in out_FMM")}
    if (!("K_z" %in% names(out_FMM))){stop("K_z not a variable in out_FMM")}


    J = out_FMM$J_Y_z
    K = out_FMM$K_z

    P = list(pi = NULL, mu = NULL, S = NULL)
    P$pi = out_FMM$pi_z
    P$mu = out_FMM$mu_z[1:J, 1:K]
    P$S = out_FMM$S_z[1:J,1:J,1:K]

    return(P)

}
