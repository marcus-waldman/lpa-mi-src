#' Converts from Plist to get_FMM list needed to make input files
#'
#' @param Plist_x 
#' @return 
#' @export
#' @examples

Plist2getFMM<-function(Plist_x){

    get_FMM_x = list(mu_z = NULL, 
                     S_z = NULL, 
                     pi_z = NULL, 
                     K_z = NULL, 
                     J_Y_z = NULL, 
                     J_Xinc_z = 0, 
                     J_Xcom_z = 0, 
                     J = NULL, 
                     MD_z = NA, 
                     rho_YX_z = NA, 
                     dX_z = NA)
    
    
    get_FMM_x$mu_z = Plist_x$mu
    get_FMM_x$S_z = Plist_x$S
    get_FMM_x$pi_z = Plist_x$pi
    get_FMM_x$K_z = dim(Plist_x$S)[3]
    get_FMM_x$J_Y_z<-get_FMM_x$J<-nrow(Plist_x$S[,,1])

    return(get_FMM_x)

}
