#' Function to check Rcond and update it
#'
#' @param Plist (list)
#'
#' @return (list)
#' @export
#'
#' @examples

check_rcond_S<-function(Plist){

  if ( !("S"%in%names(Plist)) ){
    stop("S is not a variable in Plist. Check Plist$S exists and is an array of covariance matrices by class.")
  }
  if (length(dim(Plist$S))!=3){
    stop("Plist$S is not an array")
  }

  KK = dim(Plist$S)[3]
  poor_rcond = FALSE
  for (kk in 1:KK){
    rc_kk = rcond(Plist$S[,,kk])
    if (rc_kk<1E-4){
      poor_rcond = TRUE
      diag_kk = diag(Plist$S[,,kk])
      max_kk = max(diag_kk)
      if(max_kk<1E-3){
        max_kk = 1E-2
        Plist$S[which.max(diag_kk), which.max(diag_kk),kk] = max_kk
      }
      diag_kk[which.min(diag_kk)] = max_kk*1E-3
      for(j in 1:length(diag_kk)){
        Plist$S[j,j,kk] = diag_kk[j]
      }
    }
  }
  return(list(Plist=Plist, poor_rcond = poor_rcond))
}
