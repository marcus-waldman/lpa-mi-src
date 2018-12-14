#' Format Mplus Parameters to Q list
#'
#'
#' @param params_df (data.frame) Unstandardized parameters extracted from MplusAutomation
#' @return (list) with pi, mu, and S as arguments defining the mixture model
#' @export
#' @examples
#' Mplus2Q(params_df)

Mplus2Qlist<-function(params_df){

    K = max(as.numeric(subset(params_df, LatentClass !=  "Categorical.Latent.Variables")$LatentClass))
    J = length(subset(params_df, paramHeader == "Means" & LatentClass == "1")$est)

    Q = list(pi = NULL, mu = NULL, S = NULL)


    # Get mixing proportions
    gamma_vec = subset(params_df, paramHeader == "Means" & LatentClass == "Categorical.Latent.Variables")$est
    pi_vec = gamma2pi(gamma_vec)
    Q$pi = pi_vec


    # Get means
    mu_mat = mat.or.vec(nr = J, nc = K)+NA
    for(k in 1:K){
      mu_mat[1:J,k] = subset(params_df, paramHeader == "Means" & LatentClass == as.character(k))$est
    }
    Q$mu = mu_mat

    # Get Sigma
    S_array = array(dim = c(J,J,K))
    for (k in 1:K){
      for (j in 1:J){
        inds_jk = with(params_df, which(paramHeader == "Variances" & LatentClass == paste0(k) & param == paste0("Y",j)))
        S_array[j,j,k] = params_df$est[inds_jk]
      }
    }
    for(k in 1:K){
      for (j1 in seq(1,K-1)){
        for (j2 in seq(j1+1, K)){
          inds_j1j2k = with(params_df, which(paramHeader == paste0("Y",j1,".WITH") & param == paste0("Y",j2) & LatentClass == paste0(k)))
          S_array[j1,j2,k] <- S_array[j2,j1,k] <- params_df$est[inds_j1j2k]
        }
      }
    }
    Q$S = S_array

    return(Q)

}
