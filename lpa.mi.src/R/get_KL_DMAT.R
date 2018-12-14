#' Obtain K-L Divergence discrepancy matrix across classes using population values
#'
#' @param params () Parameter estimates obtained from Mplus
#' @param z (integer)
#' @param data_conditions (data.frame)
#' @return (matrix) of size K-by-K with the KL divergences
#' @export
#' @examples
#' param_complete = get_KL_DMAT(params, z, data_conditions)


get_KL_DMAT<-function(params, z, data_conditions){

  require(gaussDiff)

      # Get number of class indicators and number of classes based on data condtion z
      J_Y = data_conditions$J_Y[z]
      K = data_conditions$K[z]

      # Get population parameters based on the data condition z
      out_get_FMM<-get_FMM_params(z,data_conditions)


      # Identify the pairwise KL divergences across each population and sample/estimated classess
      DMAT = mat.or.vec(K,K) + NA
      for (k in 1:K){

        params_k = subset(params, LatentClass == k)
        mu_k = params_k$est[params_k$paramHeader=="Means"]
        S_k = diag(params_k$est[params_k$paramHeader == "Variances"])
        S_k[upper.tri(S_k)]<-S_k[lower.tri(S_k)]<-params_k$est[endsWith(params_k$paramHeader, ".WITH")]

        nvar = sum(params_k$est[params_k$paramHeader == "Variances"]<=0)
        if (nvar>0){
          out_list = list(params = params, class_switched = NA, converged = FALSE)
          return(out_list)
        }

        kl_vec = rep(NA,K)

        for (kk in 1:K){

          kl_vec[kk] = normdiff(mu1 = mu_k, sigma1 = S_k,
                                mu2 = out_get_FMM$mu_z[1:J_Y, kk], sigma2 = out_get_FMM$S_z[1:J_Y,1:J_Y,kk],
                                method = c("KL"))
          kl_vec[kk] = kl_vec[kk] +  normdiff(mu1 = out_get_FMM$mu_z[1:J_Y, kk], sigma1 = out_get_FMM$S_z[1:J_Y,1:J_Y,kk],
                                              mu2 = mu_k, sigma2 = S_k,
                                              method = c("KL"))
        }

        DMAT[,k] = kl_vec
      }

    return(DMAT)
}
