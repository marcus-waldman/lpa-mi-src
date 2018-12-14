#' Samples from mixture of MVNs
#'
#' This function draws samples from a finite mixture of MVNs
#' @param N (integer) Number of samples to draw
#' @param pi_vec (numeric vector) with K elements containing component probabilities
#' @param mu_mat (matrix) of size JxK, where J is number of items
#' @param S_array (array) of size JxJxK with within class covariance matrix.
#' @return (matrix) of size NxJ  containing draws from finite ixture of MVNs
#' @export
#' @examples
#' rmixmvrnorm(N, pi_vec, mu_mat, S_array)

rmixmvrnorm<-function(N, pi_vec, mu_mat, S_array){


  require(MASS)

  J = nrow(mu_mat)
  K = length(pi_vec)

  # Simulate class membership
  kappa_smart = t(rmultinom(n = N, size = 1, prob = pi_vec))
  kappa_i = apply(kappa_smart, MARGIN = 1, FUN = which.max)

  # Simulate outcomes
  Y_i = mat.or.vec(nr = N, nc = J) + NA
  for (k in 1:K){
    inds_k = which(kappa_i == k)
    N_k = length(inds_k)
    if (N_k>0){
      Y_i[inds_k, 1:J] = mvrnorm(n = N_k, mu = as.vector(mu_mat[1:J,k]), Sigma = as.matrix(S_array[1:J,1:J,k]))
    }
  }

  return(Y_i)

}
