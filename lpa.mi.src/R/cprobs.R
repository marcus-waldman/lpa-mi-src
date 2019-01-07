#' Calculate posterior class probabilities
#'
#' This function calculates the posterior class probabiltieis from a finite mixture of multivariate normal distributions.
#' @param Y_i (matrix) of size NxJ of observed values
#' @param pi_vec (numeric vector) with K elements containing component probabilities
#' @param mu_mat (matrix) of size JxK, where J is number of items
#' @param S_array (array) of size JxJxK with within class covariance matrix.
#' @return (matrix) of size NxK containing the posterior class porbabilties
#' @export
#' @examples
#' cprobs(Y_i, pi_vec, mu_mat, S_array)


dmixmvrnorm<-function(Y_i, pi_vec, mu_mat, S_array, log = FALSE){
    require(mixtools)
    require(pracma)

    K = length(pi_vec)
    J = ncol(Y_i)
    N = nrow(Y_i)

    f_mat = mat.or.vec(nr = N, nc = K) + NA
    num_i = mat.or.vec(nr = N, nc = K-1) + NA
    for (k in 1:K){
      f_mat[1:N, k] = dmvnorm(y = Y_i, mu = as.vector(mu_mat[1:J,k]), sigma = as.matrix(S_array[1:J,1:J,k]))
      if (k<K){
        num_i[1:N,k] = f_mat[1:N,k]
      }
    }

    d_i = repmat(f_mat%*%as.numeric(pi_vec), n = 1, m = K-1)

    cprob_i = num_i/d_i
    Kprob_i = 1- apply(cprob_i,1,"sum")
    cprob_i = cbind(cprob_i, Kprob_i)


    return(cprob_i)

}
