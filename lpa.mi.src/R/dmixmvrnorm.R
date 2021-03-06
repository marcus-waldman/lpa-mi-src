#' Calculate (log) density from a mixture of MVNs
#'
#' This function calculates the (log) density from a finite mixture of multivariate normal distributions.
#' @param Y_i (matrix) of size NxJ of observed values
#' @param pi_vec (numeric vector) with K elements containing component probabilities
#' @param mu_mat (matrix) of size JxK, where J is number of items
#' @param S_array (array) of size JxJxK with within class covariance matrix.
#' @param log (logical) Whether the log density should be return. Defaults to FALSE.
#' @return (vector) of size N containing (log) density values
#' @export
#' @examples
#' dmixmvrnorm(Y_i, pi_vec, mu_mat, S_array, log = TRUE)


dmixmvrnorm<-function(Y_i, pi_vec, mu_mat, S_array, log = FALSE){
    require(mixtools)

    K = length(pi_vec)
    J = ncol(Y_i)
    N = nrow(Y_i)

    f_mat = mat.or.vec(nr = N, nc = K) + NA
    for (k in 1:K){
      f_mat[1:N, k] = dmvnorm(y = Y_i, mu = as.vector(mu_mat[1:J,k]), sigma = as.matrix(S_array[1:J,1:J,k]))
    }

    out_i = as.vector(f_mat%*%as.numeric(pi_vec))
    if (log == TRUE){out_i = log(out_i)}

    return(out_i)

}
