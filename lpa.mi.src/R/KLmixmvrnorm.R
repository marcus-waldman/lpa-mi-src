#' Kullback-Leibler divergence between MVN mixtures P & Q (e.g. D_KL(P||Q) )
#'
#'
#' This function uses Monte Carlo integration to approximate the K-L divergence from (the estimated) distribution Q to (the population) distribution P.
#' @param P (list) with the following named elements defining the population:
#'                  (1) pi-(numeric vector) with K elements containing component probabilities
#'                  (2) mu-(matrix) of size JxK, where J is number of items
#'                  (3) S-(array) of size JxJxK with within class covariance matrix.
#' @param Q (list) with the same named elements as above, except with the estimates instead of the population values.
#' @param n (integer) defining the number of Monte Carlo integration points to use (defaults to 5E5).
#' @return (data.frame) with the estimate (est) and standard error (se) of the K-L approximation.
#' @export
#' @examples
#' KLmixmvrnorm(P, Q, n = 1E6)

KLmixmvrnorm<-function(P, Q, n = 1E4){

  # Simulate obserations from P for Monte Carlo integratoin
  Y_i = rmixmvrnorm(N = n, pi_vec = P$pi, mu_mat = P$mu, S_array = P$S)

  # Calculate respective densities
  logp_i = dmixmvrnorm(Y_i = Y_i, pi_vec = P$pi, mu_mat = P$mu, S_array = P$S, log = TRUE)
  logq_i = dmixmvrnorm(Y_i = Y_i, pi_vec = Q$pi, mu_mat = Q$mu, S_array = Q$S, log = TRUE)

  # Calculate the difference in the log densities
  delta_i = logp_i - logq_i

  # Obtain KL divergence estimate and standard error
  KL_est = mean(delta_i)
  KL_se = sd(delta_i)/sqrt(n)

  KL_df = data.frame(est = KL_est, se = KL_se)

  return(KL_df)
}
