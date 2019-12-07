#' Conditional distribution of lower-dimensional MVN from higher-dimensional MVN
#'
#' @param mu (numeric) Mean vector of full joint MVN distribution
#' @param S (matrix) Covariance matrix of full joint MVN distribution
#' @param U (numeric) Indices of missing variables
#' @param V (numeric) Indices of observed variables
#'
#' @return A list with the following elements multivariate regression covefficients (beta) and the residual variance ovariance matrics S_eps
#' @export
#'
#' @examples get_conditional(mu,S,U,V)


get_conditional<-function(mu,S,U,V){

  #U are/is the index ofthe missing variable(s) (numeric)
  #V are/is the index of the observed variable(s) (numeric)
  # Formulas taken from BDA pg. 582
  # mu = mu_mat[1,]
  # S = S_array[,,1]
  # U = c(1,2) # 3 #unobserved
  # V = 3 # c(1,2) # observed

  if (sum(U%in%V)>0){stop("U and V must be disjoint")}
  E_U = mu[U]
  cov_UV = S[U,V]
  var_V = S[V,V]
  i_var_V = solve(var_V)
  E_V = mu[V]

  t_beta_vec = cov_UV%*%i_var_V
  t_beta_vec_cons = E_U - t_beta_vec%*%E_V
  beta = rbind(t(t_beta_vec_cons), t(t_beta_vec))


  var_U = S[U,U]
  cov_VU = S[V,U]
  S_eps = var_U - t_beta_vec%*%cov_VU



  return(list(beta = beta, S_eps = S_eps))

}
