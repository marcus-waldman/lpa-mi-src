#' Convert conditional log-odds to probabilities
#'
#' 
#' @param gamma_vec (numeric vector) of size K-1 with conditional log odds
#' @return (numberic vector) of size K with corresponsding probabilites
#' @export
#' @examples
#' gamma2pi(c(0,0))

gamma2pi<-function(gamma_vec){

    K = length(gamma_vec)+1
    pi_K = 1/(1 + sum(exp(gamma_vec)))
    
    pi_vec = rep(NA, K)
    pi_vec[seq(1,K-1)] = pi_K*exp(gamma_vec)
    pi_vec[K] = pi_K
    
    return(pi_vec)
}

