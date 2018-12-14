gamma_vec = c(1,1) #K-1 size vector with conditional log odds

K = length(gamma_vec)+1
pi_K = 1/(1 + sum(exp(gamma_vec)))

pi_vec = rep(NA, K)
pi_vec[seq(1,K-1)] = pi_K*exp(gamma_vec)
pi_vec[K] = pi_K

# return(pi_vec)