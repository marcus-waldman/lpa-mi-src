# Inputs Y_i (NxJ), pi_vec (size K), mu_mat (JxK), S_array(JxJxK)


require(mixtools)

K = length(pi_vec)
J = ncol(Y_i)
N = nrow(Y_i)



f_mat = mat.or.vec(nr = N, nc = K) + NA 
for (k in 1:K){
  f_mat[1:N, k] = dmvnorm(y = Y_i, mu = as.vector(mu_mat[1:J,k]), sigma = as.matrix(S_array[1:J,1:J,k]))
}

logf_i = as.vector(log(f_mat%*%pi_vec))

# return logf_i