#' Calculate posterior class probabilities
#'
#' This function calculates the posterior class probabiltieis from a finite mixture of multivariate normal distributions.
#' @param Y_i (matrix) of size NxJ of observed values
#' @param pi_vec (numeric vector) with K elements containing component probabilities
#' @param mu_mat (matrix) of size JxK, where J is number of items
#' @param S_array (array) of size JxJxK with within class covariance matrix.
#' @param list_mdpattern (list) Optional. Generated from get_mdpatterns function.
#' @return (matrix) of size NxK containing the posterior class porbabilties
#' @export
#' @examples
#' cprobs(Y_i, pi_vec, mu_mat, S_array)


cprobs<-function(Y_i, pi_vec, mu_mat, S_array, list_mdpattern = NULL, log = FALSE){

    require(Rfast)
    require(pracma)
    require(dynr)
    require(mixtools)
    require(pkgmaker)

    K = length(pi_vec)
    J = ncol(Y_i)
    N = nrow(Y_i)

    if (is.null(list_mdpattern)){list_mdpattern = get_mdpatterns(df = Y_i)}

    fp_mat = mat.or.vec(nr = N, nc = K) + NA
    for (k in 1:K){
      for (r in list_mdpattern$table_r$idR){
        ix_r = which(list_mdpattern$patterns_i$idR == r)
        if(r==0){idV=1:J}else{idV=list_mdpattern$V_r[[r]]}
        if(length(ix_r)==1){Y_r = as.numeric(Y_i[ix_r,idV])}else{Y_r= as.matrix(Y_i[ix_r,idV])}

        S_rk = S_array[idV,idV,k]
        if(!isReal(S_rk)){
          if(nrow(S_rk)==1){
            S_rk = as.numeric(S_rk)
            if (S_rk<1E-5){S_rk = 1E-5}
          } else{
            S_rk = as.matrix(S_rk)
            if (rcond(S_rk)<1E-5){
              l = dynr::dynr.ldl(S_rk)
              e_vec = diag(l)
              e_vec[which.min(e_vec)] = max(e_vec)*1E-5
              d = diag(e_vec)
              diag(l)<-1
              S_rk = l%*%d%*%t(l)
            }
          }
        } else {
          if (S_rk<1E-5){S_rk = 1E-5}
        }

        #fp_mat[ix_r, k] = Rfast::dmvnorm(x = Y_r, mu = as.numeric(mu_mat[idV,k]),
        #                          sigma = S_rk, logged = FALSE)*pi_vec[k]
        fp_mat[ix_r, k] = mixtools::dmvnorm(y = Y_r, mu = as.numeric(mu_mat[idV,k]),
                                            sigma = S_rk)*pi_vec[k]

      } #end for r
    } # end for k
    denom_mat = t(repmat(apply(fp_mat,1,"sum"), n = K, m = 1))
    cprob_i = fp_mat/denom_mat
    return(cprob_i)

}
