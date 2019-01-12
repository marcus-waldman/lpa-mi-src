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
        fp_mat[ix_r, k] = Rfast::dmvnorm(x = Y_r, mu = as.numeric(mu_mat[idV,k]),
                                  sigma = as.matrix(S_array[idV,idV,k]), logged = FALSE)*pi_vec[k]
      } #end for r
    } # end for k
    denom_mat = t(repmat(apply(fp_mat,1,"sum"), n = K, m = 1))
    cprob_i = fp_mat/denom_mat
    return(cprob_i)

}
