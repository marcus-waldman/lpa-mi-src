#' Calculate posterior class probabilities given a latent class regression
#'
#' This function calculates the posterior class probabiltieis from a finite mixture of multivariate normal distributions.
#' @param Y_i (matrix) of size NxJ of observed values
#' @param X_i (matrix) Nx(P+1) design matrix. First column all 1's
#' @param pi_vec (numeric vector) with K elements containing component probabilities
#' @param beta_array (array) of size Jx(P+1)xK, where J is number of items, P is the number of predictors
#' @param S_array (array) of size JxJxK with within class residual covariance matrix.
#' @param list_mdpattern (list) Optional. Generated from get_mdpatterns function.
#' @param fast (logical) If TRUE, the use Rfast::dmvnorm. Otherwise use mixtools::dmvnorm
#' @return (matrix) of size NxK containing the posterior class porbabilties
#' @export
#' @examples
#' cprobs_lcr(Y_i, X_i pi_vec, beta_array, S_array)


cprobs_lcr<-function(Y_i, X_i, pi_vec, beta_array, S_array, list_mdpattern = NULL, fast = FALSE){

    require(Rfast)
    require(pracma)
    require(dynr)
    require(mixtools)
    require(pkgmaker)

    K = length(pi_vec)
    J = ncol(Y_i)
    N = nrow(Y_i)

    if(unique(X_i[,1])!=1){stop("First column in X_i must be all ones for the intercept.")}

    if (is.null(list_mdpattern)){list_mdpattern = get_mdpatterns(df = as.data.frame(Y_i))}

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



        for(i in ix_r){

          if (fast==TRUE){
          fp_mat[i, k] = Rfast::dmvnorm(x = Y_r[i,],
                                    mu = as.matrix(X_i[i, ]%*%t(beta_array[idV, ,k])),
                                    sigma = S_rk, logged = FALSE)*pi_vec[k]
          
          } else {
          fp_mat[i, k] = mixtools::dmvnorm(y = Y_r[i,],
                                              mu = as.matrix(X_i[i, ]%*%t(beta_array[idV, ,k])),
                                              sigma = S_rk)*pi_vec[k]
          }
        } #end for i


      } # end for r
    } # end for k
    denom_mat = t(repmat(apply(fp_mat,1,"sum"), n = K, m = 1))
    cprob_i = fp_mat/denom_mat
    return(cprob_i)

}
