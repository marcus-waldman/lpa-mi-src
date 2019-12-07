#' Format Mplus Parameters to Q list
#'
#'
#' @param params_df (data.frame) Unstandardized parameters extracted from MplusAutomation
#' @return (list) with pi, mu, and S as arguments defining the mixture model
#' @export
#' @examples
#' Mplus2Q(params_df)

Mplus2Qlist<-function(params_df){

  require(stringr)
    params_df$param = as.character(params_df$param)
    params_df$paramHeader = as.character(params_df$paramHeader)
    params_df$LatentClass = as.character(params_df$LatentClass)
    if("value" %in% names(params_df)){
      if ("est" %in% names(params_df)){stop("Both est and value cannot be supplied")}
      names(params_df)[names(params_df)=="value"] = "est"
    }
    if (!("est" %in% names(params_df))){stop("est or value not found but must be supplied")}

    K = max(as.numeric(subset(params_df, LatentClass !=  "Categorical.Latent.Variables")$LatentClass))
    J = length(subset(params_df, paramHeader == "Means" & LatentClass == "1")$est)
    Q = list(pi = NULL, mu = NULL, S = NULL)

    names_vec = subset(params_df, paramHeader == "Means" & LatentClass == "1")$param


    # Get mixing proportions
    if(K>1){
      gamma_vec = subset(params_df, paramHeader == "Means" & startsWith(LatentClass, "Categorical"))$est
      pi_vec = gamma2pi(gamma_vec)
      Q$pi = pi_vec
    } else {
      Q$pi = 1
    }


    # Get means
    mu_mat = matrix(NA, nrow = J, nc = K)
    for(k in 1:K){
      mu_mat[1:J,k] = subset(params_df, paramHeader == "Means" & LatentClass == as.character(k))$est
    }
    Q$mu = mu_mat

    # Get Sigma
    S_array = array(dim = c(J,J,K))
    for (k in 1:K){
      for (j in 1:J){
        inds_jk = with(params_df, which(paramHeader == "Variances" & LatentClass == k & param == names_vec[j]))
        S_array[j,j,k] = params_df$est[inds_jk]
      }
    }

    idx_cov = which(str_detect(string = params_df$paramHeader, ".WITH"))
    if (length(idx_cov)>0){
      cov_df = params_df[idx_cov, ]
      var_x = as.character(str_split(cov_df$paramHeader, pattern = ".WITH", simplify = TRUE)[,1])
      var_y = cov_df$param
      cov_df = transform(cov_df, k = as.numeric(LatentClass), j1 = match(var_x, names_vec), j2 = match(var_y, names_vec))
      for(ii in 1:nrow(cov_df)){
        S_array[cov_df$j1[ii], cov_df$j2[ii], cov_df$k[ii]] <- S_array[cov_df$j2[ii], cov_df$j1[ii], cov_df$k[ii]]<- cov_df$est[ii]
      }
    }
    S_array[is.na(S_array)] = 0
    Q$S = S_array

    return(Q)

}
