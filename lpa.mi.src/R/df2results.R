#' df2results
#'
#' What does this function do?
#' @param df
#' @param pop_params
#' @param dtype Defaults to ""
#' @param p Defaults to NA
#' @param rep Defaults to NA
#' @param z Defaults to NA
#' @param pm Defaults to NA
#' @param pva Defaults to NA
#' @return
#' @export
#' @examples
#' lavaan2resultsa(obj_lavaan, pop_params)

df2results<-function(df, pop_params, dtype = "", p = NA, rep = NA, z = NA, pm = NA, pva = NA){

  J_Y_z = pop_params$J_Y_z
  K_z = pop_params$K_z

  intercept = paste("Y", 1:J_Y_z, "~1", sep = "")
  vcovar = mat.or.vec(nr = J_Y_z, nc = J_Y_z)
  for (i in 1:J_Y_z){
    for (j in 1:i){
      vcovar[i,j] = paste("Y",j,"~~","Y",i, sep = "")
    }
  }
  covariances = vcovar[lower.tri(vcovar, diag = TRUE)]

  params = as.character(c(intercept,covariances))
  rep_df = data.frame(class = sort(rep(1:K_z, length(params))))
  rep_df = transform(rep_df, parameter = rep(params, K_z))
  rep_df$parameter = as.character(rep_df$parameter)

  rep_df = transform(rep_df, estimate = NA, se = NA, value = NA)
  for(k in unique(rep_df$class)){

    inds_k = which(df$class ==k)
    rep_df$estimate[which(endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = apply(df[inds_k,1:J_Y_z],2,mean, na.rm = TRUE)

    vc_k = cov(df[inds_k, 1:J_Y_z], use = "pairwise.complete.obs")
    rep_df$estimate[which(!endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = vc_k[lower.tri(vc_k, diag = TRUE)]


    mu_z = pop_params$mu_z
    rep_df$value[which(endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = pop_params$mu_z[1:J_Y_z,k]

    S_z = pop_params$S_z[1:J_Y_z,1:J_Y_z,k]
    rep_df$value[which(!endsWith(rep_df$parameter,"~1") & rep_df$class==k)] = S_z[lower.tri(S_z, diag = TRUE)]
  }

  rep_df = transform(rep_df, data_type = dtype, p = p, rep = rep, z = z, pm = pm, pva = pva)
  return(rep_df)

}

lavaan2results<-function(obj_lavaan, pop_params, dtype = "",
                         p = NA, rep = NA, z = NA,  pm = NA, pva = NA){
  J_Y_z = pop_params$J_Y_z
  summ_lavaan = summary(obj_lavaan);
  est_lavaan = summ_lavaan$est; se_lavaan = summ_lavaan$se;
  lhs_lavaan = summ_lavaan$lhs; op_lavaan = summ_lavaan$op;
  rhs_lavaan = summ_lavaan$rhs; gr_lavaan = summ_lavaan$group;
  rep_df = data.frame(
    class = summ_lavaan$group,
    parameter = paste(summ_lavaan$lhs, summ_lavaan$op, summ_lavaan$rhs, sep = ""),
    estimate = est_lavaan,
    se = se_lavaan
  )
  rep_df = transform(rep_df, value = NA)
  for(k in unique(rep_df$class)){

    mu_z = pop_params$mu_z
    rep_df$value[which(op_lavaan=="~1" & rep_df$class==k)] = pop_params$mu_z[1:J_Y_z,k]

    S_z = pop_params$S_z[1:J_Y_z,1:J_Y_z,k]
    var_k = diag(S_z)
    cov_k = as.vector(S_z[lower.tri(S_z)])
    rep_df$value[which(op_lavaan=="~~" & rep_df$class==k)] = c(var_k,cov_k)
  }
  rep_df = transform(rep_df, data_type = dtype, p = p, rep = rep, z = z, pm = pm, pva = pva)
  return(rep_df)
}

model_lavaan = '
# Means
Y1~1
Y2~1
Y3~1

# Variances
Y1~~Y1
Y2~~Y2
Y3~~Y3

# Covariance
Y1~~Y2
Y1~~Y3
Y2~~Y3
'
