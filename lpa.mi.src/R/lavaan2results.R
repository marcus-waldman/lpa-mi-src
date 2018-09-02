#' lavaan2results
#'
#' What does this function do?
#' @param obj_lavaan
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
