#' Stratified Amelia Imputation
#'
#' This function obtain imputations using Amelia, stratified by observed group membership
#' @param x (data.frame) Incomplete data set
#' @param strata  (vector of integers) Group membership indicator
#' @param ... additional arguments passed to the amelia function 
#' @return out_list (list) with each element being an imputed data set.
#' @export
#' @examples
#' stratamelia(x, strata, m = 1)


stratamelia<-function(x,strata,m,...){

  require(Amelia)

  
  imparray = array(NA, dim = c(nrow(x), ncol(x), m))
  
  for (k in unique(strata)){
    inds_k = which(strata == k)
    x_k = x[inds_k, ]
    args_k = list(x = x_k, m = m, ...)
#args_k = list(x = x_k, m = m)
    obj_k = do.call("amelia", args_k)
    imparray[inds_k,,] = array(unlist(obj_k$imputation), dim = c(nrow(x_k), ncol(x_k), m))
  }
  longimp_df = data.frame(x)
  longimp_df = transform(longimp_df, subpop = strata, .imp = 0)
  for (mm in 1:m){
    tmp_df = data.frame(imparray[,,mm])
    names(tmp_df) = names(x)
    tmp_df = transform(tmp_df, subpop = strata, .imp = mm)
    longimp_df = rbind(longimp_df, tmp_df)
  }
  mids_obj = as.mids(longimp_df)
  
  return(mids_obj)

}