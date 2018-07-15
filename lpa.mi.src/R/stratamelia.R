#' Stratified Amelia Imputation
#'
#' This function obtain imputations using Amelia, stratified by observed group membership
#' @param x (data.frame) Incomplete data set
#' @param strata  (vector of integers) Group membership indicator
#' @param ... additional arguments passed to the amelia function 
#' @return out_list (list) with each element being an imputed data set.
#' @export
#' @examples
#' stratamelia(x, strata, m = 50)


stratamelia<-function(x,strata,m = 100,...){

  require(Amelia)

  strata_vec = sort(unique(strata))
  
  out_list = rep(list(mat.or.vec(nr = nrow(x), nc = ncol(x))+NA),m)
  
  
  for (k in strata_vec){
    inds_k = which(strata == k)
    obj_amelia = do.call(amelia, args = list(x = x, m = m, ...))
    #head(obj_amelia$m)
    
    for (mm in seq(1,m)){
      out_list[[mm]][inds_k, ] = obj_amelia$imputations[[mm]]
    }
    
  }
  
  return(out_list)
  
  
}