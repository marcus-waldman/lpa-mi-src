#' Compile fitted models to an outlist for complete and observed data only
#'
#' @param out_ftc (list) from out_ftc()
#' @param out_switched (list) from resolve_label_switch
#' @param Data (character) indicating 
#' @param p (integer)
#' @param z (integer)
#' @param rep (integer)
#' @param pm (integer) defaults to NA
#' @param pva (integer) defaults to NA
#'
#' @return
#' @export
#'
#' @examples
#' 
to_outlist<-function(out_ftc, out_switched, Data, p, z, rep, pm = NA, pva = NA){
  
  outlist_complete = list(Data = Data,p = p, z = z, rep = rep, pm = NA, pva = NA, MM_pva = NA, 
                          converged = out_ftc$Converged, 
                          rcond = out_ftc$Rcond, 
                          starts = out_ftc$Starts,
                          problem = out_ftc$problem,
                          summary = out_switched$out_readModels$summaries,
                          parameters = out_switched$out_readModels$parameters$unstandardized, 
                          switch_list = list(class_switch = out_switched$class_switched, map_from = out_switched$map_from, map_to = out_switched$map_to, 
                                             input_from = out_ftc$out_readModels$input, input_to = out_switched$out_readModels$input),
                          imputation_list = list(NULL)
  )
  
  return(outlist_complete)
  
}