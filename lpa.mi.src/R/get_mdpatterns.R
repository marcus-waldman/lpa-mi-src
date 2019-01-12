#' Missing data patterns
#'
#' @param df (data.frame) of size N-by-J with missing values
#'
#' @return List with the following elements 
#'          (a) pattern_i - (data.frame) giving the missing data binary pattern (e.g., Rpatt) and numeric pattern identifier
#'          (b) table_r - (data.frame) providing missing data patterns present in the provided data frame
#'          (c) U_r - (list) A list of size equal to the number of patterns in table_r (minus 1 due to excluding complete data). Each 
#'                          elements contains the indices corresponding to missing values for a part
#'          (d) V_r = (list) Similar to U_r, except with observed data indices  
#' @export
#'
#' @examples
get_mdpatterns<-function(df){
  Rdf = with(df, data.frame(id = 1:nrow(df), 
                            Rpatt =apply(df, 1, FUN = function(x){paste0(as.numeric(!is.na(x)), collapse = "")})
                            ))
  tableRdf = data.frame(R = sort(unique(Rdf$Rpatt), decreasing = TRUE)) %>% transform(idR = seq(0,length(R)-1)) 
  row.names(tableRdf) <- with(tableRdf, as.character(idR))
  list_idU = lapply(2:nrow(tableRdf), function(x){simplify(gregexpr("0",tableRdf$R[x]))}) 
  list_idV = lapply(2:nrow(tableRdf), function(x){simplify(gregexpr("1",tableRdf$R[x]))})  
  
  
  Rdf = transform(Rdf, idR = with(tableRdf,idR[match(Rdf$Rpatt, R)]))
  
  return(list(patterns_i = Rdf, table_r = tableRdf, U_r = list_idU, V_r = list_idV))
}