#' Psuedo-class draw LC membership
#'
#' @param CPROBS_df (data.frame) N-by-K data frame with posterior class probabilities for each individual
#'
#' @return class_vec (vector) Vector with N psuedoclass draws for each individual
#' @export
#'
#' @examples
get_pseudoc_draw<-function(CPROBS_df){
  
# Input: 
#    - CPROBS_df (dataframe) N-by-K data frame with posterior class probabilities for each individual
# Output:
#    - class_vec (vector) Vector with N psuedoclass draws for each individual

class_vec = apply(CPROBS_df, 1, function (p_i) which(rmultinom(n=1, size = 1, prob = p_i)==1))

return(class_vec)

}