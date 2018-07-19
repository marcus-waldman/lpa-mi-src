#' Create the within 
#'
#' This function obtain imputations using Amelia, stratified by observed group membership
#' @param mu_k - (vector) of length J_Y_z (where J_Y_z is the number of indicators)
#' @param S_k - (Matrix) of size J_Y_z by J_Y_z
#' @param usev_vec (vector) of charcters with the names of the latent class indicators
#' @param Model_Covariance - (logical) independence model or not. Defaults to FALSE.
#' @param symbol_txt ("character") Either "*" or the "(at)" symbol. The latter for fixing parameters.
#' @return Returns a character vector with Mplus code.
#' @export
#' @examples
#' create_withinclass_Mplus_syntax(mu_k, S_k, usev_vec)

naive_withinclass_Mplus_syntax<-function(mu_k, S_k, usev_vec, Model_Covariance = FALSE, symbol_txt = "*"){
  
  # Syntax for the means with starting values
    means_txt = paste("[", usev_vec[1],symbol_txt, mu_k[1],"];", sep = "")
    for (t in seq(2,length(usev_vec))){
      means_txt = c(means_txt, 
                    paste("[", usev_vec[t],symbol_txt, mu_k[t],"];", sep = ""))
    }
    
  # Syntax for the variances with starting values
    vars_txt = paste("Y1*", S_k[1,1], ";", sep = "")
    for (t in seq(2,length(usev_vec))){
      vars_txt = c(vars_txt, 
                   vars_txt = paste("Y",t,symbol_txt, S_k[t,t], ";", sep = "")
      )
    }
    
  # Syntax for the covariances with starting values or held constant, as appropriate
    symbol = ifelse(Model_Covariance, symbol_txt,"@")
    covs_txt = NULL
    for(t1 in seq(1,length(usev_vec)-1)){
      for (t2 in seq(t1+1, length(usev_vec))){
        covs_txt = c(covs_txt, paste(usev_vec[t1], " WITH ", usev_vec[t2], symbol, S_k[t1,t2],";",sep = ""))
      }
    }
    
  # Output the within class Mplus modeling syntax
    within_class_txt = c(means_txt, 
                         vars_txt,
                         covs_txt)
    
    return(within_class_txt)
}