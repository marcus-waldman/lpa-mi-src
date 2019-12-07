#' Extract class proportions from Bayesian fitted LC model in Mplus .out file.
#'
#' Extract class proportions from Mplus .out file.
#' @param file (character) The filename of the Mplus .out file where the starting values are located
#' @param path (character) The path where the file is located
#' @return (character) vector with Model sytax at starting values from the .out file.
#' @export
#' @examples
#' extract_class_proportions(file = "example.out")

extract_class_proportions <- function(file, path=getwd()){
  
  require(stringr)
  require(tidyverse)
  
    hi = readLines(con = paste(path,"/",file, sep = ""))
  
    i_o = which(hi == "Class Proportions")+2
    i_f = which(hi == "MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES")-3
    
    bye = hi[i_o:i_f]
    list_split = str_split(str_squish(bye), " ")
    
    K = length(list_split)
    
    proportion_df = data.frame(paramHeader = rep("Proportion",K), 
                               param = rep("Class",K)) %>% 
                    transform(est = NA, posterior_sd = NA, pval = NA, lower_2.5ci = NA, upper_2.5ci = NA, sig = NA, LatentClass = 1:K)
    
    for(kk in 1:K){
      for (jj in 3:7)
      proportion_df[kk,jj] = list_split[[kk]][jj]
    }
    proportion_df$est = as.numeric(proportion_df$est)
    proportion_df$posterior_sd = as.numeric(proportion_df$posterior_sd)
    proportion_df$pval = as.numeric(proportion_df$pval)
    proportion_df$lower_2.5ci = as.numeric(proportion_df$lower_2.5ci)
    proportion_df$upper_2.5ci = as.numeric(proportion_df$upper_2.5ci)
    proportion_df$sig = as.numeric(proportion_df$sig)
    
    proportion_df$LatentClass = as.character(proportion_df$LatentClass)
    
    proportion_df = transform(proportion_df, sig = ifelse(pval<0.05,TRUE, FALSE))
    
    intercepts_df = data.frame(paramHeader = rep("Means",K-1), 
                               param = paste0("C#",seq(1,K-1))) %>%
                    transform(est = log(proportion_df$est[-K]/proportion_df$est[K]), 
                              posterior_sd = NA, 
                              pval = NA, 
                              lower_2.5ci = NA, 
                              upper_2.5ci = NA, 
                              sig = NA, 
                              LatentClass = "Categorical.Latent.Variables")
    
    out_df = rbind(proportion_df, intercepts_df)
    
    return(out_df)

}
