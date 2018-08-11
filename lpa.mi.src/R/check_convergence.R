#' Check convergence
#'
#' 
#' @param file (integer) 
#' @param folder
#' @param starts_txt
#' @return (logical) 
#' @export
#' @examples


check_convergence<-function(file, folder_wd = getwd(), starts_txt = "0;"){

#folder_wd = "C:/Users/marcu/Dropbox/Dissertation Proposal/Paper 1/Simulation Code/Impute LPA Sim - known K -  v1_0/Error diagnosis/p1/Complete data"
#file = "Naive dfcom p1 z1 rep15 v2.out"
  
    setwd(folder_wd)
    
    hi = readLines(con = file)
    i_0 = which(hi ==  "THE MODEL ESTIMATION TERMINATED NORMALLY" )
    i_1 = which(hi ==  "THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED.  RERUN WITH AT LEAST TWICE THE" )
   
    Converge = FALSE
    if (length(i_0)>0){
      if(starts_txt=="0;"){
        
        Converge = TRUE
        
      } else{
        
        if(length(i_1)==1){
          
          Converge = TRUE
          
          }
        
      }
      
    }

    return(Converge)
}
