#' Resolve label switching with population values known.
#'
#'
#' @param mu_est_mat (matrix) 
#' @param S_est_array (array)
#' @param z (integer)
#' @param data_conditions (data.frame)
#' @param parameters_df (data.frame) Defaults to Null. If non-Null, switches the labels in the data frame.
#' @param name (character) Deflau
#' @param t_max (integer)
#' @param approximate (integer)
#' @param round_cplex (integer)
#' @param trace_cplex (integer)
#' @return (list) A list with the following things: 
#'           (a) class_switched - (logical) identifies if class switching had occurred.
#'           (b) map_from - (vector) always just 1:K showing class definitions from sample
#'           (c) map_to - (NULL or vector) same as map_from if no class switching. 
#'                      Otherwise it is the corresponding population-based classes from sample classes
#'           (d) parameters_df (data.frame). Only included if parameters_df input is non-null. Same as input, except labels are switch. 
#'
#' @export
#' @examples
#'
#'
#'
#'
resolve_label_switch<-function(mu_est_mat, S_est_array, z, data_conditions,
                               parameters_df = NULL, name = "glpk", t_max = 30, approximate = 0, round_cplex = 0, trace_cplex = 0){
                                        
      require(gaussDiff); require(designmatch); require(plyr)
      
      # Get number of class indicators and number of classes based on data condtion z
      J_Y = data_conditions$J_Y[z]
      K = data_conditions$K[z]
      
      # Get population parameters based on the data condition z
      out_get_FMM<-get_FMM_params(z,data_conditions)    
  
  
      # Identify the pairwise KL divergences across each population and sample/estimated classess
      DMAT = mat.or.vec(K,K) + NA  
      for (g in 1:K){
        
        mu_g = mu_est_mat[1:J_Y,g]
        S_g = S_est_array[1:J_Y,1:J_Y,g]
        
        kl_vec = rep(NA,K)
        
        for (gg in 1:K){
          
          kl_vec[gg] = normdiff(mu1 = mu_g, sigma1 = S_g,
                                mu2 = out_get_FMM$mu_z[1:J_Y, gg], sigma2 = out_get_FMM$S_z[1:J_Y,1:J_Y,gg], 
                                method = c("KL"))
          kl_vec[gg] = kl_vec[gg] +  normdiff(mu1 = out_get_FMM$mu_z[1:J_Y, gg], sigma1 = out_get_FMM$S_z[1:J_Y,1:J_Y,gg],
                                              mu2 = mu_g, sigma2 = S_g,
                                              method = c("KL"))
        }
        
        DMAT[,g] = kl_vec
      }
      
      #Identify if class switching occured
      class_switched = FALSE;
      if(prod(apply(DMAT,2,"which.min") == c(1:K))==0){class_switched = TRUE}
      
      # If class switching occured, relable classes according to population labels
      map_from = 1:K
      map_to = 1:K
      if (class_switched){
         obj.bmatch = bmatch(t_ind = c(rep(1,K), rep(0,K)),dist_mat = t(DMAT), 
                             solver = list(name = name, t_max = t_max, approximate = approximate, round_cplex = round_cplex,   trace_cplex = trace_cplex), 
                             total_groups = K)
         map_to =  obj.bmatch$c_id-K
      }
      
      # Create and return out list
      out_list = list(class_switched = class_switched, 
                      map_from = map_from, map_to = map_to, converged = TRUE)
      
      
      if (is.null(parameters_df)==FALSE){
        parameters_df$LatentClass = mapvalues(parameters_df$LatentClass, from = map_from, to = map_to)
        out_list$parameters_df = parameters_df
      }
      
      return(out_list)
      
}
