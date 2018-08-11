resolve_class_switch_bmatch_ver2<-function(mu_est_mat, S_est_array, z, data_conditions,
                                      name = "glpk", t_max = 30, approximate = 0, round_cplex = 0, trace_cplex = 0){

  # Inputs:
  #     - mu_est_mat
  #     - S_est_array
  #     - z: (integer) The idnentifier for the data conditions
  #     - data_conditions: (data.frame)
  #     - name, t_max, approximate... all "solver" options for bmatch function
  # Outputs: 
  #     - A list with the following things: 
  #           (a) class_switched - (logical) identifies if class switching had occurred.
  #           (b) map_from - (vector) always just 1:K showing class definitions from sample
  #           (c) map_to - (NULL or vector) same as map_from if no class switching. 
  #                       Otherwise it is the corresponding population-based classes from sample classes
  
                                        
      require(gaussDiff); require(designmatch)
      
      # Get number of class indicators and number of classes based on data condtion z
      J_Y = data_conditions$J_Y[z]
      K = data_conditions$K[z]
      
      # Get population parameters based on the data condition z
      out_get_FMM<-get_FMM_params(z,data_conditions)    
  
  
      # Identify the pairwise KL divergences across each population and sample/estimated classess
      DMAT = mat.or.vec(K,K) + NA  
      for (k in 1:K){
        
        mu_k = mu_est_mat[1:J_Y,k]
        S_k = S_est_array[1:J_Y,1:J_Y,k]
        
        kl_vec = rep(NA,K)
        
        for (kk in 1:K){
          
          kl_vec[kk] = normdiff(mu1 = mu_k, sigma1 = S_k,
                                mu2 = out_get_FMM$mu_z[1:J_Y, kk], sigma2 = out_get_FMM$S_z[1:J_Y,1:J_Y,kk], 
                                method = c("KL"))
          kl_vec[kk] = kl_vec[kk] +  normdiff(mu1 = out_get_FMM$mu_z[1:J_Y, kk], sigma1 = out_get_FMM$S_z[1:J_Y,1:J_Y,kk],
                                              mu2 = mu_k, sigma2 = S_k,
                                              method = c("KL"))
        }
        
        DMAT[,k] = kl_vec
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
      return(out_list)
      
}
