#' Identify label switching occured using KL divergence
#'
#' @param nn (integer) row in tracker_df
#' @param Qlist_x (list) List with sample population estimates for LPA model
#' @param Plist_x (list) List with population values for LPA model
#' @param name (character) name of discrete optimization solver
#' @param t_max (numeric) maximum time
#' @param approximate 0 or 1. Defaults to 0 indicating not ot use an approximate solution.
#' @param round_cplex See bmatch
#' @param trace_cplex See bmatch
#'
#' @return
#' @export
#'
#' @examples KL_check_switching(mu_est_mat, S_est_array, z,)
KL_check_switching<-function(nn, Qlist_x, Plist_x,
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


      require(gaussDiff)



      # Get number of class indicators and number of classes based on data condtion z
      J_Y = nrow(Plist_x$mu)
      K = dim(Plist_x$S)[3]

      # Identify the pairwise KL divergences across each population and sample/estimated classess
      DMAT = mat.or.vec(K,K) + NA
      for (k in 1:K){

        mu_k = Qlist_x$mu[1:J_Y,k]
        S_k = Qlist_x$S[1:J_Y,1:J_Y,k]

        kl_vec = rep(NA,K)

        for (kk in 1:K){

          kl_vec[kk] = normdiff(mu1 = mu_k, sigma1 = S_k,
                                mu2 = Plist_x$mu[1:J_Y, kk], sigma2 = Plist_x$S[1:J_Y,1:J_Y,kk],
                                method = c("KL"))
          #kl_vec[kk] = kl_vec[kk] +  normdiff(mu1 = Plist_x$mu[1:J_Y, kk], sigma1 = Plist_x$S[1:J_Y,1:J_Y,kk],
          #                                    mu2 = mu_k, sigma2 = S_k,
          #                                    method = c("KL"))
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
        require(designmatch)
         obj.bmatch = bmatch(t_ind = c(rep(1,K), rep(0,K)),dist_mat = t(DMAT),
                             solver = list(name = name, t_max = t_max, approximate = approximate, round_cplex = round_cplex,   trace_cplex = trace_cplex),
                             total_groups = K)
         map_to =  obj.bmatch$c_id-K
      }
      if(identical(as.integer(map_to),1:K)){class_switched = FALSE}

      # Create and return out list
      df_map_to = data.frame(t(map_to))
      names(df_map_to) = paste0("c",1:K)
      df_switched = data.frame(nn = nn, switched = class_switched, singular = FALSE)
      out_df = cbind(df_switched, df_map_to)


      return(out_df)

}
