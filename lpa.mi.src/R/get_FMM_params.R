#'Population parameters for the LPA/FMM model.
#'
#' Obtain the population parameters given the data
#' @param z  (integer) condition identifier (number) for the complete data (as specified by data_conditions)
#' @param data_conditions  (data.frame) simulation conditions pertaining to the complete data
#' @param t_rotate (numeric) Rotation angle in radians. Defaults to no rotation (i.e. t_rotate = 0)
#' @return out_list  (list) with the following elements: 
#'        (A) mu_z - (J-by-K matrix) with the means for the j-th variable in the k-th class.
#'        (B) S_z -  (J-by-J-by-K array) for the k-th class's covariance matrix
#'        (C) pi_z - (vector with K elements) marginal probabilties for the k-th class
#'        (D) K_z - (integer) specifying the number of classes in z-th data condition
#'        (E) J_Y_z - - (integer) specifying the number of latent class indicators in z-th data condition
#'        (F) J_Xcom_z - - (integer) specifying the number complete data missing data corrleates in z-th data condition
#'        (G) J_Xinc_z - - (integer) specifying the number incomplete data missing data correlates in z-th data condition
#'        (H) J - - (integer) specifying the total number of variables in the joint distribution in z-th data condition
#'        (I) MD_z - (numeric) specifying the class separation in z-th data condition
#'        (J) rho_YX_z - - (numeric) specifying the missing data correlates correlation in z-th data condition
#' @export
#' @examples
#' get_FMM_params(z,data_conditions)

get_FMM_params<-function(z,data_conditions, t_rotate = 0){
  

  # Last revised: 12/16/2017
  #
  # Inputs: 
  #     z - (integer) condition identifier (number) for the complete data (as specified by data_conditions)
  #     data_conditions - (data.frame) simulation conditions pertaining to the complete data
  # Outputs:
  #     out_list with the following elements:
  #        A) mu_z - (J-by-K matrix) with the means for the j-th variable in the k-th class.
  #        B) S_z -  (J-by-J-by-K array) for the k-th class's covariance matrix
  #        C) pi_z - (vector with K elements) marginal probabilties for the k-th class
  #        D) K_z - (integer) specifying the number of classes in z-th data condition
  #        E) J_Y_z - - (integer) specifying the number of latent class indicators in z-th data condition
  #        F) J_Xcom_z - - (integer) specifying the number complete data missing data corrleates in z-th data condition
  #        G) J_Xinc_z - - (integer) specifying the number incomplete data missing data correlates in z-th data condition
  #        H) J - - (integer) specifying the total number of variables in the joint distribution in z-th data condition
  #        I) MD_z - (numeric) specifying the class separation in z-th data condition
  #        J) rho_YX_z - - (numeric) specifying the missing data correlates correlation in z-th data condition


    ##### Needed Information from the data_conditions data frame #####
    
      # Get the number of mixture components
      K_z = data_conditions$K[z]
      
      # Get the number of indicators classes
      J_Y_z = data_conditions$J_Y[z]
      #if (J_Y_z!=2){
      #  stop(paste("Currently there must be two latent class indicators. Data condition ", z))
      #}
      
      J_Xcom_z = data_conditions$J_Xcom[z]
      J_Xinc_z = data_conditions$J_Xinc[z]
      J = J_Y_z + J_Xcom_z + J_Xinc_z
      
      # Get the marginal probabilities
      pi_z = data_conditions[z, startsWith(names(data_conditions),"pi")]; pi_z = pi_z[,1:K_z]
      if(sum(pi_z)!=1){
        stop(paste("Marginal probabilities do not sum to 1 for data condition", z))
      }
      
      # Get the Mahalanobis distances
      MD_z = data_conditions$MD[z]
      
      # Get the with class covariance matrix
      rho_YX_z = data_conditions$rho_YX[z]
      C_modifies_YX_z = data_conditions$C_modifies_YX[z]
      if (C_modifies_YX_z==TRUE){ kvec_rho_YX=seq(-1*rho_YX_z,rho_YX_z,len = K_z) }
      if (C_modifies_YX_z==FALSE){ kvec_rho_YX=rep(rho_YX_z, K_z) }
      
    
    
    #### Construct the population parameters using the information in the previous section #####
      # Mean vectors of the latent class indicators 
      mu_all = get_mu_all(J = J_Y_z, t_rotate = t_rotate)
      mu_Y_z = MD_z*mu_all[,round(seq(1,ncol(mu_all), len = K_z),0)]
      
      # Mean vecotrs of the covariates
      mu_X_z = mat.or.vec(nr = J_Xcom_z + J_Xinc_z, nc = K_z)
      
      #Combine the covariate and indicator column vectors
      mu_z = rbind(mu_Y_z, mu_X_z)
      
      # Calculate the class covariance matrices
      S_z = array(0, dim = c(J, J, K_z))
      for(k in 1:K_z){
        S_z[1:J,1:J,k] = diag(J)
        S_z[1:J_Y_z,-seq(1,J_Y_z),k] <- S_z[-seq(1,J_Y_z),1:J_Y_z,k] <- kvec_rho_YX[k]
      }
  
    
      
      
      
       # Create the out list
      out_list = list(mu_z = mu_z, 
                      S_z = S_z,
                      pi_z = pi_z,
                      K_z = K_z,
                      J_Y_z = J_Y_z,
                      J_Xinc_z = J_Xinc_z,
                      J_Xcom_z = J_Xcom_z,
                      J = J,
                      MD_z = MD_z,
                      rho_YX_z = rho_YX_z)
      
      return(out_list)
      
}