#' Simulate complete LPA data for simulation study
#'
#' This function simulates complete data.
#' @param z (integer) condition identifier for the complete data (as specified by data_conditions)
#' @param data_conditions (data.frame) simulation conditions pertaining to the complete data
#' @param rep (integer) Replication number (defaults to 1).
#' @param save_it (logical) if TRUE, then it saves the data set and the following also must be specificed
#' @param temp_wd_rep_vec (character vector). Required if save_it==TRUE.  processor-specific temporary directory
#' @return out_list  (list) with the following elements
#'       (A) dfcom - (data.frame) complete data corresponding to the z-th condition number
#'       (B) mu - (J-by-K matrix) with the means for the j-th variable in the k-th class.
#'       (C) S - (J-by-J-by-K array) for the k-th class's covariance matrix
#'       (D) pi - (vector) with K elements. marginal probabilties for the k-th class
#'       (E) dffolderfiles (data.frame) with the files and folders of the saved data
#' @export
#' @examples
#' get_complete_data(z,data_conditions,save_it = FALSE)

get_complete_data<-function(z,data_conditions, rep = NA, save_it = FALSE, temp_wd_rep_vec = NULL){

#Stuff for error diagnosing:
# rep = 1
# z = 1
# save_it = TRUE



      temp_wd_rep_z = paste0(temp_wd_rep_vec[rep],"/z",z)

      require(MASS); require(psych)

      if (save_it == TRUE){require(MplusAutomation)}
      if(save_it == TRUE & (is.na(rep) | is.null(temp_wd_rep_z) | is.na(rep))){
        stop("rep & temp_wd_rep_z arguments cannot be NA or NULL with save_it = TRUE")
      }

      # Get population-level parameters for the mixture model
      out_get_FMM<-get_FMM_params(z=z,data_conditions=data_conditions)
      mu_z = out_get_FMM$mu_z; S_z = out_get_FMM$S_z; pi_z = out_get_FMM$pi_z; K_z = out_get_FMM$K_z
      J_Y_z = out_get_FMM$J_Y_z; J_Xinc_z = out_get_FMM$J_Xinc_z; J_Xcom_z = out_get_FMM$J_Xcom_z
      J = out_get_FMM$J; MD_z = out_get_FMM$MD_z;
      rho_YX_z = out_get_FMM$rho_YX_z; dX_z = out_get_FMM$dX_z

      # Generate latent classe memberships
      N_z = data_conditions$N[z]
      Nk_z = floor(pi_z*N_z)
      class_z = rep(NA,N_z)
      for(k in 1:K_z){
        lower = ifelse(k == 1, 1, (k-1)*Nk_z[1,k-1]+1)
        upper = lower + Nk_z[1,k] - 1
        class_z[lower:upper] = k
      }
      inds.na = which(is.na(class_z))
      if(!is.null(inds.na)){
        draws = rmultinom(length(inds.na),1,as.vector(pi_z))
        ii = 0
        for (i in inds.na){
          ii = ii+1
          class_z[i] = which(draws[,ii]==1)
        }
        rm(draws); rm(ii)
        class_z = sort(class_z)
      }


      # Generate the complete data
      dfcom_z = data.frame(mat.or.vec(nr = N_z, nc = J)+NA)
      names(dfcom_z)[1:J_Y_z] = paste("Y",1:J_Y_z, sep = "")
      if (J_Xinc_z>0){names(dfcom_z)[seq(J_Y_z+1,J_Y_z+J_Xinc_z, by = 1)]=paste("Xinc",1:J_Xinc_z,sep = "")}
      if (J_Xcom_z>0){names(dfcom_z)[seq(J_Y_z+J_Xinc_z+1, J_Y_z+J_Xinc_z+J_Xcom_z,by = 1)]= paste("Xcom",1:J_Xcom_z,sep = "")}
      dfcom_z = transform(dfcom_z, subpop = class_z)
      for (k in 1:K_z){
        inds_k = which(dfcom_z$subpop==k)
        n_k = length(inds_k)
        temp =  mvrnorm(n = n_k, mu = mu_z[,k], Sigma = S_z[,,k])
        dfcom_z[inds_k,1:J] = temp
      }
      dfcom_z = transform(dfcom_z, rep = rep, data_condition = z)

      # Do descriptive statistics
      univariates_z = describeBy(dfcom_z[,seq(1,J+1)], group = "subpop", fast = TRUE)
      bivariates_z = list(NULL)
      for (k in 1:K_z){
        tmp = subset(dfcom_z, subpop == k)
        bivariates_z[[k]] = cor(tmp[,1:J])
      }
      descriptives_z = list(univariates = univariates_z, bivariates = bivariates_z)

      # Create the out_list for return and save, if needed
      out_list = list(dfcom = dfcom_z, mu = mu_z, S = S_z, pi = pi_z, dffolderfiles = NULL,
                      descriptives = descriptives_z)

      if (save_it == TRUE){

        dat_folder ="Complete data"
        dat_fname = paste0("dfcom rep", rep , " z", z,".dat")

        out_list$dffolderfiles = data.frame(rep = rep, z = z,
                                            folders = paste0("rep",rep,"/z",z,"/",dat_folder),
                                            files = dat_fname,
                                            pm = NA, pva = NA, m = NA)


        # Save a ".dat" file for Mplus
        nms = names(dfcom_z)
        jkeep = which(startsWith(nms,"Y") | startsWith(nms,"X") | nms=="subpop")
        prepareMplusData(dfcom_z, keepCols = jkeep,
                         filename = paste(temp_wd_rep_z,dat_folder,dat_fname, sep = "/"), inpfile = FALSE,
                         overwrite = TRUE)
      }

      return(out_list)
}
