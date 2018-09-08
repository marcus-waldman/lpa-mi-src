#' Apply missing data mechanism to simulate observed LPA indicator data
#'
#' This function accepts complete data and applies the missing data mechanism.
#' @param z  (integer) condition identifier (number) for the complete data (as specified by data_conditions)
#' @param df  (data.frame) "complete data" from get_complete_data() function
#' @param pctmiss_vec  (vector) percent missing simulation conditions
#' @param data_conditions  (data.frame) simulation conditions pertaining to the complete data
#' @param rep (integer) Replication number (defaults to 1).
#' @param p (integer) Processor number (defaults to 1).
#' @param save_it (logical) if TRUE, then it saves the data set and the following also must be specificed
#' @param temp_wd_p_vec (character vector). Required if save_it==TRUE.  processor-specific temporary directory
#' @return out_list  (list) with the following elements:
#'       (A) list_obsdf_z  (list) with length(pctmiss_vec) elements, each corresponding to a
#'                          (data.frame) of the observed with %missing values for the z-th condition number
#'       (B) dffolderfiles (data.frame) with the files and folders of the saved data
#' @export
#' @examples
#' get_obs_data(z,df,data_conditions, save_it = FALSE)

get_obs_data<-function(z, df,pctmiss_vec,data_conditions, rep = NA, p = NA,
                       save_it = FALSE, temp_wd_p_vec = NULL){

  # Last revised: 01/29/2018
  #
  # Inputs:
  #     z - (integer) condition identifier (number) for the complete data (as specified by data_conditions)
  #     df - (data.frame) "complete data" from get_complete_data.R function
  #     pctmiss_vec - (vector) percent missing simulation conditions
  #     data_conditions - (data.frame) simulation conditions pertaining to the complete data
  #     save_it (logical) - if TRUE, then it saves the data set and the following also must be specificed
  #       A) rep - (integer) replication number.
  #       B) p - (integer) processor number
  #       C) temp_wd_p - (character) processor-specific temporary directory
  # Output:
  #     out_list - (list) with the following elements
  #       A)  list_obsdf_z - (list) with length(pctmiss_vec) elements, each corresponding to a
  #                          (data.frame) of the observed with %missing values for the z-th condition number
  #       B) dffolderfiles (data.frame) with the files and folders of the saved data

  temp_wd_p =temp_wd_p_vec[p]

  if (save_it==TRUE){require(MplusAutomation)}


  if(save_it == TRUE & (is.na(rep) | is.null(temp_wd_p) | is.na(p))){
    require(MplusAutomation)
    stop("rep, p, & temp_wd_p arguments cannot be NA or NULL with save_it = TRUE")
  }

# Error diagnostic stuff
#df = list_get_complete$dfcom

  # Pre-processing: Extract information and variables from data_conditions
  J_Y_z = data_conditions$J_Y[z]
  J_Xinc_z = data_conditions$J_Xinc[z]
  J_Xcom_z = data_conditions$J_Xcom[z]
  J = J_Y_z + J_Xinc_z + J_Xcom_z
  N_z = data_conditions$N[z]

  # Calculate an underlying latent variable (Phi) that determines missingness based on rank orderings (rankPhi)
  Phi = as.vector(df[,startsWith(names(df),"Xcom")] +
                    mvrnorm(n = N_z, mu = rep(0,J_Xcom_z), Sigma = 0.1*diag(J_Xcom_z)))
  rankPhi = rank(Phi, ties.method = "random")/N_z


  # For each % missing condition, identify observations that will have missing data.
  # NOTE: The way this is specified, it is the percent of observations missing at least ONE
  # data point.
  list_inds.miss = NULL
  for (pm in seq(1,length(pctmiss_vec))){
    if (pm==1){
      list_inds.miss[[pm]] =  which(rankPhi<=pctmiss_vec[pm])
    } else {
      list_inds.miss[[pm]] = setdiff(which(rankPhi<=pctmiss_vec[pm]), list_inds.miss[[pm-1]])
    }

  }

  # Construct the missingness-type vector for each observation (R_df). Values mean:
  # Inf - Not missing under any % missing condition
  # 1 - missing under the first % missing condition
  # 2 - missing under the first AND second % missing condition
  # etc, etc.
  R_df = data.frame(mat.or.vec(nr = N_z, nc = J_Y_z+J_Xinc_z)+Inf);
  names(R_df) = paste("R_",names(df)[seq(1,J_Y_z+J_Xinc_z,by = 1)], sep = "")
  for(pm in seq(1,length(pctmiss_vec))){
    inds.miss = list_inds.miss[[pm]]
    for (i in inds.miss){
      flag = 0
      while(flag == 0){
        njmiss = sample(seq(1,J_Y_z+J_Xinc_z-1, by = 1),1, FALSE)
        cols_miss = sort(sample(seq(1,J_Y_z+J_Xinc_z, by = 1),njmiss,FALSE))
        flag = ifelse(sum((1:J_Y_z)%in%cols_miss)==J_Y_z,0,1)
      }
      R_df[i,cols_miss] <- length(pctmiss_vec) - pm  + 1
    }
  }


  # Construct a list with the observed data for each missingness condition, with identified missing values from R_df set to missing
  # Save observed data, as required.
  out_list = list(list_obsdf = NULL, dffolderfiles = NULL)
  list_obsdf_z = NULL
  for(pm in seq(1,length(pctmiss_vec))){

    obsdf_z = df
    for (j in seq(1,J_Y_z+J_Xinc_z,by = 1)){
      obsdf_z[R_df[,j]<=pm,j] = NA
    }
    obsdf_z = transform(obsdf_z, pctmiss = pctmiss_vec[pm])
    obsdf_z = cbind(obsdf_z, R_df)

    list_obsdf_z[[pm]] = obsdf_z

    # Create the out_list for return and save, if needed

    if (save_it == TRUE){

      out_list$dffolderfiles = rbind(out_list$dffolderfiles,
                                     data.frame(folders = paste0("Observed data/pm",pm),
                                                files = paste("obsdf p", p," z",z," rep",rep, " pm", pm,".dat",sep = ""),
                                                data_condition = z, m = NA) )

      # Save a copy in RData format
      # save(obsdf_z,
      #      file = paste0(temp_wd_p,"/Observed data/pm",pm,"/obsdf p", p, " z", z," rep", rep, " pm", pm,".RData"))

      # Save a ".dat" file for Mplus
      nms = names(obsdf_z)
      jkeep = which(startsWith(nms,"Y") | startsWith(nms,"X") | nms=="subpop")
      prepareMplusData(obsdf_z, keepCols = jkeep,
                       filename = paste0(temp_wd_p,"/Observed data/pm",pm,"/obsdf p", p, " z", z, " rep", rep, " pm", pm,".dat"), inpfile = FALSE,
                       overwrite = TRUE)
    }# End save_it == TRUE

  }
  out_list$list_obsdf = list_obsdf_z


  return(out_list)

}
