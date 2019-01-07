#' Apply missing data mechanism to simulate observed LPA indicator data
#'
#' This function accepts complete data and applies the missing data mechanism.
#' @param z  (integer) condition identifier (number) for the complete data (as specified by data_conditions)
#' @param df  (data.frame) "complete data" from get_complete_data() function
#' @param pctmiss_vec  (vector) percent missing simulation conditions
#' @param data_conditions  (data.frame) simulation conditions pertaining to the complete data
#' @param rep (integer) Replication number (defaults to 1).
#' @param save_it (logical) if TRUE, then it saves the data set (i) temp_wd_rep_vec & (ii) rep also must be specificed
#' @param temp_wd_rep_vec (character vector). Required if save_it==TRUE. Replication-specific temporary directory
#' @param pm_vec (integer) giving the correpsonding the percent missing condition. If NA, then defaults to 1:length(pctmiss_vec)
#' @return out_list  (list) with the following elements:
#'       (A) list_obsdf_z  (list) with length(pctmiss_vec) elements, each corresponding to a
#'                          (data.frame) of the observed with %missing values for the z-th condition number
#'       (B) dffolderfiles (data.frame) with the files and folders of the saved data
#' @export
#' @examples
#' get_obs_data(z,df,data_conditions, save_it = FALSE)

get_obs_data<-function(z, df,pctmiss_vec,data_conditions, rep = NA,
                       save_it = FALSE, temp_wd_rep_vec = NULL, pm_vec = NA){


  require(MASS)
  require(tidyverse)

  temp_wd_rep_z = paste0(temp_wd_rep_vec[rep],"/z",z)

  if (save_it==TRUE){require(MplusAutomation)}


  if(save_it == TRUE & (is.na(rep) | is.null(temp_wd_rep_z) | is.na(rep))){
    require(MplusAutomation)
    stop("rep, & temp_wd_rep_z arguments cannot be NA or NULL with save_it = TRUE")
  }

  if(length(pm_vec)==1){
    if(is.na(pm_vec)){
      pm_vec = 1:length(pctmiss_vec)
    }
  }


  if(length(pm_vec)!=length(pctmiss_vec)){
    stop("length(pm_vec)!=length(pctmiss_vec)")
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
  if(!is.null(ncol(Phi))){stop("Code not currently set up to handle J_Xcom_z>1")}
  rankPhi = rank(Phi, ties.method = "random")/N_z


  # For each % missing condition, identify observations that will have missing data.
  # NOTE: The way this is specified, it is the percent of observations missing at least ONE
  # data point.
  list_inds.miss = NULL
  for (i in seq(1,length(pctmiss_vec))){
    if (i==1){
      list_inds.miss[[i]] =  which(rankPhi<=pctmiss_vec[i])
    } else {
      list_inds.miss[[i]] = setdiff(which(rankPhi<=pctmiss_vec[i]), list_inds.miss[[i-1]])
    }
  }

  # Construct the missingness-type vector for each observation (R_df). Values mean:
  # Inf - Not missing under any % missing condition
  # 1 - missing under the first % missing condition
  # 2 - missing under the first AND second % missing condition
  # etc, etc.
  notmi<-function(x){sum(is.infinite(x))/length(x)}; nobs_fn<-function(x){sum(is.infinite(x))}
  midesc_list = list(rates = list(NULL), counts = list(NULL))

  R_df = data.frame(mat.or.vec(nr = N_z, nc = J_Y_z+J_Xinc_z)+Inf);
  names(R_df) = paste("R_",names(df)[seq(1,J_Y_z+J_Xinc_z,by = 1)], sep = "")

  for(i in seq(1,length(pctmiss_vec))){
    inds.miss = list_inds.miss[[i]]
    for (j in inds.miss){
      flag = 0
      while(flag == 0){
        njmiss = sample(seq(1,J_Y_z+J_Xinc_z-1, by = 1),1, FALSE)
        cols_miss = sort(sample(seq(1,J_Y_z+J_Xinc_z, by = 1),njmiss,FALSE))
        flag = ifelse(sum((1:J_Y_z)%in%cols_miss)==J_Y_z,0,1)
      }
      R_df[j,cols_miss] <- length(pctmiss_vec) - i  + 1
    }

    df_rates_z_pm = R_df%>%transform(subpop = df$subpop)%>%group_by(subpop)%>%summarise_all(funs(notmi))%>%data.frame()
    df_counts_z_pm = R_df%>%transform(subpop = df$subpop)%>%group_by(subpop)%>%summarise_all(funs(nobs_fn))%>%data.frame()
    names(df_rates_z_pm)<-names(df_counts_z_pm)<-c("subpop",names(df)[seq(1,J_Y_z+J_Xinc_z,by = 1)])
    df_rates_z_pm = df_rates_z_pm %>%transform(z=z,rep=rep,pm=pm_vec[i],pctmiss=pctmiss_vec[i])
    df_counts_z_pm = df_counts_z_pm %>%transform(z=z,rep=rep,pm=pm_vec[i],pctmiss=pctmiss_vec[i])
    midesc_list$rates[[i]] = df_rates_z_pm; midesc_list$counts[[i]] = df_counts_z_pm
    names(midesc_list$rates)[[i]]<-names(midesc_list$counts)[[i]]<-paste0("pm",pm_vec[i])
  }
  # Construct a list with the observed data for each missingness condition, with identified missing values from R_df set to missing
  # Save observed data, as required.
  out_list = list(list_obsdf = NULL, dffolderfiles = NULL, obs_desc = midesc_list)
  list_obsdf_z = NULL
  for(i in seq(1,length(pctmiss_vec))){

    obsdf_z = df
    for (j in seq(1,J_Y_z+J_Xinc_z,by = 1)){
      obsdf_z[R_df[,j]<=i,j] = NA
    }
    obsdf_z = transform(obsdf_z, pctmiss = pctmiss_vec[i])
    obsdf_z = cbind(obsdf_z, R_df)

    list_obsdf_z[[i]] = obsdf_z
    names(list_obsdf_z)[[i]] = paste0("pm",pm_vec[i])

    # Create the out_list for return and save, if needed

    if (save_it == TRUE){

      dat_folder =paste0("Observed data/pm",pm_vec[i])
      dat_fname = paste0("obsdf rep", rep , " z", z, " pm", pm_vec[i],".dat")

      out_list$dffolderfiles = rbind(out_list$dffolderfiles,
                                     data.frame(rep = rep, z = z,
                                                folders =  paste0("rep",rep,"/z",z,"/",dat_folder),
                                                files = dat_fname,
                                                pm = pm_vec[i], pva = NA, m = NA) )


      # Save a ".dat" file for Mplus
      nms = names(obsdf_z)
      jkeep = which(startsWith(nms,"Y") | startsWith(nms,"X") | nms=="subpop")
      prepareMplusData(obsdf_z, keepCols = jkeep,
                       filename = paste(temp_wd_rep_z,dat_folder,dat_fname, sep = "/"), inpfile = FALSE,
                       overwrite = TRUE)
    }# End save_it == TRUE

  }
  out_list$list_obsdf = list_obsdf_z


  return(out_list)

}
