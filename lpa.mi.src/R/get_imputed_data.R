#' Obtain (multiply) imputed datasets
#'
#' Conduct (multiple) imputation on the observed data for the LPA-MI simulation study.
#' @param z (integer) condition identifier for the complete data (as specified by data_conditions)
#' @param list_get_obs (list) List obtained from the get_obs_data() function. Contains the following:
#'      A) dfcom - (data.frame) complete data corresponding to the z-th condition number
#       B) mu - (J-by-K matrix) with the means for the j-th variable in the k-th class.
#       C) S (J-by-J-by-K array) for the k-th class's covariance matrix
#       D) pi (vector with K elements) marginal probabilties for the k-th class
#       E) dffolderfiles (data.frame) with the files and folders of the saved data
#' @param list_get_complete (list) List obtained from the get_complete_data() function. Contains the following:
#'      A) dfcom - (data.frame) complete data corresponding to the z-th condition number
#       B) mu - (J-by-K matrix) with the means for the j-th variable in the k-th class.
#       C) S (J-by-J-by-K array) for the k-th class's covariance matrix
#       D) pi (vector with K elements) marginal probabilties for the k-th class
#       E) dffolderfiles (data.frame) with the files and folders of the saved data
#' @param methods_list (list) List w/ pva elements. Designates the imputation methods to try.
#' @param list_lcfcs_svals (list) Starting values for the LC-FCS proposed method.
#' @param data_conditions (data.frame) simulation conditions pertaining to the complete data
#' @param pctmiss_vec  (vector) Vector w/ pm elements. Designates percent missing simulation conditions
#' @param rep (integer) Replication number (defaults to 1).
#' @param save_it (logical) if TRUE, then it saves the data set and the following also must be specificed
#' @param temp_wd_rep_vec (character vector). Required if save_it==TRUE. Replication-specific temporary directory
#' @param pm_vec (integer) giving the correpsonding the percent missing condition. If NA, then defaults to 1:length(pctmiss_vec)
#' @param ... Additional arguments to the imputation function
#' @return  out_list - (list) with the following elements
#'       (A) obj_call - (mids) named mids object for imputed data each
#'           (i)  pm-th percent missing
#'           (ii) pva-th imputation method (see methods_list for details)
#'       (B) dffolderfiles (data.frame) with the files and folders for each of the save imputed data sets
#' @export
#' @examples
#' get_imputed_data(z,list_get_obs, list_get_complete, methods_list, data_conditions, save_it = FALSE)

get_imputed_data<-function(z, list_get_obs, list_get_complete, methods_list, data_conditions, pctmiss_vec,
                           rep = NA, save_it = FALSE, temp_wd_rep_vec = NULL, pm_vec = NA, ...){

  require(tidyverse)
  require(Amelia)
  require(fastDummies)
  require(mice)
  require(miceadds)
  require(mitools)
  if (save_it==TRUE){require(MplusAutomation)}

  # Extract imputation algorithm identifier
  pva_vec = methods_list$pva_vec #imputation algorithm identifier
  # Assign percent missing identifier if not specified
  if(length(pm_vec)==1){
    if(is.na(pm_vec)){
      pm_vec = 1:length(pctmiss_vec)
    }
  }

  # Check that pva_vec is the same as the number of elements
  if(!identical(length(methods_list$procedure), length(methods_list$pva_vec))){
    stop("pva_vec does not equal elements in methods_list.")
  }

  # Check that user specifed pct missing identifier with same number of elements as vector of proportions
  if(length(pm_vec)!=length(pctmiss_vec)){
    stop("length(pm_vec)!=length(pctmiss_vec)")
  }

  # Check if that an exported .dat file can be saved to a folder b/c such a folder is specified

  if(save_it == TRUE){
    if(is.na(rep) | is.null(temp_wd_rep_vec)){
      stop("rep & temp_wd_rep_vec arguments cannot be NULL with save_it = TRUE")
    } else {
      temp_wd_rep_z = paste0(temp_wd_rep_vec[rep],"/z",z)
    }
  }

  # Pre-processing, get needed variables for later
  J_Y = data_conditions$J_Y[z]; J_Xcom = data_conditions$J_Xcom[z]; J_Xinc = data_conditions$J_Xinc[z]
  J = J_Y + J_Xcom + J_Xinc

  # Allocate space for the out_list
  out_list = list(obj_call=NULL, dffolderfiles = NULL)
  tmp_list = NULL
  tmpj_list = list(NULL)
  # Create container or tmp_list
  for(jj in 1:length(pva_vec)){
    tmpj_list[[jj]] = NA
    names(tmpj_list)[jj] = paste0("pva",pva_vec[jj])
  }

  #ii = 1
  #jj = 1
  # Loop over A) pm =  1:length(pctmiss_vec), 2) pva = 1:Nprocedures (pva sounds for plausible values algorithsm)
  for (ii in 1:length(pctmiss_vec)){
    tmp_list[[ii]] = tmpj_list
    names(tmp_list)[ii] = paste0("pm",pm_vec[ii])
    for (jj in 1:length(pva_vec)){
      # Extract the observed data
      df_ii = list_get_obs$list_obsdf[[ii]]; obsdf = df_ii[,1:J]

      # Get the needed arguments for the imputation procedure
      what_tmp = methods_list$procedure[[jj]]
      args_tmp = methods_list$args[[jj]]

      if (what_tmp == "mice" | what_tmp == "parlmice"){
        args_tmp$data = obsdf
        # Note that data.init is commented out because it was leading to the same pvs for each imputation. Reason unknown
        #args_tmp$data.init = list_get_complete$dfcom[,1:J]
        if (args_tmp$method == "bygroup"){
          args_tmp$group = df_ii$subpop
        }
        if ("blocks" %in% names(args_tmp)){
          args_tmp$blocks = list(block1 = names(obsdf)[1:J])
        }
      } # end if (what_tmp == mice)

      if (what_tmp == "stratamelia"){
        args_tmp$x = obsdf
        args_tmp$strata = df_ii$subpop
      }

      if(what_tmp == "amelia"){
        args_tmp$x = obsdf
        if(!is.null(args_tmp$cs)){
          names_obs = names(list_get_obs$list_obsdf[[ii]])
          col_subpop = which(names_obs == args_tmp$cs)
          args_tmp$x = transform(args_tmp$x, subpop = as.factor(df_ii$subpop))
        }
        if(!is.null(args_tmp$empri)){
          args_tmp$empri = args_tmp$empri*nrow(obs_df)
        }
      }

      if(what_tmp == "amelia-FE"){
        args_tmp$x = dummy_columns(df_ii, select_columns = "subpop", remove_first_dummy = TRUE) %>%
          select(starts_with("Y"),starts_with("X"),starts_with("subpop_"))
        what_tmp = "amelia"
      }

      if (what_tmp=="EMs_LPA"){
        args_tmp$obsdf = list_get_obs$list_obsdf[[1]]
        args_tmp$z = z
        args_tmp$data_conditions = data_conditions
        args_tmp$rep = rep
      }

      # Run the imputation procedure
      obj_call = do.call(what = what_tmp, args = args_tmp)
      print(summary(obj_call))

      # Transform the data, as needed
      if (what_tmp == "amelia"){
        longimp_df = obsdf
        longimp_df = transform(longimp_df, subpop =df_ii$subpop, .imp = 0)
        if (obj_call$m == 1){
          tmp_df = data.frame(obj_call$imputations$imp1) %>% select(starts_with("Y"),starts_with("X"))
          tmp_df = transform(longimp_df, subpop = as.integer(df_ii$subpop),.imp = 1)
          longimp_df = rbind(longimp_df, tmp_df)
        } else {
          for (mm in 1:obj_call$m){
            tmp_df = obj_call$imputations[[mm]] %>% data.frame() %>% select(starts_with("Y"),starts_with("X"))
            tmp_df = transform(tmp_df, subpop = as.integer(df_ii$subpop), .imp = mm)
            longimp_df = rbind(longimp_df, tmp_df)
          } #end for mm
        } #end if (obj_call$m == 1)
        mids_obj = as.mids(longimp_df)
        obj_call = mids_obj
      }

      tmp_list[[ii]][[jj]] = obj_call


      if (save_it == TRUE){
        require(MplusAutomation)

        # Save a  ".dat" file of the imputed data set for Mplus
        list_mice = mice::complete(obj_call, "all")
        for (m in 1:length(list_mice)){
          list_mice[[m]] = transform(list_mice[[m]], subpop = list_get_complete$dfcom$subpop)
        }

        dat_folder =paste0("Imputed data/pm",pm_vec[ii],"/pva",pva_vec[jj])
        dat_fname = paste0("impdf rep", rep , " z", z, " pm", pm_vec[ii], " pva", pva_vec[jj])

        prepareMplusData(list_mice,
                         filename = paste0(temp_wd_rep_z,"/",dat_folder,"/",dat_fname), inpfile = FALSE,
                         overwrite = TRUE, imputed = TRUE)


        for (m in 1:length(list_mice)){

          out_list$dffolderfiles = rbind(out_list$dffolderfiles,
                                         data.frame(rep = rep, z = z,
                                                    folders = dat_folder,
                                                    files = paste0(dat_fname,"_imp_",m, ".dat"),
                                                    pm = pm_vec[ii], pva = pva_vec[jj], m = m) )

        } # end for (m in ...)
      } #end if(saveit==TRUE)

    } # END jj = 1,...,length(pva_vec)
  } #END ii = 1,...,length(pm_vec)

  out_list$obj_call=tmp_list

  return(out_list)
}
