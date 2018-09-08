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
#' @param p (integer) Processor number (defaults to 1).
#' @param save_it (logical) if TRUE, then it saves the data set and the following also must be specificed
#' @param temp_wd_p_vec (character vector). Required if save_it==TRUE.  processor-specific temporary directory
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
                           rep = NA, p = NA, save_it = FALSE, temp_wd_p_vec = NULL,...){

    require(Amelia)
    require(mice)
    require(miceadds)
    require(mitools)


  # Inputs: list_get_obs, methods_list, z, data_conditions, pctmiss_vec, save_it, rep, p, temp_wd_p


# list_get_complete = out_comp; list_get_obs = out_obs;
# save_it = TRUE; p = 1; ; rep = 1;

  temp_wd_p = temp_wd_p_vec[p];

  if (save_it==TRUE){require(MplusAutomation)}


    if(save_it == TRUE & (is.na(rep) | is.null(temp_wd_p) | is.na(p))){
      stop("rep, p, & temp_wd_p arguments cannot be NULL with save_it = TRUE")
    }


    # Pre-processing, get needed variables for later
    J_Y = data_conditions$J_Y[z]; J_Xcom = data_conditions$J_Xcom[z]; J_Xinc = data_conditions$J_Xinc[z]
    J = J_Y + J_Xcom + J_Xinc

    # Allocate space for the out_list
    out_list = list(obj_call=NULL,implist = NULL, dffolderfiles = NULL)
    tmp_list = list(NULL)
    vv = 0
    # Loop over A) pm =  1:length(pctmiss_vec), 2) pva = 1:Nprocedures (pva sounds for plausible values algorithsm)
    for (pm in 1:length(pctmiss_vec)){
     for (pva in 1:length(methods_list$procedure)){
          obsdf = list_get_obs$list_obsdf[[pm]]; obsdf = obsdf[,1:J]


          # Get the needed arguments for the imputation procedure
          what_tmp = methods_list$procedure[[pva]]
          args_tmp = methods_list$args[[pva]]
          if (what_tmp == "mice"){
            args_tmp$data = obsdf
    # Note that data.init is commented out because it was leading to the same pvs for each imputation. Reason unknown
            #args_tmp$data.init = list_get_complete$dfcom[,1:J]
            if (args_tmp$method == "bygroup"){
              args_tmp$group = list_get_complete$dfcom$class
            }
            if ("blocks" %in% names(args_tmp)){
              args_tmp$blocks = list(block1 = names(obsdf)[1:J])
            }
          }
          if (what_tmp == "stratamelia"){
            args_tmp$x = obsdf
            args_tmp$strata = list_get_complete$dfcom$subpop
          }
          if(what_tmp == "amelia"){
            args_tmp$x = obsdf
          }
          # Run the imputation procedure
          obj_call = do.call(what = what_tmp, args = args_tmp)

          # Transform the data, as needed
          if (what_tmp == "amelia"){
            longimp_df = obsdf
            longimp_df = transform(longimp_df, .imp = 0)
            if (obj_call$m == 1){
              tmp_df = data.frame(obj_call$imputations$imp1)
              tmp_df = transform(longimp_df, .imp = 1)
              longimp_df = rbind(longimp_df, tmp_df)
            } else {
              for (mm in 1:obj_call$m){
                tmp_df = data.frame(obj_call$imputations[,,mm])
                names(tmp_df) = names(x)
                tmp_df = transform(tmp_df, subpop = strata, .imp = mm)
                longimp_df = rbind(longimp_df, tmp_df)
              } #end for mm
            } #end if (obj_call$m == 1)
            mids_obj = as.mids(longimp_df)
            obj_call = mids_obj
          }

          vv = vv+1
          tmp_list[[vv]] = obj_call
          names(tmp_list)[vv] = paste0("mids_pm",pm,"pva",pva)

          # Add in subpop
          #obj_call = cbind.mids(obj_call, subpop = list_get_complete$dfcom$subpop)
          #print("hello world")
         #print(class(obj_call))


          if (save_it == TRUE){
            require(MplusAutomation)

            out_list$dffolderfiles = rbind(out_list$dffolderfiles,
                                           data.frame(folders = paste0("Imputed data/pm",pm,"/pva",pva),
                                                files = paste("impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat",sep = ""),
                                                data_condition = z, m = NA) )
            # # Save a copy of the obj_call
            # save(obj_call,
            #      file = paste0(temp_wd_p,"/Imputed data/pm",pm,"/pva",pva,"/mids p", p," z",z," rep",rep, " pm", pm, " pva", pva,".RData")
            # )


            # Save a  ".dat" file of the imputed data set for Mplus

            colMax = ifelse(what_tmp=="stratamelia",J+1,J)
            list_mice = mice::complete(obj_call, "all")
            for (m in 1:length(list_mice)){
              list_mice[[m]] = transform(list_mice[[m]], subpop = list_get_complete$dfcom$subpop)
            }
              #nms = names(list_get_obs$obsdf_z)
              #jkeep = which(startsWith(nms,"Y") | startsWith(nms,"X") | nms=="subpop")
            prepareMplusData(list_mice, #keepCols = jkeep,
                             filename = paste0(temp_wd_p,"/Imputed data/pm",pm,"/pva",pva,"/impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat"), inpfile = FALSE,
                             overwrite = TRUE, imputed = TRUE)


            dat_ms = read.delim(paste0(temp_wd_p,"/Imputed data/pm",pm,"/pva",pva,"/impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat"), header = FALSE)
            for (m in 1:nrow(dat_ms)){
              out_list$dffolderfiles = rbind(out_list$dffolderfiles,
                                             data.frame(folders = paste0("Imputed data/pm",pm,"/pva",pva),
                                                        files = paste("impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,"_imp_",m,".dat",sep = ""),
                                                        data_condition = z, m = m) )
            }
          }

      } # END pva = 1,...
    } #END pm = 1,...

    out_list$obj_call=tmp_list

    return(out_list)
}
