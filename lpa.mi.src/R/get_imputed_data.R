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
#'       (A) obj_call - (mids) mids object for imputed data for the
#'           (i)   z-th data condition number, 
#'           (ii)  pm-th percent missing
#'           (iii) pva-th imputation method (see methods_list for details)
#'       (B) dffolderfiles (data.frame) with the files and folders for each of the save imputed data sets
#' @export
#' @examples
#' get_imputed_data(z,list_get_obs, list_get_complete, methods_list, data_conditions, save_it = FALSE)

get_imputed_data<-function(z, list_get_obs, list_get_complete, methods_list, data_conditions, pctmiss_vec, 
                           rep = NA, p = NA, save_it = FALSE, temp_wd_p_vec = NULL,...){

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

        
# pm = 1
# pva = 1


    # Loop over A) pm =  1:length(pctmiss_vec), 2) pva = 1:Nprocedures (pva sounds for plausible values algorithsm)
    for (pm in 1:length(pctmiss_vec)){
     for (pva in 1:length(methods_list$procedure)){
          obsdf = list_get_obs$list_obsdf[[pm]]; obsdf = obsdf[,1:J]
          
          
          # Get the needed arguments for the imputation procedure
          what_tmp = methods_list$procedure[[pva]]
          args_tmp = methods_list$args[[pva]]
          if (what_tmp == "mice"){
            args_tmp$data = obsdf
            args_tmp$data.init = list_get_complete$dfcom[,1:J]
            if (args_tmp$method == "bygroup"){
              args_tmp$group = list_get_complete$dfcom$class
            }
            if ("blocks" %in% names(args_tmp)){
              args_tmp$blocks = list(block1 = names(obsdf)[1:J])
            }
          }
          if (what_tmp == "stratamelia"){
            args_tmp$x = obsdf
            args_tmp$strata = list_get_complete$dfcom$class
          }
          # Run the imputation procedure
          obj_call = do.call(what = what_tmp, args = args_tmp)
          out_list$obj_call[[pm]][[pva]] = obj_call
          
          
          if (save_it == TRUE){
            require(MplusAutomation)
            
            out_list$dffolderfiles = rbind(out_list$dffolderfiles, 
                                           data.frame(folders = "Imputed data", 
                                                files = paste("impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat",sep = ""), 
                                                data_condition = z, m = NA) )
            # Save a copy of the obj_call
            save(obj_call,
                 file = paste(temp_wd_p,"/Imputed data/mids p", p," z",z," rep",rep, " pm", pm, " pva", pva,".RData",sep = ""))
            
            
            # Save a  ".dat" file of the imputed data set for Mplus
            
            colMax = ifelse(what_tmp=="stratamelia",J+1,J)
            list_mice = mice::complete(obj_call, "all")
            invisible(
              prepareMplusData(list_mice, keepCols = 1:colMax,
                             filename = paste(temp_wd_p,"/Imputed data/impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat",sep = ""), inpfile = FALSE, 
                             overwrite = TRUE, imputed = TRUE)
              )

            dat_ms = read.delim(paste(temp_wd_p,"/Imputed data/impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,".dat",sep = ""), header = FALSE)
            for (m in 1:nrow(dat_ms)){
              out_list$dffolderfiles = rbind(out_list$dffolderfiles, 
                                             data.frame(folders = "Imputed data", 
                                                        files = paste("impdf p", p," z",z," rep",rep, " pm", pm, " pva", pva,"_imp_",m,".dat",sep = ""), 
                                                        data_condition = z, m = m) )
            }
          }
    
      } # END pva = 1,...
    } #END pm = 1,...


    return(out_list)
}