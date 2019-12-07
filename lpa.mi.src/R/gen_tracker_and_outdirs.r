#' Generate Directories for Mplus .out Files and File Tracker
#'
#' @param outwd (character) where to store the out files
#' @param datwd (character) where the dat files are located
#' @param rdatawd (character) where the Rdata file with the list resides
#' @param data_type  (character) any of 'Complete data', 'Observed data', 'Stacked data' or 'Imputed data'
#' @param rep_vec (numeric) replications to track
#' @param z_vec (numeric) data conditions to track
#' @param data_conditions  (data.frame) must include the variable z
#' @param fix_vec (logical) Either FALSE (free estimation) or 1 (fixed esitmation) or c(0,1) (i.e, both).
#' @param pm_vec (numeric) Defaults to NA. Percent missing identifier.
#' @param pva_vec (numeric) Defaults to NA. Imputation algorithm identfier
#' @param M (integer) Defaults to NA. Number of imputations.
#' @param dk_vec  (numeric). Defaults to 0.
#' @param starts0 (integer). Defaults to 32.
#' @param mkdir_out (logical). Defaults to TRUE. Whether or not to make the folder for where the outfiles will be stored.
#'
#' @return tracker_df (data.frame)
#' @export
#'
#' @examples get_tracker_and_outdir(outwd, datwd, data_type, rep_vec, z_vec, data_conditions)
gen_tracker_and_outdirs<-function(outwd, datwd, rdatawd, data_type, rep_vec, z_vec, data_conditions,
                         fix_vec=FALSE, pm_vec = NA, pva_vec = NA, M = NA, dk_vec = 0, starts0 = 32,
                         mkdir_out = TRUE){

    require(tidyverse)
    require(pbapply)

    # Make initial checks for any possible errors in inputs
    if (!data_type%in%c("Imputed data", "Observed data", "Complete data", "Stacked data")){stop("data_type must be 'Imputed data', 'Observed data', 'Complete data', or 'Stacked data'")}
    if (length(rep_vec)==1){warning("Provided rep_vec is one number. Check this is what is intended")}
    if (!"z"%in%names(data_conditions)){stop("z needs to be specified as a variable in data_conditions")}
    if (sum(fix_vec %in% c(0,1))!=length(fix_vec)){stop("fix_vec must only be 0 or 1")}

    print(" ")
    print(" ")
    print(paste0("Generating directories and tracker for the ", data_type, "."))

    # Ensure inputs are provided consistent with data_type
    if (data_type != "Complete data"){
      if (length(pm_vec)==0){stop(paste0("pm_vec must be specified with data_type = ",data_type ))}
      if (data_type != "Observed data" ){
        if (length(pva_vec)==0){stop(paste0("pva_vec must be specified with data_type = ", data_type))}
        if(data_type == "Imputed data" & is.na(M)){stop(paste0("M must be specified with data_type = ", data_type))}
      } else { #observed data
        pva_vec = NA #set to NA even if specified
      }
    } else { #Complete data conditions
      pm_vec <- pva_vec <- NA #set to NA even if specified
    }


#    # Initiate the tracker
   m_vec <-NA;
   if (data_type == "Imputed data"){m_vec = seq(1,M)}


          # Add the number of classes fit to the dat files to the tracker

    tracker_df = expand.grid(fix = fix_vec, rep = rep_vec, z = z_vec, pm = pm_vec, pva = pva_vec, dk = dk_vec, m =m_vec)
    tracker_df = transform(tracker_df, kfit = data_conditions$K[tracker_df$z] + tracker_df$dk)
    if(sum(tracker_df$kfit<0)>1){
      warning("dk specified so that negative classes fit. Removing these from tracker")
      tracker_df = subset(tracker_df, kfit>=0)
    }

    # Rearrange variables, sort, and add starts, and other diagnostics
    tracker_df = tracker_df %>%  arrange(fix,rep, z, pm, pva, dk, m) %>%  transform(nn = 1:nrow(.))
    tracker_df = tracker_df[,c("nn","fix","rep","z","pm","pva","m","kfit","dk")] %>% transform(starts0 = starts0)
    tracker_df = tracker_df %>% transform(data_type = data_type, datwd = datwd, outwd = outwd, rdatawd = rdatawd)

    # Create the datfolder
    tracker_df = transform(tracker_df, datfolder = paste0("/rep", rep,"/z",z,"/",data_type))
    if (data_type != "Complete data"){ #observed, stacked, or imputed
      tracker_df = transform(tracker_df, datfolder = paste0(datfolder,"/pm",pm))
      if (data_type != "Observed data"){ # stacked or imputed
        tracker_df = transform(tracker_df, datfolder = paste0(datfolder,"/pva",pva))
      }
    }

    # Create the outfolder
    tracker_df = transform(tracker_df, outfolder = paste0(datfolder,"/k",kfit))

    # Create rdata folder
    tracker_df = transform(tracker_df, rdatafolder = "rdata-files.zip")

    # Specify the dat file name to load into Mplus
    if(data_type == "Complete data"){
      tracker_df = transform(tracker_df, datfile = paste0("dfcom rep",rep," z",z,".dat"))
    }
    if (data_type == "Observed data"){
      tracker_df = transform(tracker_df, datfile = paste0("obsdf rep",rep," z",z," pm", pm, ".dat"))
    }
    if (data_type == "Stacked data"){
      tracker_df = transform(tracker_df, datfile = paste0("stackdf rep",rep," z",z," pm", pm, " pva", pva,".dat"))
    }
    if(data_type == "Imputed data"){
      tracker_df = transform(tracker_df, datfile = paste0("impdf rep",rep," z",z," pm", pm, " pva", pva, "_imp_",m,".dat"))
    }


    # Specify the out file name to write out the Mplus fitting results
    if(data_type == "Complete data"){
      tracker_df = transform(tracker_df, outfile = paste0(ifelse(fix==0,"Free","Fixed"),"-Naive-LPA rep",rep," z",z, " k", kfit))
    }
    if (data_type == "Observed data"){
      tracker_df = transform(tracker_df, outfile = paste0(ifelse(fix==0,"Free","Fixed"),"-Naive-LPA rep",rep," z",z," pm", pm, " k", kfit))
    }
    if (data_type == "Stacked data"){
      tracker_df = transform(tracker_df, outfile = paste0(ifelse(fix==0,"Free","Fixed"),"-Naive-LPA rep",rep," z",z," pm", pm, " pva", pva, " k", kfit))
    }
    if(data_type == "Imputed data"){
      tracker_df = transform(tracker_df, outfile = paste0(ifelse(fix==0,"Free","Fixed"),"-Naive-LPA rep",rep," z",z," pm", pm, " pva", pva, " m",m, " k", kfit))
    }

   # Specify the
    if(data_type == "Complete data"){
      tracker_df = transform(tracker_df, rdatafile = paste0("list-complete rep", rep, " z", z, ".RData"))
    }
    if (data_type == "Observed data"){
      tracker_df = transform(tracker_df, rdatafile = paste0("list-observed rep", rep," z",z," pm", pm, ".RData"))
    }
    if (data_type == "Stacked data"){
      tracker_df = transform(tracker_df, rdatafile = NA)
    }
    if(data_type == "Imputed data"){
      tracker_df = transform(tracker_df, rdatafile = paste0("list-imputed rep", rep, " z",z," pm", pm, " pva", pva,".RData"))
    }


    # Create the directories
    if(mkdir_out==TRUE){
      dir_df = subset(tracker_df, (m == 1 | is.na(m)))
      tmp = with(dir_df,
                 pblapply(X = 1:nrow(dir_df),
                          FUN = function(x){dir.create(paste0(as.character(outwd[x]),"/",as.character(outfolder[x])), recursive = TRUE)}
                 )
                )
    }

return(tracker_df)
}
