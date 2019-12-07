#' Reads in .out files and saves the results
#'
#' @param tracker_df (data.frame) the tracker
#' @param estwd (character) working directory for the estimates
#' @param cl Optional cluster object .
#' @param args_readModels (list) Additional arguments to MplusAutomation::readModels()
#'
#' @return tracker_df
#' @export
#'
#' @examples
tracker2_save_readMplus<-function(tracker_df, estwd, cl = NULL, args_readModels = NULL){
  require("parallel")
  require("MplusAutomation")

  tracker_df = transform(tracker_df, estwd = as.character(estwd), estfolder = as.character(outfolder),
                         estfile = paste0("est-", gsub(" ","-",tolower(data_type)), "-",outfile,".RData"))

  x_vec = with(tracker_df, which(converged==TRUE))
  tracker_df$estwd[-x_vec] = as.character(NA)
  tracker_df$estfolder[-x_vec] = as.character(NA)
  tracker_df$estfile[-x_vec] = as.character(NA)


  print("Creating nessary est-files directories:")
  tmp<-pblapply(X = x_vec,
                FUN = function(x){
                  track_x = tracker_df[x,];
                  with(track_x,
                       dir.create(path = paste0(estwd,estfolder),
                                  recursive = T,
                                  showWarnings = F)
                  );
                }
    )






  print(" ")
  print(" ")
  print(paste0("Importing estimates from .out files and saving the resulting MplusAutomation::readModels() lists for ", length(x_vec), " models..."))

  if (is.null(cl)){
  tmp = pblapply(X = x_vec,
                           FUN = function(x){
                                  require(MplusAutomation);
                                  out_Mplus_x = with(tracker_df[x, ], MplusAutomation::readModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                                                                 filefilter = as.character(outfile),
                                                                                 what = ifelse(is.null(args_readModels$what),"all", args_readModels$what),
                                                                                 quiet = ifelse(is.null(args_readModels$quiet), TRUE, args_readModels$quiet)));

                                  list_estimates = list(tracker_row = tracker_df[x,], out_Mplus = out_Mplus_x)
                                  save(list_estimates, file = with(tracker_df[x,],paste0(estwd,"/",estfolder,"/",estfile)));
                               },
                           cl = NULL
                           )
  } else {
    print(cl)
    pb <- pbapply::timerProgressBar(max = length(x_vec), style = 1, width = getOption("width")/4)
    progress <- function(x){setTimerProgressBar(pb, x)}
    opts <- list(progress = progress)
    tmp<-foreach(x = x_vec, .packages= c("MplusAutomation"), .inorder = TRUE, .options.snow = opts) %dopar% {
      out_Mplus_x = with(tracker_df[x, ], MplusAutomation::readModels(target = paste0(as.character(outwd),"/",as.character(outfolder)),
                                                                      filefilter = as.character(outfile),
                                                                      what = ifelse(is.null(args_readModels$what),"all", args_readModels$what),
                                                                      quiet = ifelse(is.null(args_readModels$quiet), TRUE, args_readModels$quiet)));

      list_estimates = list(tracker_row = tracker_df[x,], out_Mplus = out_Mplus_x)
      save(list_estimates, file = with(tracker_df[x,],paste0(estwd,"/",estfolder,"/",estfile)));
      return(NULL)
    }
    close(pb)
  }
  rm(tmp)



  return(tracker_df)

}
