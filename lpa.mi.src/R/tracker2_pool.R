#' Pool analysis files given tracker information
#'
#' @param tracker_df (data.frame) the tracker
#' @param pool_wd - (character) Where the Rdata file with the pooled estimates and statistics should be saved
#' @param data_conditions (data.frame) the data conditoins
#' @param M_max - (integer) Maximum number of imputations to use
#' @param cl (cluster) A cluster variable from the SNOW package.
#' @param ... Optional argument to send to pool_subroutine (e.g., M_max, M_min, calculate.ARIU, and/or calculate.KL)
#'
#' @return tracker_df - updated with saved pooled information
#' @export
#'
#' @examples
tracker2_pool<-function(tracker_df, pool_wd, data_conditions,
                        cl = NULL, M_max = Inf, ...){


    #tracker2_pool
    # Inputs
    # M_max - Maximum number of imputations to use
    # M_min - Minimum number of imputation to report results. Defaults to 5
    # pool_wd -
    # data_conditions - data frame with the data conditoins
    # methods_list -
    # calcARIU - TRUE
    # calcKL - TRUE
    # cl -

    # Output
    # pooled_parameters_df - contains final parameter estimates as well as the deviance from populatoin values and whether properly covered
    # pooled_summary_df - contins means of LL, AIC, BIC, entropy, aBIC. If the the corresponing set to true, also contains:
    # (a) ARIV (Raghunathan & Rubin, 1991)
    # (b) ARIU (my alternatrive to ARIV which ARIV should approximate with proper imputations)
    # (c) KL (K-L divergence to population estimate) including the standard error of the estimate from the Monte Carlo integration


    require(plyr)
    require(tidyverse)
    require(pbapply)
    require(lpa.mi.src)
    require(data.table)
    require(snow)
    require(doSNOW)
    require(foreach)

    ##### Create the relevant directories ####
    # Generate pick1_df, which is simply a subset of the tracker for the imputed data sufficient giving one row for where all the est results are stored
    pick1_df = ddply(subset(tracker_df, data_type == "Imputed data" & converged == TRUE), .(rep,z,pm,pva,kfit), summarize,
               m = min(m),
               MM_max = sum(converged))
    pick1_df = merge(x = pick1_df,
               y = tracker_df[,c("rep","z","pm","pva","kfit","m","outfolder","estwd","estfolder")], by = c("rep","z","pm","pva","kfit","m"), all.x = TRUE, all.y = FALSE, sort = FALSE)
    # Add where hte pooling results will be stored and tidy up
    pick1_df = transform(pick1_df, poolwd = pool_wd,
                       poolfolder = outfolder,
                       poolfile = paste0("pooled-imputed-data-Free-Naive-LPA M", M_max, " rep", rep, " z",z, " pm",pm," pva", pva, " k",kfit,".RData"))
    pick1_df = pick1_df[, -match(c("outfolder","m"),names(pick1_df))]
    pick1_df = transform(pick1_df, M_sufficient = FALSE, data_type = "Imputed data", converged = TRUE)
    # Create the directories
    tmp = lapply(X = 1:nrow(pick1_df), FUN = function(x){with(pick1_df[x,], dir.create(paste0(poolwd,"/",poolfolder), recursive = TRUE))})
    rm(tmp)


    ##### Conduct Pooling ####
    print(paste0("Pooling for replication ", tracker_df$rep[1],":"))
    if (is.null(cl)){
      tmp = pblapply(X = 1:nrow(pick1_df),
                     FUN = function(x){
                       tracker_x = do.call(what = "pool_subroutine",
                                           args = list(oneline_df = pick1_df[x,],
                                                       data_conditions = data_conditions,
                                                       M_max = M_max,
                                                       ...));
                                       return(tracker_x);
                     },
                     cl = NULL)
    } else {
      print(cl)
      pb <- pbapply::timerProgressBar(max = nrow(pick1_df), style = 1, width = getOption("width")/4)
      progress <- function(x){setTimerProgressBar(pb, x)}
      opts <- list(progress = progress)
      tmp<-foreach(nn = 1:nrow(pick1_df),
                   .packages = c("lpa.mi.src","MplusAutomation","MASS","mixtools"),
                   .inorder = TRUE, .options.snow = opts) %dopar% {
                     tracker_x = do.call(what = "pool_subroutine",
                                         args = list(oneline_df = pick1_df[nn,],
                                                     data_conditions = data_conditions,
                                                     M_max = M_max, ...
                                                     ));
                     return(tracker_x);
      }
      close(pb)
    }

    new_pick1_df = data.table::rbindlist(tmp)

    tracker_df = merge(x = tracker_df, y = new_pick1_df, by = c("rep","z","pm","pva","kfit","data_type","converged","estwd","estfolder"), sort = FALSE, all.x = TRUE)
    return(tracker_df)
}
