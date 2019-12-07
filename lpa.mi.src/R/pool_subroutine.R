#' Pooling subroutine called from tracker2_pool
#'
#' @param oneline_df - (data.frame) containing a single line in from pick1_df in tracker2_pool for only one of the imputed files
#' @param data_conditions - (data.frame)
#' @param M_max - (integer). The number of imputations to use
#' @param M_min - (integer). Optional. Defaults to a requirement of at least 5 imputations to report results.
#' @param calculate.ARIU (logical) Defaults to TRUE. Whether to report the average relative increase in variance with reference to the complete data set
#' @param calculate.KL (logical) Defaults to TRUE. Whether to calculate an estimate of the KL divergence between the model implied by the the pooled point estimates and the population model.
#' @param points.montecarlo (integer) Defaults to 5E4
#' @param deviance_coverage (logical) Defaults to TRUE. Will calculate deviance and coverage of each parameter from population values. Only available if kfit == 3
#' @param save_details (logical) Defaults to FALSE. If true, then also save the individual point estimates, within imputation asymptotic cov matrix, and between imputation estimate
#'
#' @return (data.frame) with updated line in tracker, including location of saved file analyzing results
#' @export
#'
#' @examples
pool_subroutine<-function(oneline_df, data_conditions, M_max,
                          M_min = 5, calculate.ARIU = TRUE, calculate.KL = TRUE, points.montecarlo=5E4,
                          deviance_coverage = TRUE, save_details = FALSE){

  require(lpa.mi.src)
  require(tidyverse)

  if(nrow(oneline_df)>1){stop("oneline_df can only be one row")}

  if(deviance_coverage==TRUE & oneline_df$kfit!=3){stop("deviance_coverage can only be set to TRUE with kfit=3.")}


  # Get the population values
  z_x = oneline_df$z
  if (deviance_coverage==TRUE){
    popvalues_df = lpa.mi.src::Plist2Mplus(z = z_x, data_conditions = data_conditions)
  }

  # Obtain a list of readModel resultsof length <= M_max
  rep_x = oneline_df$rep; pm_x = oneline_df$pm; pva_x = oneline_df$pva; kk_x = oneline_df$kfit
  est_rdata_vec = list.files(path = paste0(oneline_df$estwd,"/",oneline_df$estfolder), full.names = TRUE)
  if (length(est_rdata_vec)>M_max){est_rdata_vec=est_rdata_vec[1:M_max]}
  list_est_rdata = lapply(X = est_rdata_vec, FUN = function(x){load(x);return(list_estimates)})
  names(list_est_rdata) = paste0("mm",1:length(list_est_rdata))

  # Read in the point-value parameter estimates
  list_params = lapply(X = 1:length(list_est_rdata),
                       FUN = function(x){
                        if (length(list_est_rdata[[x]]$out_Mplus$errors)==0 ){
                           tmp_x = list_est_rdata[[x]]
                           tech1_x = get_tech1(out_readModels =  tmp_x$out_Mplus)
                           tmp_df = list_est_rdata[[x]]$out_Mplus$parameters$unstandardized
                           if (oneline_df$kfit==1){tmp_df = transform(tmp_df, LatentClass = 1)}
                           params_x =  tmp_df %>%
                             merge(tech1_x, by = c("paramHeader", "param", "LatentClass"), sort = FALSE) %>%
                             subset(!is.na(tech1)) %>% arrange(tech1)
                           return(params_x)
                         }
                       })
  xvec_okay = setdiff(1:length(list_est_rdata), which(sapply(list_params, is.null)) )
  if(length(xvec_okay)<length(list_est_rdata)){
    list_params = list_params[-which(sapply(list_params, is.null))]
  }
  names(list_params) = paste0("mm",1:length(list_params))

  # Witt_dfn Imputation Point Estimates
  # Obtain a Q-by-MM matrix giving the point estimates, where
  #       Q is the number of parameters in the model
  #       MM is the number of imputations (upto the maximum M_max)
  theta_mat = sapply(X = 1:length(list_params),
                     FUN = function(x){return(list_params[[x]]$est)})
  Q = nrow(theta_mat)
  MM = ncol(theta_mat)

  if(MM>=M_min){

            oneline_df$M_sufficient = TRUE

            # Pooled Point Estimates
            t_bar = apply(theta_mat, 1, "mean", na.rm = TRUE)

            # Witt_dfn Imputation Standard Errors
            # Obtain the within imputation asymptotic variance covariance matrices.
            # Save as a Q-by-Q-by-MM array
            V_array = sapply(X = xvec_okay,
                             FUN = function(x){
                               if (length(list_est_rdata[[x]]$out_Mplus$errors)==0 ){
                                 chi = list_est_rdata[[x]]$out_Mplus$tech3$paramCov; chi[is.na(chi)] = 0;
                                 tchi = t(chi); diag(tchi) = 0;
                                 V_m = chi+tchi;
                                 return(V_m);
                               }
                             },
                             simplify = "array"
            )

            # Pooled Within-Imputation Asymptotic Covariance Matrix
            Vbar = apply(V_array, c(1,2),"mean", na.rm = TRUE)

            # Between Imputation Deviance
            # Saved as a Q-by-Q-by-MM array
            B_array = sapply(X = 1:length(xvec_okay),
                             FUN = function(x){
                               dt_m = as.numeric(theta_mat[,x])-t_bar;
                               sumB = dt_m%*%t(dt_m); #See Schafer, 1997, p. 113
                               return(sumB)
                             },
                             simplify = "array"
            )

            # Estimate of Between-Imputation Variance (see Schafer, 1997, p. 113, for formula)
            sumB = apply(B_array, c(1,2), "sum", na.rm = TRUE)
            B = sumB/(MM-1)

            # Total variance estimate (approximation of posterior variance)
            Vtot = Vbar + B + B/MM

            # Obtain parameter estimates in Mplus form
            pooled_parameters_df = list_params$mm1[,c("paramHeader","param","LatentClass")] %>%
              transform(est = t_bar, se = sqrt(diag(Vtot)))
            if (deviance_coverage == TRUE){
                pooled_parameters_df = merge(pooled_parameters_df, popvalues_df, by = c("paramHeader","param","LatentClass"),
                                             all.x = TRUE, all.y = FALSE, sort = FALSE) %>%
                  transform(deviance = est - value,
                            covered = ifelse(value>=(est-1.959964*se) & value<(est+1.959964*se),TRUE,FALSE))
            } else {
                pooled_parameters_df = transform(pooled_parameters_df,
                                                 deviance = NA,
                                                 covered = NA)
            }

            # Estimate of the average relative increase in variance proved by Raghunathan & Rubin, 1991
            trBV = sum(diag(B%*%solve(Vbar))) #Fromula see Enders Equation 8.21
            ARIV = (1+(1/MM))*trBV/Q

            # Obtain the average of the ratio of between imputation variation to the complete data variation
            ARIU = NA
            if(calculate.ARIU == TRUE){
              complete_rdata_vec = list.files(paste0(oneline_df$estwd,"/rep",rep_x,"/z",z_x,"/Complete data/k",kk_x), full.names = TRUE)
              if (length(complete_rdata_vec)==1){
                load(complete_rdata_vec);
                chi = list_estimates$out_Mplus$tech3$paramCov; chi[is.na(chi)] = 0;
                tchi = t(chi); diag(tchi) = 0;
                U = chi+tchi;
                ARIU = sum(diag(B%*%solve(U)))/Q
              } # end if length(complete_rdata_vec) == 1
            }

            # Obtain the KL divergence as necessary
            KL_df = data.frame(est = NA, se = NA, points_montecarlo = NA)
            if(calculate.KL == TRUE & deviance_coverage == TRUE){
              Plist_x = lpa.mi.src::get_Plist(z = z_x, data = data_conditions)
              Qlist_x = lpa.mi.src::Mplus2Qlist(params_df = pooled_parameters_df %>% dplyr::select("paramHeader","param","LatentClass","est"))
              KL_df = lpa.mi.src::KLmixmvrnorm(P = Plist_x, Q = Qlist_x, n = points.montecarlo)
            }

            pooled_summaries_df = data.frame(Q = Q, MM = MM, ARIV = ARIV, ARIU = ARIU,
                                             KL = KL_df$est, se_KL = KL_df$se, points_montecarlo = points.montecarlo)

            details_list = NULL
            if(save_details == TRUE){
              details_list = list(theta_mat = theta_mat,
                                  V_array = V_array,
                                  B_array = B_array,
                                  Vbar = Vbar,
                                  B = B,
                                  list_params = list_params)
            }

            # Save the list out
            pooled_list = list(parameters = pooled_parameters_df,
                               additional = pooled_summaries_df,
                               details = details_list)


            save(pooled_list,
                 file = paste(oneline_df$poolwd,oneline_df$poolfolder,oneline_df$poolfile, sep = "/"))



  } else { #Define what to do if MM<M_min

    oneline_df$M_sufficient = FALSE
    oneline_df$poolfile = as.factor(NA)


  }  #end if(MM>=M_min)

  return(oneline_df)
}# end pool_routine
