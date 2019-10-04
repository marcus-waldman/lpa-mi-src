#' Imputation of MAR LPA Model using EM with Sampling
#'
#'
#' @param obsdf (data.frame) Observed data
#' @param M (integer) Number of imputations to construct
#' @param z (integer) Data condition number
#' @param data_conditions (data.frame) Simulation conditions pertaining to the complete data
#' @param starts0 (integer) Defaults to 20. Number of intial starts for Mplus when fitting FMM models.
#' @param bayes_boot (logical) Defaults to TRUE. Whether or not to perform Bayesian Bootstrap.
#' @param boot_2x (logical) Defaults to FALSE. Whether or not to implement a double Bootstrap.
#' @param tempdir_mplus (character) Defaults to that returned by getwd(). Temporary directory to store .dat, .inp, and .out files. Deleted upon successful termination.
#' @param rep (integer) Defaults to NULL. Replication number for temprorary directory naming purpose. If one not provided, then a random identifier is constructed.
#' @param itermax (integer) Defaults to 2. Number of iterations
#' @param multiplier (integer) Defaults to 64. Number to multiply starts by at each iteration.
#' @param cl Optional SOCKcluster object.
#'
#' @return (mids) An imputed data set of time (mids)
#' @export
#'
#' @examples
#' EM_with_sampling(obsdf, M, z, data_conditions)
EMs_LPA<-function(obsdf, M, z, data_conditions,
                               starts0 = 20, bayes_boot = TRUE, boot_2x = FALSE,
                               tempdir_mplus = getwd(), rep = NULL, itermax = 2, multiplier = 64, cl = NULL){



    require("dplyr")
    require("tidyverse")
    require("data.table")
    if (bayes_boot == TRUE){require("gtools")}
    require("MASS")
    require("mice")
    require("MplusAutomation")
    require("ids")
    require("lpa.mi.src")
    require("pbapply")
    require("foreach")
    require("matrixcalc")


    N_Processors = max(c(length(cl),1))

    # Add hash to temporary filename
    if(is.null(rep)){rep = ids::random_id(n = 1, bytes = 2)}
    mwid = ids::random_id(n = M, bytes = 3)

    # Get population parameters
    pop_params_z = lpa.mi.src::get_FMM_params(z = z, data_conditions = data_conditions)
    N_z = data_conditions$N[z]
    J_z = pop_params_z$J
    K_z = pop_params_z$K_z
    J_Y_z = pop_params_z$J_Y_z
    J_Xcom_z = pop_params_z$J_Xcom_z
    J_Xinc_z = pop_params_z$J_Xinc_z


    # Construct sampling weights based on bootstrapping and add to observed data
    if (bayes_boot == TRUE){
      Wmat_mn = N_z*gtools::rdirichlet(M, rep(1,N_z))
      if(boot_2x == TRUE){
        Wmat_mn = t(sapply(1:M, function(x){return( N_z*gtools::rdirichlet(1,Wmat_mn[x,]))}))
      }
      wdf = data.frame(t(Wmat_mn));
    } else {

       get_wgt_boot <- function(x){
          idx_boot = sample(1:N_z, N_z, replace = TRUE)
          if (boot_2x == TRUE){
            idx_boot = sample(idx_boot, N_z, replace = TRUE)
          }
          tab_boot = table(idx_boot)
          wgt_boot = rep(0,N_z)
          wgt_boot[as.integer(names(tab_boot))] = as.integer(tab_boot)
          return(wgt_boot)
       }
       wdf = data.frame(sapply(1:M, "get_wgt_boot"))
    }
    names(wdf) = paste0("W",mwid)
    wobsdf = cbind(obsdf[,1:J_z], wdf)

    # Create A new file and write out the data.
    path_tmpdir = paste0(tempdir_mplus,"/temp-EMs-FMM-rep", rep, "-z",z)
    dir.create(path = path_tmpdir, recursive = TRUE)
    obj_Mplus = prepareMplusData(df = wobsdf,
                                 filename = paste0(path_tmpdir,"/EMs_FMM rep",rep," z", z,".dat"))

    #### Fit bootstrapped models until convergence with the assistance of a tracker (trackdf) ####
      trackdf = data.frame(path = path_tmpdir,
                              stem = paste0("abayes rep",rep, " z", z),
                              m = mwid,
                              starts0 = starts0,
                              converged = rep(FALSE, M))

      S_compare_fn<-function(S_z,K_z){
        k_df = expand.grid(k1 = seq(1,K_z-1), k2 =  seq(2,K_z)) %>% subset(k1<k2)
        S1_z = lapply(1:nrow(k_df), FUN = function(x){S_z[,,k_df$k1[x]]})
        S2_z = lapply(1:nrow(k_df), FUN = function(x){S_z[,,k_df$k2[x]]})
        return(identical(S1_z, S2_z))
      }
      with_fn <- function(kk = NULL){
        prefix_txt = paste0(names(obsdf)[seq(J_Y_z+1, J_z)], " WITH ", names(obsdf)[1:J_Y_z],"*")
        if(is.null(kk)){
          out_txt = paste0(prefix_txt,";\n")
        } else {
          if( with(pop_params_z, S_compare_fn(S_z = S_z, K_z = K_z)) == TRUE ){kk = 1}
          out_txt = paste0(prefix_txt, " (s12_",kk,");\n")
        }
        return(out_txt)
      }
      var_fn<-function(kk = NULL){
        prefix_txt = paste0(names(obsdf)[1],"-",names(obsdf[J_z]),"*")
        if(is.null(kk)){
          out_txt = paste0(prefix_txt,";\n")
        } else {
          if( with(pop_params_z, S_compare_fn(S_z = S_z, K_z = K_z)) == TRUE ){kk = 1}
          out_txt = paste0(prefix_txt, " (s11_",kk,");\n")
        }
        return(out_txt)
      }
      iter = 0; flag = 0; Mvec = 1:M; processor.cores = 1; flag_cl = 0;
      while(flag == 0){

          iter = iter+1

          print("Writing .inp files")
          tmp<-pblapply(X = Mvec, FUN = function(i){ #Write out .inp files for replicates that have not converged
              # Title:
              title_txt = c("TITLE:\n",
                            with(trackdf,paste0(stem[i], " m",m[i],"\n")))
              # Data
              data_txt = c("DATA:\n",
                           paste0("FILE = ","'EMs_FMM rep",rep," z", z,".dat';\n") )
              #Variable
              variable_txt = c("VARIABLE:",
                               paste0("NAMES = \n", paste(" ",names(wobsdf), "\n", collapse = " "), ";\n"),
                               paste0("USEV = ", names(wobsdf)[1],"-",names(wobsdf)[J_z],";\n"),
                               paste0("CLASSES = c(",K_z,");\n"),
                               paste0("WEIGHT = W",trackdf$m[i],";\n"),
                               "MISSING = .;\n")
              # Analysis
              starts0_i = trackdf$starts0[i]
              analysis_txt = c("ANALYSIS:\n",
                               "TYPE = mixture;\n",
                               "ESTIMATOR = mlr;\n",
                               paste0("STARTS = ",starts0_i," ", floor(8*starts0_i/20),";\n"),
                               paste0("PROCESSORS = ", processor.cores,";"))
              # Model
              with_txt = paste0(names(obsdf)[seq(J_Y_z+1, J_z)], " WITH ", names(obsdf)[1:J_Y_z],"*;\n")
              var_txt = paste0(names(obsdf)[1],"-",names(obsdf[J_z]),"*;\n")

              overall_txt = c("MODEL:\n",
                              "%OVERALL%",
                              with_fn(),
                              var_fn())
              model_txt = overall_txt
              c_txt = as.vector(sapply(1:K_z, FUN = function(x){c(paste0("%c#",x,"%"), with_fn(x), var_fn(x))}))
              model_txt = c(overall_txt,
                            c_txt)
              #Output
              output_txt = c("OUTPUT:\n")
              #Savedata
              savedata_txt = c("SAVEDATA:\n")
              # Combine everything and write out
              Mplus_txt = c(title_txt,
                            data_txt,
                            variable_txt,
                            analysis_txt,
                            model_txt,
                            output_txt,
                            savedata_txt)
              fileConn<-file(with(trackdf, paste0(path[i],"/",stem[i]," m",m[i],".inp")))
              writeLines(Mplus_txt, fileConn)
              close(fileConn)

              return(NULL)
          })
          rm(tmp)

          print("Applying EM algorithm to bootstrapped data...")
          if (processor.cores>1 | is.null(cl)){
            tmp<-pblapply(X = Mvec, FUN = function(i){ # Fit the models that have not reached convergence
                    runModels(target = with(trackdf, paste0(path[i],"/",stem[i], " m",m[i],".inp")),
                              showOutput = FALSE,
                              logFile = NULL)
                    return(NULL)
                    }
                )
            rm(tmp)
          } else {
            pb <- pbapply::timerProgressBar(max = length(Mvec), style = 1, width = getOption("width")/4)
            progress <- function(x){setTimerProgressBar(pb, x)}
            opts <- list(progress = progress)
            tmp <- foreach(i = Mvec,  .packages = c("MplusAutomation"),
                           .inorder = TRUE, .options.snow = opts) %dopar% {

                             runModels(target = with(trackdf, paste0(path[i],"/",stem[i], " m",m[i],".inp")),
                                       showOutput = FALSE,
                                       logFile = NULL)

                             return(NULL)

                       } #end foreach
            rm(tmp)
          } # end if/else processor.cores

          # Assess convergence of the fit models
          print("Assessing convergence")
          list_checkc = pblapply(X = Mvec,
                               FUN = function(x){check_convergence(file = with(trackdf,paste0(stem[x]," m",m[x],".out")),
                                                                   folder_wd = path_tmpdir,
                                                                    starts_txt = with(trackdf, paste0(starts0[x]," ",floor(8*starts0[x]/20), ";")),
                                                                   Rcond_min = ifelse(iter<=1, 1E-6, 1E-8))})
          checkc_df = rbindlist(list_checkc) %>%
                      transform(converged = (normal==TRUE &  poor_Rcond==FALSE))

          trackdf$converged[Mvec] = checkc_df$converged

          Mvec = which(trackdf$converged == FALSE)

          # Determine whether to continue within loop or exit
          if (length(Mvec)==0 ){
            flag = 1
          } else {
            if (iter==itermax){flag = 1}
            trackdf$starts0[Mvec] = multiplier*trackdf$starts0[Mvec]
            print(paste0(length(Mvec)," model failed to converge."))
          }


      }# end while


    # Read models in temporary directory
      print(paste0("In total, ", sum(trackdf$converged)," converged. Reading in the coverged results..."))
      trackdf = subset(trackdf, converged == TRUE)
      if(nrow(trackdf)==0){return(out_mids=NULL)}

      mplus_list = pblapply(X = 1:nrow(trackdf),
                           FUN = function(x){
                              readModels(target = path_tmpdir,
                                         filefilter = with(trackdf[x, ],paste0(stem," m",m)),
                                         what = c("parameters"), quiet = TRUE)
                              })
      names(mplus_list) = gsub(" ",".",with(trackdf,paste0(stem," m",m,".out")))


      # Check that the within class covariance matrices are all positive definite
      problem_vec = sapply(X = 1:nrow(trackdf),
             FUN = function(x){
              params_df = mplus_list[[x]]$parameters$unstandardized
              Qlist_m = Mplus2Qlist(params_df)
              return(sum(apply(Qlist_m$S,3,"rcond")<=1E-6))
             }
          )
      trackdf = transform(trackdf, problem_S = problem_vec>0)
      print(with(trackdf, paste0("Of the converged results, ", sum(problem_S), " displayed poor within-class condition numbers and will not be used to generate plausible values")))
      mplus_list = mplus_list[which(problem_vec==0)]
      trackdf = subset(trackdf, problem_S == FALSE)
      if (nrow(trackdf)==0){return(out_mids=NULL)}



    list_pvs = NULL
    # Get missing data pattern information
      out_mdpatt = get_mdpatterns(df = obsdf[,1:J_z])

    # Construct plausible values
      doPV<-function(x){
          if (x == 0){
            Yimp_m = obsdf[,1:J_z]
          } else {
            params_df = mplus_list[[x]]$parameters$unstandardized
            Qlist_m = Mplus2Qlist(params_df)

            cprobs_m = cprobs(Y_i = obsdf[,1:J_z],
                              pi_vec = Qlist_m$pi,
                              mu_mat = Qlist_m$mu,
                              S_array = Qlist_m$S,
                              list_mdpattern = out_mdpatt)
            if (sum(is.na(cprobs_m))==0){
              cdraws_vec = get_pseudoc_draw(CPROBS_df = cprobs_m)
            } else {
              stop("Could not calculate psuedo-class draws. NA values observed")
            }


            Yimp_m = obsdf[,1:J_z]
            for(k in 1:K_z){
              ivec_k = which(cdraws_vec==k)
              idRvec_k = with(out_mdpatt, unique(patterns_i$idR[ivec_k])); idRvec_k = idRvec_k[idRvec_k>0];
              for (r in idRvec_k){
                ivec_rk = with(out_mdpatt, intersect(ivec_k, which(patterns_i$idR == r)))
                Y_r = obsdf[ivec_rk,1:J_z]
                N_r = nrow(Y_r)
                idU = with(out_mdpatt, U_r[[r]]); idV = with(out_mdpatt, V_r[[r]])
                list_UV_mr = with(Qlist_m, get_conditional(mu = mu[,k], S = S[,,k],
                                                           U = idU,
                                                           V = idV)
                )
                if(!is.positive.definite(list_UV_mr$S_eps)){
                  stop("S_eps not positive definite")
                }
                X_mr = as.matrix(cbind(rep(1,N_r),Y_r[,idV]))
                E_mr = with(list_UV_mr, X_mr%*%beta)
                eps_mr = with(list_UV_mr, mvrnorm(n = N_r, mu = 0.*idU, S = S_eps))
                Yimp_m[ivec_rk,idU] = E_mr+eps_mr
              } #end for i
            } # end for k
          }
          return(Yimp_m)
      }


  print("Sampling plausible values...")
  list_pvs = pblapply(X = 1:nrow(trackdf), FUN = function(x){try(doPV(x),T)})
  ids_keep = sapply(X = 1:length(list_pvs),FUN = function(x){is.data.frame(list_pvs[[x]])}) %>% which()
  list_pvs = list_pvs[ids_keep]

  print(paste0("Plausible values could not be calculated for ",
               nrow(trackdf)-length(list_pvs), " of the ", nrow(trackdf), " elgible results."))


  # Clean up and spit out a 'mids' object
  print("Cleaning up...")
  system(paste0("rm -r '",gsub("/","\\\\",path_tmpdir),"' "))

  impdf = NULL
  if(!is.null(list_pvs)){
    for (x in 1:length(list_pvs)){
      list_pvs[[x]] = list_pvs[[x]]%>% transform(subpop = obsdf$subpop, .imp = x)
    }
    impdf = rbindlist(list_pvs)
  }

  obsdf2 = obsdf %>% transform(.imp = 0)
  keep_cols = match(names(impdf),names(obsdf2)); keep_cols = keep_cols[!is.na(keep_cols)]
  impdf = rbind(impdf, obsdf2[,keep_cols])

  out_mids = as.mids(impdf)
  return(out_mids)


} #end function
