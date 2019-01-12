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
#' @param time_max (integer) Defaults to 30 minutes maximum time of fitting all M models.
#'
#' @return (mids) An imputed data set of time (mids)
#' @export
#'
#' @examples
#' EM_with_sampling(obsdf, M, z, data_conditions)
EMs_LPA<-function(obsdf, M, z, data_conditions,
                               starts0 = 20, bayes_boot = TRUE, boot_2x = FALSE,
                               tempdir_mplus = getwd(), rep = NULL, time_max = 30){




    require("tidyverse")
    require("data.table")
    if (bayes_boot == TRUE){require("gtools")}
    require("MASS")
    require("mice")
    require("MplusAutomation")
    if (is.null(rep)){require("ids")}
    require("lpa.mi.src")

    # Add hash to temporary filename
    if(is.null(rep)){rep = ids::random_id(n = 1, bytes = 4)}

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
    names(wdf) = paste0("W",1:M)
    wobsdf = cbind(obsdf[,1:J_z], wdf)

    # Create A new file and write out the data.
    path_tmpdir = paste0(tempdir_mplus,"/temp-EMs-FMM-rep", rep, "-z",z)
    dir.create(path = path_tmpdir, recursive = TRUE)
    obj_Mplus = prepareMplusData(df = wobsdf,
                                 filename = paste0(path_tmpdir,"/EMs_FMM rep",rep," z", z,".dat"))

    #### Fit bootstrapped models until convergence with the assistance of a tracker (tracker_df) ####
      tracker_df = data.frame(path = path_tmpdir,
                              stem = paste0("abayes rep",rep, " z", z),
                              m = 1:M,
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
      iter = 0; flag = 0; Mvec = 1:M; tic = proc.time()
      while(flag == 0){
          for (i in Mvec){ #Write out .inp files for replicates that have not converged
              # Title:
              title_txt = c("TITLE:\n",
                            with(tracker_df,paste0(stem[i], " m",m[i],"\n")))
              # Data
              data_txt = c("DATA:\n",
                           paste0("FILE = ","'EMs_FMM rep",rep," z", z,".dat';\n") )
              #Variable
              variable_txt = c("VARIABLE:",
                               paste0("NAMES = \n", paste(" ",names(wobsdf), "\n", collapse = " "), ";\n"),
                               paste0("USEV = ", names(wobsdf)[1],"-",names(wobsdf)[J_z],";\n"),
                               paste0("CLASSES = c(",K_z,");\n"),
                               paste0("WEIGHT = W",tracker_df$m[i],";\n"),
                               "MISSING = .;\n")
              # Analysis
              starts0_i = tracker_df$starts0[i]
              analysis_txt = c("ANALYSIS:\n",
                               "TYPE = mixture;\n",
                               "ESTIMATOR = mlr;\n",
                               paste0("STARTS = ",starts0_i," ", floor(8*starts0_i/20),";\n"))
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
              fileConn<-file(with(tracker_df, paste0(path[i],"/",stem[i]," m",m[i],".inp")))
              writeLines(Mplus_txt, fileConn)
              close(fileConn)

          } # end for i in Mvec

          for(i in Mvec){ # Fit the models that have not reached convergence
            runModels(target = with(tracker_df, paste0(path[i],"/",stem[i], " m",m[i],".inp")),
                      showOutput = FALSE,
                      logFile = NULL)
          } # end for i in Mvec
          toc = proc.time()-tic

          # Assess convergence of the fit models
          tracker_df$converged[Mvec] = with(tracker_df,
                                            sapply(stem[Mvec],
                                                   function(x){
                                                     check_convergence(file = paste0(stem[x]," m",m[x],".out"),
                                                                       folder_wd = path_tmpdir,
                                                                       starts_txt = paste0(starts0[x]," ",floor(8*starts0[x]/20), ";"))
                                                   }
                                            )
          )
          Mvec = which(tracker_df$converged == FALSE)

          # Determine whether to continue within loop or exit
          if (length(Mvec)==0 & toc<(time_max*60)){
            flag = 1
          } else {
            tracker_df$starts0[Mvec] = 2*tracker_df$starts0[Mvec]
            print(tracker_df[Mvec, ])
          }


      }# end while


    # Read models in temporary directory
      mplus_list = readModels(target = path_tmpdir, what = c("parameters"), quiet = TRUE)

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
          cdraws_vec = get_pseudoc_draw(CPROBS_df = cprobs_m)

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
              X_mr = as.matrix(cbind(rep(1,N_r),Y_r[,idV]))
              E_mr = with(list_UV_mr, X_mr%*%beta)
              eps_mr = with(list_UV_mr, mvrnorm(n = N_r, mu = 0.*idU, S = S_eps))
              Yimp_m[ivec_rk,idU] = E_mr+eps_mr
            } #end for i
          } # end for k
        }
        return(Yimp_m %>% transform(subpop = obsdf$subpop, .imp = x))
      }

      list_pvs = lapply(X = 0:M, FUN = "doPV")

  # Clean up and spit out a 'mids' object
  system(paste0("rm -r '",gsub("/","\\\\",path_tmpdir),"' /s /q"))
  return(as.mids(rbindlist(list_pvs))
)

} #end function
