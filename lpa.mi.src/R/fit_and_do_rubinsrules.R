#' Fit mixture model to multiply imputed data and pool using Rubin's Rules (1987)
#'
#' @param dff_imputed (data.frame)
#' @param methods_list (list)
#' @param pop_params_z (list)
#' @param data_conditions (data.frame)
#' @param temp_wd_p_vec (vector) of characters
#' @param p (integer)
#' @param z (integer)
#' @param pm (integer)
#' @param pva (integer)
#' @return (list) out_list with (a) parameters, (b) ARIV, (c) Vtot, (B)
#' @export
#' @examples
#' fit_and_do_rubinrules(dff_imputed, methods_list, pop_params_z, data_conditions, temp_wd_p_vec, p, z, pm, pva)
#'

fit_and_do_rubinrules<-function(dff_imputed, methods_list, pop_params_z,
                                data_conditions, temp_wd_p_vec, p, z, pm, pva){

    require(Matrix)
    require(abind)
    require(data.table)

    M_pva = methods_list$args[[pva]]$m
    list_V_m = list(NULL)     #Q-by-Q-by_M_pva with within imputation asympototic covariance
    list_theta_m = list(NULL) #Q-by-M_pva with parameter column vectors
    list_params = list(NULL)
    mm = 0
    for (m in 1:M_pva){

      idx_pm_pva_m = which(endsWith(as.character(dff_imputed$folders),paste0("pm",pm,"/pva",pva)) &
                             endsWith(as.character(dff_imputed$files), paste0("imp_",m,".dat")))
      out_ftc_m = fit_til_convergence(p = p , z = z, temp_wd_p_vec = temp_wd_p_vec,
                                      dff_target = dff_imputed[idx_pm_pva_m,],
                                      target_wd = paste0(temp_wd_p_vec[p], "/Imputed data/pm",pm,"/pva",pva),
                                      pop_params_kk = pop_params_z,
                                      type_imputation = FALSE, m = m,
                                      readModels = TRUE, savedata = FALSE, output_txt = "tech3;")


      if (out_ftc_m$problem == FALSE & out_ftc_m$Rcond > 1E-8){


          mm = mm+1

          # Obtain TECH1  parameter and create placeholder for tech1_switched_m
          tech1_m = get_tech1(out_readModels = out_ftc_m$out_readModels ); names(tech1_m)[names(tech1_m) == "tech1"] = "tech1_from"
          tech1_switched_m = tech1_m; names(tech1_switched_m)[names(tech1_switched_m) == "tech1_from"] = "tech1_to"

          # Obtain Model parameters and merge TECH1 (tech1_from) information
          params_m = out_ftc_m$out_readModels$parameters$unstandardized
          params_m = merge(params_m, tech1_m, by = c("paramHeader", "param", "LatentClass"), sort = FALSE)

          # ODeal with label switching using the KL divergence discrepancy
          DMAT_m = get_KL_DMAT(params = params_m,
                               z = z,
                               data_conditions = data_conditions)
          out_switch = resolve_label_switch(mu_est_mat = NULL, S_est_array = NULL, z = z, data_conditions = data_conditions,
                                            parameters_df = params_m,
                                            DMAT = DMAT_m)
          params_switch_m = out_switch$parameters_df

          # Add in the TECH1 (tech1_to) information
          params_switch_m = merge(params_switch_m, tech1_switched_m, by = c("paramHeader", "param", "LatentClass"), sort = FALSE)
          params_switch_m = subset(params_switch_m, !is.na(tech1_from)) %>% arrange(tech1_from)

          # Obtain within-imputation asymptotic covariance and deal with label switching
          hi = out_ftc_m$out_readModels$tech3$paramCov; hi[is.na(hi)] = 0
          thi = t(hi); diag(thi) = 0;
          V_m = hi+thi
          params_switch_m = params_switch_m %>% arrange(tech1_to)
          V_m_switched = V_m[params_switch_m$tech1_from, params_switch_m$tech1_from]


          # Save parameters and asympoptotic covariance estimates in a matrix
          list_params[[mm]] = params_switch_m
          list_V_m[[mm]] = as.matrix(V_m_switched) #Q-by-Q-by_M_pva with within imputation asympototic covariance
          list_theta_m[[mm]] = params_switch_m$est #Q-by-M_pva with parameter column vectors


      } # end if out_ftc_com$problem == FALSE

    } # end for m


  MM_pva = length(list_theta_m)
  if (MM_pva>=5){

  #print("line 77")
      V_m_array = array(dim = c(nrow(V_m), ncol(V_m), MM_pva)) + NA
      theta_m_mat = mat.or.vec(nr = nrow(V_m), nc = MM_pva)
      for (mm in 1:MM_pva){
        V_m_array[,,mm] = as.matrix(list_V_m[[mm]])
        theta_m_mat[,mm] = as.numeric(list_theta_m[[mm]])
      }


  #print("line 86")
      Vbar = apply(V_m_array, c(1,2),"mean")
      t_bar = apply(theta_m_mat, 1, "mean")


  #print("line 91")
      # Get B
      sumB = mat.or.vec(nr = nrow(Vbar), nc = ncol(Vbar))
      for (mm in 1:MM_pva){
        dt_m = as.numeric(theta_m_mat[,mm])-t_bar
        sumB = sumB + dt_m%*%t(dt_m)
      }
      B = sumB/(MM_pva-1)


  #print("line 101")
      trBV = sum(diag(B%*%solve(Vbar)))
      Q = nrow(Vbar)
      ARIV = (1+1/M_pva)*trBV/Q


      Vtot = Vbar + (1+1/MM_pva)*B

      parameters_df = params_switch_m[,c("paramHeader","param","LatentClass")] %>%
                        transform(est = t_bar, se = sqrt(diag(Vtot))) %>%
                        transform(est_se = est/se) %>%
                        transform(pval = 1-pnorm(abs(est_se)))

  #print("line 114")
      out_list = list(problem = FALSE,
                      MM_pva = MM_pva,
                      parameters = parameters_df,
                      ARIV = ARIV,
                      Vtot = Vtot,
                      Vbar = Vbar,
                      B = B,
                      V_m_list = list_V_m,
                      theta_m_list = list_theta_m,
                      params_m_list = list_params)

  } else {
    warning("Total imputated data sets that converged was less than 5.")
    out_list = list(problem = TRUE,
                    MM_pva = MM_pva,
                    parameters = NULL,
                    ARIV = NULL,
                    Vtot = NULL,
                    Vbar = NULL,
                    B = NULL,
                    V_m_list = NULL,
                    theta_m_list = NULL,
                    params_m_list = NULL)

  }

  #print("line 124")
      return(out_list)

}
