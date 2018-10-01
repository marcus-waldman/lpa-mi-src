#'Get D statistic and D-statistic based difference ICs.
#'
#' Obtain the population parameters given the data
#' @param LL_df  (data.frame) with the number of p
#' @param M  (integer) Number of imputation
#' @param n (integer) Sample size
#' @return out_list  (list) with the following elements:
#'        (A) k_chosen - (data.frame) with the selected model based on ICs with D statistic
#'        (B) Rel_IC -  (data.frame) with relative information criteria (relative to first class IC values)
#'        (C) LL - (data.frame) LL_df updated with D statistic, D ICs, and ARIV estimate
#' @export
#' @examples
#' get_D_ICs(LL_df, M, n)

get_D_ICs<-function(LL_df, M, n){

    out_df = LL_df
    out_df = transform(out_df, Pool_Method = "D-statistic", D = NA, D_AIC = NA, D_BIC = NA)

    for (j in seq(1, nrow(LL_df)-1)){

      k_j = which(LL_df$kk == j)
      kp1_j = which(LL_df$kk == j+1)

      if (LL_df$Converged[kp1_j]==TRUE & LL_df$Converged[k_j]==TRUE){

          # From Enders(2010, p. 240-241). Using his notation
          LR_bar_j = -2*(LL_df$LL_free[k_j] - LL_df$LL_free[kp1_j])
          LR_bar_constrained_j =  -2*(LL_df$LL_fixed[k_j] - LL_df$LL_fixed[kp1_j])
          df_j = LL_df$params[kp1_j] - LL_df$params[k_j]
          ARIV_j = ((M+1)/(df_j*(M-1)))*(LR_bar_j - LR_bar_constrained_j)

          D_j = LR_bar_constrained_j /( df_j*(1+ARIV_j) )
          #Intermediary step since D not always positive
          #D_j = ifelse(D_j<0, 0, D_j)
          D_AIC_j =  df_j*(D_j - 2) #Specifically mentioned by Claeskins and Con..(2010)
          D_BIC_j = df_j*(D_j - log(n)) # Not specifically mentioned

          out_df$ARIV[j+1] = ARIV_j
          out_df$D[j+1] = D_j
          out_df$D_AIC[j+1] = -1*D_AIC_j # For intepretation of kp1 vs k rather k vs. kp1
          out_df$D_BIC[j+1] = -1*D_BIC_j
      } #end if

    } # end for j

    D_IC_df = data.frame(k = 1:nrow(LL_df),
                         R_AIC = rep(NA, nrow(LL_df)),
                         R_BIC = rep(NA, nrow(LL_df)))
    D_IC_df$R_AIC[1] <- D_IC_df$R_BIC[1] <- 0
    for (j in seq(1, nrow(LL_df)-1)){
      D_IC_df$R_AIC[j+1] = D_IC_df$R_AIC[j] + out_df$D_AIC[j+1]
      D_IC_df$R_BIC[j+1] = D_IC_df$R_BIC[j] + out_df$D_BIC[j+1]
    }

    selected_df= data.frame(t(apply(D_IC_df, 2, which.min)[-1]))
    names(selected_df) = c("D_AIC","D_BIC")
    selected_df = transform(selected_df,
                            kk = out_df$kk[1], z = out_df$z[1],
                            p = out_df$p[1], rep = out_df$rep[1],
                            pm = out_df$pm[1], pva = out_df$pva[1])

    out_list = list(k_chosen = selected_df, Rel_IC = D_IC_df, LL = out_df)
    return(out_list)

} # end function
