#'Get difference in ICs from enumeration summaries
#'
#' Obtain the population parameters given the data
#' @param result_df  (data.frame) giving summary results during enumeration with variables: (a) Parameters, (b)AIC, and (c) BIC
#' @return (data.frame) with difference ICs appended to results_df input file
#' @export
#' @examples
#' get_diffICs(LL_df)

get_diffICs<-function(result_df, token = ""){

    result_df = transform(result_df, df = NA, D = NA, D_AIC = NA, D_BIC = NA)
    for (j in 2:nrow(result_df)){


      if (result_df$Converged[j]==TRUE & result_df$Converged[j-1]==TRUE){

        df_j = result_df$Parameters[j] - result_df$Parameters[j-1]
        D_j = 2*(result_df$LL[j] - result_df$LL[j-1])/df_j
        D_AIC_j = result_df$AIC[j]-result_df$AIC[j-1]
        D_BIC_j = result_df$BIC[j]-result_df$BIC[j-1]

        result_df$df[j] = df_j
        result_df$D[j] = D_j
        result_df$D_AIC[j] = D_AIC_j
        result_df$D_BIC[j] = D_BIC_j

      }

    }

    names(result_df)[seq(ncol(result_df)-3, ncol(result_df))] = paste0(names(result_df)[seq(ncol(result_df)-3, ncol(result_df))], token)

    return(result_df)

}
