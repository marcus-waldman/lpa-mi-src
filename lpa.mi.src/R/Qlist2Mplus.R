#' Convert estimated Qlist to Mplus format
#'
#' @param Qlist (list)
#' @param z (integer)
#' @param data_conditions (data.frame)
#'
#' @return (data.frame) with estimates in Mplus parameter format
#' @export
#'
#' @examples
Qlist2Mplus<-function(Qlist, z, data_conditions){

    require(tidyverse)

    FMM_params = get_FMM_params(z = z, data_conditions = data_conditions)
    Nparams = with(FMM_params, (K_z-1) + K_z*(2*J_Y_z + choose(J_Y_z,2)))


    means_df = with(FMM_params,
                    data.frame(paramHeader = rep("Means", K_z*J_Y_z),
                               param = rep(paste0("Y",1:J_Y_z),K_z),
                               LatentClass = as.character(sort(rep(1:K_z,J_Y_z)))
                    )) %>% transform(est = as.vector(Qlist$mu))


    var_df = with(FMM_params,
                  data.frame(paramHeader = rep("Variances", K_z*J_Y_z),
                             param = rep(paste0("Y",1:J_Y_z),K_z),
                             LatentClass = as.character(sort(rep(1:K_z,J_Y_z)))
                  )) %>% transform(est =  as.numeric(apply(Qlist$S,
                                                            MARGIN = 3,
                                                            FUN = function(x){diag(x)})) )




    vec1 = NULL
    for(j1 in seq(1,FMM_params$J_Y_z-1)){
      vec1 = c(vec1, with(FMM_params, rep(j1,J_Y_z-j1)))
    }
    vec1 = with(FMM_params, rep(vec1, K_z))

    vec2 = NULL
    for(j2 in seq(2,FMM_params$J_Y_z)){
      vec2 = c(vec2, with(FMM_params, seq(j2,J_Y_z)))
    }
    vec2 = with(FMM_params, rep(vec2, K_z))




    cov_df = with(FMM_params,
                  data.frame(paramHeader = paste0("Y",vec1,".WITH"),
                             param = paste0("Y",vec2),
                             LatentClass = as.character(sort(rep(1:K_z,choose(K_z,2))))
                  ))%>% transform(est = as.vector(apply(Qlist$S,MARGIN = 3,
                                                          FUN = function(x){x[lower.tri(x, diag = FALSE)]}
                                                          ))) %>%
                  arrange(LatentClass)

    pivec = as.numeric(Qlist$pi)
    lor_vec = with(FMM_params, log(pivec[-K_z]/pivec[K_z]))
    lor_df = with(FMM_params,
                  data.frame(paramHeader = rep("Means", K_z-1),
                             param = paste0("C#",seq(1,K_z-1)),
                             LatentClass = as.character(rep("Categorical.Latent.Variables", K_z-1))
                  )) %>% transform(est = lor_vec)

    mplus_df = rbind(means_df, var_df, cov_df, lor_df)


    return(mplus_df)
}
