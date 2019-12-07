#' Get tech1 estimates identifiers
#'
#' @param out_readModels Object from MplusAutomations readModels() function.
#' @return tech1_df (data.frame) with tech1 identifier for Mplus parameters
#' @export
#' @examples
#' tech1_df = get_tech1(out_readModels)


get_tech1<-function(out_readModels){

      tech1_out = out_readModels$tech1$parameterSpecification
      K_out = ifelse( length(out_readModels$class_counts)>0 , nrow(out_readModels$class_counts[[1]]), 1)
      params_unstd = out_readModels$parameters$unstandardized
      if (K_out==1){params_unstd = transform(params_unstd, LatentClass = 1)}
      tech1_df = params_unstd[,c("paramHeader", "param", "LatentClass")] %>% transform(tech1 = NA)

      for(k in 1:K_out){
        tech1_k = if(K_out>1){ tech1_k = tech1_out[[k]] }else{ tech1_k=tech1_out }
        nu_k = data.frame(tech1_k$nu)
        theta_k = tech1_k$theta

        J_k = length(nu_k)
          for(j in 1:J_k){
            inds_jk = which( tech1_df$paramHeader=="Means" & tech1_df$LatentClass == k & tech1_df$param == names(nu_k)[j] )
            tech1_df$tech1[inds_jk] = as.numeric(as.vector(nu_k)[j])
          }

          for(j in 1:J_k){
            inds_jk = which( tech1_df$paramHeader=="Variances" & tech1_df$LatentClass == k & tech1_df$param == names(nu_k)[j] )
            tech1_df$tech1[inds_jk] = as.numeric(theta_k[j,j])
          }

        for(j1 in 2:J_k){
          for(j2 in seq(1,j1)){
            inds_jk = which( tech1_df$paramHeader == paste0(names(nu_k)[j2],".WITH") & tech1_df$LatentClass == k & tech1_df$param == names(nu_k)[j1] )
            tech1_df$tech1[inds_jk] = as.numeric(theta_k[j1,j2])
          }
        }

      }

      clog_vec = as.vector(tech1_out$LATENT.CLASS.REGRESSION.MODEL.PART$alpha.c)
      for (k in seq(1,K_out-1)){
        inds_jk = which(tech1_df$paramHeader == "Means" & tech1_df$param==paste0("C#",k))
        tech1_df$tech1[inds_jk] = as.numeric(clog_vec[k])
      }

      tech1_df$tech1[tech1_df$tech1==0] = NA

      return(tech1_df)

}
