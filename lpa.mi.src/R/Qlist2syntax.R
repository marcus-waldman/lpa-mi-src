#' Qlist2syntax.R
#'
#' Qlist to svals syntax
#' @param Qlist- (list)
#' @param symbol - (character) Defaults to "(at)"
#' @return (character vector)
#' @export
#'
#' @examples

Qlist2syntax<-function(Qlist, symbol = "@"){

    K = length(Qlist$pi)
    J = nrow(Qlist$mu)

    overall_txt = ""
    if(K>1){
      overall_txt = c(overall_txt,
                      "%OVERALL%")
      for(kk in seq(1,K-1)){
        overall_txt = c(overall_txt,
                        paste0("[c#",kk,symbol,round(log(Qlist$pi[kk]/Qlist$pi[K]),5),"];")
        )
      }
    }

    modelC_txt = ""
    for (kk in 1:K){
      if (K>1){
        modelC_txt = c(modelC_txt,"",
                       paste0("%c#",kk,"%"))
      }
      for(jj in 1:J){
        modelC_txt = c(modelC_txt,
                       paste0("[Y",jj,symbol,round(Qlist$mu[jj,kk],5),"];"))
      }
      for (jj in 1:J){
        if(Qlist$S[jj,jj,kk]<1E-3){Qlist$S[jj,jj,kk] = Qlist$S[jj,jj,kk] + 1E-3}
        modelC_txt = c(modelC_txt,
                       paste0("Y",jj,symbol,round(Qlist$S[jj,jj,kk], 5), ";"))
      }
      for (j2 in seq(2,J)){
        for (j1 in seq(1,j2-1)){
          modelC_txt = c(modelC_txt,
                         paste0("Y",j1," WITH Y",j2,symbol,round(Qlist$S[j1,j2,kk],5),";"))
        }
      }
    }

    svals_txt = c(overall_txt, modelC_txt)

    return(svals_txt)
}
