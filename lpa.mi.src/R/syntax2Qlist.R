#' Convert Mplus MODEL: syntax to a Qlist
#'
#' @param syntax_x (character)
#' @param Plist_x (list) Optional. Defaults to NULL. If not specified, then extracts using get_Plist
#' @param C_letter (character) Optional. Defaults to "C".
#' @param Y_letter (character) Optional. Defaults to "Y"
#' @param fix_S (logical) Defaults to true.
#'
#' @return (list)
#' @export
#'
#' @examples
syntax2Qlist<-function(syntax_x, Plist_x ,C_letter = "C",Y_letter="Y",fix_S = TRUE){
  library(stringr)
  library(lpa.mi.src)

  # Turn everything into upper-case
  syntax_x = toupper(syntax_x); C_letter = toupper(C_letter); Y_letter = toupper(Y_letter);

  # Identify the indices for each class block
  th = c(which(str_detect(string = syntax_x, pattern = paste0("%",C_letter,"#"))),
         length(syntax_x)+1)
  K_x = length(th)-1 #Number of classes in the syntax

  # Segment the syntax
  template_list = list(class = NA, syntax = NA, cov = NA, mean = NA, var = NA)
  syntax_list<-lapply(X = 1:K_x, FUN = function(k){
    list_k = template_list;
    syntax_k = toupper(str_trim(syntax_x[seq(th[k], th[k+1]-1)]));

    idx_with_k = which(str_detect(string = syntax_k, pattern = "WITH"));
    if (length(idx_with_k)==0){stop("WITH statments not found and required to be explicitly specified in syntax_x")}

    idx_mean_k = which(str_detect(string = syntax_k, pattern = "\\["));

    idx_var_k = setdiff(which(str_detect(string = syntax_k, pattern = Y_letter)),
                        c(idx_with_k, idx_mean_k));
    if (length(idx_with_k)==0){stop("Indicator variances not found. Check the indicator variable letter is specified correctly.")}

    list_k$class = k
    list_k$syntax = syntax_k;
    list_k$cov = syntax_k[idx_with_k];
    list_k$mean = syntax_k[idx_mean_k];
    list_k$var = syntax_k[idx_var_k]

    return(list_k);
  })


  # Generate a template to store the estimates from the syntax above
  Qlist_x = Plist_x
  Qlist_x$mu = NA*Qlist_x$mu
  Qlist_x$S = NA*Qlist_x$S
  Qlist_x$pi = NA*Qlist_x$pi

  J_x = nrow(Qlist_x$S[,,1])

  # GEt the probabilities
  trimmed_x = str_trim(syntax_x)
  trimmed_x = str_remove(string=trimmed_x, pattern = " ")
  idx_clog = which(str_detect(string = trimmed_x, pattern = paste0("\\[",C_letter,"#")))
  tmp = trimmed_x[idx_clog]
  tmp = str_trim(tmp)
  for (k in seq(1,K_x-1)){
    tmp = str_remove(string = tmp, pattern = paste0(C_letter,"#",k,"\\*"))
  }
  tmp = str_remove(string = tmp, pattern = c("\\["))
  tmp = str_remove(string = tmp, pattern = c("\\];"))
  tmp = str_trim(tmp)
  Qlist_x$pi = gamma2pi(as.numeric(tmp))
  rm(tmp)

  for (k in 1:K_x){

    list_k = syntax_list[[k]]



    # means
    tmp = list_k$mean
    tmp = str_remove(string = tmp, pattern = c("\\["))
    tmp = str_remove(string = tmp, pattern = c("\\];"))
    for (j in 1:length(tmp)){
      tmp = str_remove(tmp, paste0("Y",j))
    }
    tmp = str_remove(string = tmp, pattern = c("\\*"))
    tmp = str_remove(string = tmp, pattern = c("\\@"))
    tmp = str_trim(tmp)
    Qlist_x$mu[,k] = as.numeric(tmp)
    rm(tmp)


    # cov
    tmp = list_k$cov
    tmp = str_remove(string = tmp, pattern = c("WITH"))
    tmp = str_remove(string = tmp, pattern = c(";"))
    for (j in 1:length(tmp)){
      tmp = str_remove(tmp, paste0("Y",j))
    }
    tmp = str_remove(string = tmp, pattern = c("\\*"))
    tmp = str_remove(string = tmp, pattern = c("\\@"))
    tmp = str_trim(tmp)
    a = 0
    for(i in seq(1,J_x-1)){
      for (j in seq(i+1,J_x)){
        a = a+1
        Qlist_x$S[i,j,k] <- Qlist_x$S[j,i,k]<-as.numeric(tmp[a])
      }
    }
    rm(tmp); rm(a);

    # var
    tmp = list_k$var
    tmp = str_remove(string = tmp, pattern = c(";"))
    for (j in 1:length(tmp)){
      tmp = str_remove(tmp, paste0(Y_letter,j))
    }
    tmp = str_remove(string = tmp, pattern = c("\\*"))
    tmp = str_remove(string = tmp, pattern = c("\\@"))
    tmp = str_trim(tmp)
    for (j in 1:length(tmp)){
      Qlist_x$S[j,j,k] = as.numeric(tmp[j])
    }
    rm(tmp)


  }

  # Update the poor conditioning
  if(fix_S==TRUE){
    poor_rcond = TRUE
    while(poor_rcond==TRUE){
      list_check_rcond = check_rcond_S(Plist = Qlist_x)
      Qlist_x = list_check_rcond$Plist
      poor_rcond = list_check_rcond$poor_rcond
    }
  }


  return(Qlist_x)
}
