#' Index from results list
#'
#' @param x  (character vector) with names from a results list 
#' @param head (character) how x should starte (e.g. "list_complete")
#' @param p Optional (integer) with processor number 
#' @param z Optional (integer) with data condition
#' @param rep Optional (integer) with replication
#' @param pm Optional (integer) with percent missing
#' @param pva Optional (integer) with imputation algorithm
#' @param sub Substitute periods if there are spaces the elements of x. Done for correct search. Defaults to TRUE. 
#'
#' @return vector of integers with indices in x that match.
#' @export
#'
#' @examples
ix<-function(x, head, p=NULL,z=NULL,rep=NULL,pm=NULL,pva=NULL, sub = TRUE){

  
  if(!is.null(pm)&head=="list_complete"){stop("head = list_complete, but pm specified")}
  if(!is.null(pva)&head=="list_complete"){stop("head = list_complete, but pva specified")}
  if(!is.null(pm)&head=="list_observed"){stop("head = list_observed, but pm specified")}
  if(!is.null(pva)&head=="list_observed"){stop("head = list_observed, but pva specified")}
  ix_p <- ix_z <- ix_r <- ix_v <- ix_a<-NULL
  
  if(sub==TRUE){x=gsub(pattern = " ", replacement = ".", x = x)}
  
  ix = which(startsWith(x,head))
  
  ix_p = if(!is.null(p)){ix_p=grep(paste0("p",p,"."), x, fixed=T)}
  if(!is.null(ix_p)){ix = intersect(ix,ix_p)}
  
  ix_z = if(!is.null(z)){grep(paste0("z",z,"."), x, fixed=T)}
  if(!is.null(ix_z)){ix = intersect(ix,ix_z)}
  
  ix_r = if(!is.null(rep)){grep(paste0("rep",rep,"."), x, fixed=T)}
  if(!is.null(ix_r)){ix = intersect(ix,ix_r)}
  
  ix_m = if(!is.null(pm)){grep(paste0("pm",pm, "."), x, fixed=T)}
  if(!is.null(ix_m)){ix = intersect(ix,ix_m)}
  
  ix_a = if(!is.null(pva)){grep(paste0("pva",pva,"."), x, fixed=T)}
  if(!is.null(ix_a)){ix = intersect(ix,ix_a)}
  
  
  return(ix)
  
}