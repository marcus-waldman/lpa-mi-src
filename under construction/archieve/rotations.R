#' Check convergence
#'
#' Scans an .out file from an Mplus model and looks for signs of convergence.
#' @param t - (numeric) Rotation angle in radians. Defaults to pi/4.
#' @param dim (integer) Rank of rotation matrix. Defaults to 2 for rotating in R2.
#' @return (Matrix) A rotation matrix of size dim-by-dim.
#' @export
#' @examples

rotation_matrix<-function(t = pi/4,dim = 2){
  #get the rotation matrix, ROTATE
  # Inputs - (a) angle (in radians) and (b) dimension (2 or 3)
  
  #t = pi/4
  #dim = 3
  
  if (abs(t)==0 | abs(t)>pi/2){
    stop("Rotation angle in Radians must be between 0<theta<pi/2")
  }
  
  if (dim == 2){
    ROTATE = matrix(c(cos(t), -sin(t), 
                      sin(t), cos(t)), nrow = dim, byrow = TRUE)
  } else if (dim == 3){
    stop("Rotation matrix can only currently be produced in R2")
    # Rodriques' formula
    ROTATE = matrix(c( 1, 1-cos(t)-sin(t), sin(t) + 1-cos(t), 
                       sin(t)+1-cos(t), 1, -sin(t) + 1-cos(t), 
                       -sin(t)+1-cos(t), sin(t)+1-cos(t), 1), 
                    nrow = dim, byrow = TRUE)
  } else {
    stop("Rotation matrix can only currently be produced in R2")
  }
  
  return(ROTATE)
}
