#' Rotation matrix
#'
#' Construct a rotation matrix that rotates data in R2 (or R3, eventually)
#' @param t - (numeric) Rotation angle in radians. Defaults to pi/4.
#' @param dim (integer) Rank of rotation matrix. Defaults to 2 for rotating in R2.
#' @return (Matrix) A rotation matrix of size dim-by-dim.
#' @export
#' @examples
#' A = get_rotation_matrix(t = pi/3, dim = 2)

get_rotation_matrix<-function(t = pi/4,dim = 2){
  #get the rotation matrix, ROTATE
  # Inputs - (a) angle (in radians) and (b) dimension (2 or 3)
  
  #t = pi/4
  #dim = 3
  
  if (dim<2){
    stop("Must rotate in more than two dimensions")
  }
  if (abs(t)==0 | abs(t)>pi/2){
    stop("Rotation angle in Radians must be between 0<theta<pi/2")
  }
  
  A = matrix(c(cos(t), -sin(t), 
               sin(t), cos(t)), nrow = 2, byrow = TRUE)
  
  if (dim>2){
    tmp = diag(dim)
    for (j in seq(1,dim, by = 1)){
#      print(j)
      if (j<dim){
        jvec_i = c(j,j+1)
      } else {
        jvec_i = c(1,j)
      }
      A_i = diag(dim)
      A_i[jvec_i,jvec_i] = A
      tmp = A_i%*%tmp
    }
    ROTATE = tmp
  } else {
    ROTATE = A
  }
  return(ROTATE)
}
