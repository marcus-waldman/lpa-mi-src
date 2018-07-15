#' Mean locations of spheres
#'
#' Returns locations of mean vectors for spheres all separated by unit Euclidean distance
#' @param J  (integer) Dimension of space R^J
#' @return (matrix) Matrix of mean column vectors for each unit sphere.
#' @export
#' @examples
#' get_mu_all(J)

get_mu_all<-function(J){

  require(stats)  #(>= 3.5.0)

   loss<-function(y, X){
     loss = 0
     for(j in seq(1,ncol(X))){
       loss = loss + ((sum((y-X[,j])^2)-1))^2
     }
     return(loss)
   }

  K_all = J+1;
  mu_all =  mat.or.vec(nr = J, nc = K_all)
  mu_all[1,2] = 1

  if (J>2){
    for(k in 3:K_all){
      X = mu_all[1:(k-1),1:(k-1)]
      mu_start = rep(0,(k-1)); mu_start[k-1] = 1;
      obj_optim = optim(mu_start, fn = loss, gr = NULL, X, method = "L-BFGS", lower = rep(0,k-1), control = list(ndeps = rep(sqrt(.Machine$double.eps),k-1)))
      mu_all[1:(k-1),k] = obj_optim$par
    }
  }

  return(mu_all)
}
