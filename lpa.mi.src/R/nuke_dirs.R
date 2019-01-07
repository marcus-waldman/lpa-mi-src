#' Permanently remove folers and subfolders
#'
#' @param directories  (vector) Character strings of directories
#' @param ...  other arguments passed to the unlink() function.
#' @return
#' @export
#' @examples

nuke_dirs<-function(directories,...){

  message("Permanently removing directories... ")
  for (i in 1:length(directories)){
    if(!endsWith(directories[i],"/")){
      nukedir_i = paste0(directories[i],"/")
    } else {
      nukedir_i = directories[i]
    }
    do.call(what = "unlink", args = list(x = nukedir_i, ...))
    message(paste0("   ", nukedir_i," successfully nuked."))
  }

}
