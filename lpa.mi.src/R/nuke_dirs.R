#' Permanently remove folers and subfolders
#'
#' @param directories  (vector) Character strings of directories
#' @return
#' @export
#' @examples

nuke_dirs<-function(directories){

  message("Permanently removing directories... ")
  for (i in 1:length(directories)){
    if(!endsWith(directories[i],"/")){
      nukedir_i = paste0(directories[i],"/")
    } else {
      nukedir_i = directories[i]
    }
    unlink(nukedir_i, recursive = TRUE)
    message(paste0("   ", nukedir_i," successfully nuked."))
  }

}
