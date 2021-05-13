#' @export
print.DTAmcResults <- function(x, info=FALSE){
  message("Diagnostic Test Accuracy:")
  if(!info){
    n <- names(x)
    attributes(x) <- NULL
    names(x) <- n
  }
  print.default(x)
}