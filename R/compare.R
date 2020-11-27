#' Compare predictions and labels
#'
#' @param predictions 
#' @param labels 
#' @param partition 
#' @param names 
#'
#' @return
#' @export
#'
#' @examples
compare <- function(predictions,
                    labels, 
                    partition = TRUE,
                    names = c(specificity=0, sensitivity=1)){
  if(!(is.matrix(predictions) | is.data.frame(predictions)) ){
    predictions <- matrix(predictions)
  }
  stopifnot(is.numeric(labels))
  stopifnot(nrow(predictions) == length(labels))
  
  comp <- matrix(rep(labels, ncol(predictions)), byrow=FALSE, ncol=ncol(predictions))
  out <- as.data.frame(1 * (predictions == comp))
  
  if(partition){
    out <- split(out, labels)
    if(!is.null(names)){
      names(out) <- names(names)[match(names(out), names)]
    }
  }
  return(out)
}