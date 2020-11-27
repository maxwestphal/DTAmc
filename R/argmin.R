#' argmin
#'
#' @param x 
#' @param rdm 
#'
#' @return
#' @export
#'
#' @examples
argmin <- function(x, rdm=FALSE){
  am <- which(x == min(x))
  # deterministic output, if required, otherwise randomize in case of ties
  if(!rdm){ 
    return(min(am))
  }
  if(length(am) == 1){
    return(am)
  }else{
    return(sample(am, 1))
  }
}

#' p argmin
#'
#' @param ... 
#' @param args 
#' @param rdm 
#'
#' @return
#' @export
#'
#' @examples
pargmin <- function(..., args=list(), rdm=FALSE){
  stopifnot(is.list(args))
  args <- c(list(...), args)
  stopifnot(do.call(all.equal, lapply(args, length)))
  apply(as.data.frame(args), 1, argmin, rdm=rdm)
}

#' argmax
#'
#' @param x 
#' @param rdm 
#'
#' @return
#' @export
#'
#' @examples
argmax <- function(x, rdm=FALSE) {
  argmin(-x, rdm)
}

#' p argmax
#'
#' @param ... 
#' @param args 
#' @param rdm 
#'
#' @return
#' @export
#'
#' @examples
pargmax <- function(..., args=list(), rdm=FALSE){
  pargmin(-x, args, rdm)
}




