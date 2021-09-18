#' Define a contrast (matrix) to specify exact hypothesis system
#'
#' @param type character, either "Raw", "Dunnett" or "Tukey")
#' @param comparator either integer (index of comparator) or character (name of comparator)
#'
#' @return \code{DTAmc_contrast} object to be passed to \code{\link{study_dta}}
#' @export
#'
#' @examples 
#' define_contrast("Dunnett", 1)
#' @importFrom multcomp contrMat
define_contrast <- function(type = c("raw", "Dunnett", "Tukey"), comparator=NA){
  type <- match.arg(type)
  
  fun <- function(data=NULL,
                  n=colnames(data[[1]])
                  ){
    stopifnot(!is.null(data) | !is.null(n))
    
    x <- rep(1, length(n))
    names(x) <- n
  
    if(type == "raw"){
      K <- diag(x)
      rownames(K) <- colnames(K) <- names(x)
    }
    if(type == "Dunnett"){
      stopifnot(comparator %% 1 == 0 & comparator >= 1)
      if(is.numeric(comparator)){
        stopifnot(comparator %in% 1:length(n))
        b <- comparator
        comparator <- n[b]
      }
      if(is.character(comparator)){
        stopifnot(comparator %in% n)
        b <- which(comparator == n)
      }
      K <- multcomp::contrMat(x, type = type, base=b)
    }
    if(type == "Tukey"){
      K <- multcomp::contrMat(x, type = type, base=1)
    }
    
    class(K) <- "matrix"
    attr(K, "type") <- type
    attr(K, "comparator") <- ifelse(type %in% c("Dunnett"), comparator, NA)
    return(K)
  }
  class(fun) <- append(class(fun), "DTAmc_contrast")
  return(fun)
}
