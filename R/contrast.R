#' define a contrast (matrix) to specify exact hypothesis system
#'
#' @param type c("Raw", "Dunnett", "Tukey")
#' @param comparator either integer (index of comparator, cannot be larger than ncol(data[[1]])) 
#' or character (name of rule to use for )
#'
#' @return
#' @export
#'
#' @examples
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
  class(fun) <- append(class(fun), "contrastFun")
  return(fun)
}




## TODO: DELETE OLD STUFF?
# preproc_contrast <- function(x, data){
#   G <- length(data)
#   m <- ncol(data[[1]])
#   if("contrastFun" %in% class(x)){
#     return(x)
#   }
#   if(is.list(x)){
#     stopifnot(length(x) == G)
#     stopifnot(all(sapply(x, is.matrix)))
#     stopifnot(all(sapply(x, ncol) == m))
#     out <- x
#   }
#   if(is.matrix(x)){
#     stopifnot(ncol(x) == m)
#     out <- replicate(length(data), x)
#   }
#   fun <- function(...){
#     out
#   }
#   return(fun)
# }


# define_contrast("raw")(n=paste0("rule", 1:3))
# define_contrast("Dunnet", 1)(n=paste0("rule", 1:3), 2)
# define_contrast("Tukey")(n=paste0("rule", 1:3))

#data <- sample_data_lfc(m=3)
#preproc_contrast(define_contrast()(data=data), data=data)()[[1]]


# n <- c(10,20,30,40)
# names(n) <- paste("group", 1:4, sep="")
# 
# K <- multcomp::contrMat(n, "Dunnet", comparator=2) 
# class(K) <- "matrix"
# 
# K
# 
# multcomp::contrMat(n, type = "Tukey")
# multcomp::contrMat(n, type = "Sequen")
# multcomp::contrMat(n, type = "AVE")
# multcomp::contrMat(n, type = "Changepoint")
# multcomp::contrMat(n, type = "Williams")
# multcomp::contrMat(n, type = "Marcus")
# multcomp::contrMat(n, type = "McDermott")

