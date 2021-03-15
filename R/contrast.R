#' define a contrast (matrix) to specify exact hypothesis system
#'
#' @param type c("Raw", "Dunnett", "Tukey")
#' @param base 
#'
#' @return
#' @export
#'
#' @examples
define_contrast <- function(type = c("Raw", "Dunnett", "Tukey"), base=1){
  type <- match.arg(type)
  stopifnot(base %% 1 == 0 & base >= 1)

  
  fun <- function(data=NULL, n=colnames(data[[1]]), G=length(data)){
    stopifnot(!is.null(data) | !is.null(n))
    if(is.numeric(base)){
      stopifnot(base %in% 1:length(n))
      b <- base
      base <- which(base == n)
    }
    if(is.character(base)){
      stopifnot(base %in% n)
      b <- which(base == n)
    }
    if(!is.null(n)){
      x <- rep(1, length(n))
      names(x) <- n
    }
    if(type == "Raw"){
      K <- diag(x)
      rownames(K) <- colnames(K) <- names(x)
    }else{
      K <- multcomp::contrMat(x, type = type, base=base)
    }
    class(K) <- "matrix"
    attributes(x) <- NULL
    return(replicate(G, K, simplify=FALSE))
  }
  class(fun) <- append(class(fun), "contrastFun")
  attr(fun, "type") <- type
  attr(fun, "base") <- ifelse(type %in% c("Dunnett"), base, NA)
  return(fun)
}

preproc_contrast <- function(x, data){
  G <- length(data)
  m <- ncol(data[[1]])
  if("contrastFun" %in% class(x)){
    return(x)
  }
  if(is.list(x)){
    stopifnot(length(x) == G)
    stopifnot(all(sapply(x, is.matrix)))
    stopifnot(all(sapply(x, ncol) == m))
    out <- x
  }
  if(is.matrix(x)){
    stopifnot(ncol(x) == m)
    out <- replicate(length(data), x)
  }
  fun <- function(...){
    out
  }
  return(fun)
}


# define_contrast("Raw")(n=paste0("rule", 1:3), G=2)
# define_contrast("Dunnet", 1)(n=paste0("rule", 1:3), 2)
# define_contrast("Tukey")(n=paste0("rule", 1:3))

#data <- sample_data_lfc(m=3)
#preproc_contrast(define_contrast()(data=data), data=data)()[[1]]


# n <- c(10,20,30,40)
# names(n) <- paste("group", 1:4, sep="")
# 
# K <- multcomp::contrMat(n, "Dunnet", base=2) 
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

