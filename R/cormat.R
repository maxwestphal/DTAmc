#' Create an equicorrelation matrix
#'
#'
#' @param m integer, dimension
#' @param rho numeric, correlation parameter in (0,1)
#' @param d binary vector of length m, whereby TRUE/FALSE (alternatively 1/0)
#' indicate active/inactive components of underlying random vector. 
#'
#' @return \eqn{R_{ij} = \rho, i\neq j}
#' @export
#'
#' @examples
cormat_equi <- function(m, rho, d=TRUE){
  R <- matrix(rho, m, m)
  stopifnot(length(d) %in% c(1, m))
  R <- diag(d, m) %*% R %*% diag(d, m)
  diag(R) <- 1
  return(R)
}

#' Create an AR(1) correlation matrix
#'
#' @param m integer, dimension
#' @param rho numeric, correlation parameter in (0,1)
#' @param d binary vector of length m, whereby TRUE/FALSE (alternatively 1/0)
#' indicate active/inactive components of underlying random vector. 
#'
#' @return \eqn{R_{ij} = \rho^{|i-j|}}
#' 
#' @export
#'
#' @examples
cormat_ar1 <- function(m, rho, d=TRUE){
  M <- matrix(rho, m, m)
  R <- M^(abs(col(M) - row(M)))
  stopifnot(length(d) %in% c(1, m))
  R <- diag(d, m) %*% R %*% diag(d, m)
  diag(R) <- 1
  return(R)
}


## TODO: finish documentation
## TODO: need stuff below??

### FUNCTION: corr2R
# corr2R <- function(s, args=list()){
#   l <- string2list(s)
#   a <- c(l[-1], args)
#   R <- do.call(paste0("cormat_", l$type), list(a=a)) #TODO: correct?
#   return(R)
# }

### FUNCTION: string2list
# string2list <- function(s, sep1="_", sep2="=", convert=function(x)x){
#   l <- strsplit(s, split=sep1)[[1]]
#   ll <- strsplit(l, split=sep2)
#   a <- lapply(ll, function(x)x[2])
#   names(a) <- lapply(ll, function(x)x[1])
#   a <- lapply(a, convert)
#   return(a)
# }

### FUNCTION: corrmat_equi
# cormat_equi <- function(a){
#   a <- lapply(a, as.numeric)
#   R <- matrix(a$rho, a$S, a$S)
#   if(!is.null(a$d)){
#     R <- diag(a$d, length(a$d)) %*% R %*% diag(a$d, length(a$d))
#   }
#   diag(R) <- 1
#   return(R)
# }
# 
# ### FUNCTION: corrmat_ak
# cormat_ak <- function(a){
#   a <- lapply(a, as.numeric)
#   M <- matrix(a$rho, a$S, a$S)
#   R <- M^(abs(col(M) - row(M)))
#   if(!is.null(a$d)){
#     R <- diag(a$d, length(a$d)) %*% R %*% diag(a$d, length(a$d))
#   }
#   diag(R) <- 1
#   return(R)
# }