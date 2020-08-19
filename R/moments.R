data_moments <- function(d) {
  list(n = nrow(d), moments = t(d) %*% d)
}

prior_moments <- function(m, pars = c(2, 1, 1 / 2)) {
  A <- matrix(pars[3], m, m)
  diag(A) <- rep(pars[2], m)
  list(n = pars[1], moments = A)
}

add_moments <- function(pmom, dmom){
  return(mapply("+", pmom, dmom, SIMPLIFY = FALSE))
}

mom2est <- function(mom){
  return(diag(mom[[2]]) / mom[[1]])
}

mom2cov <- function(mom){
  n <- mom[[1]]
  A <- mom[[2]]
  a <- diag(A)
  return((n*A - (a %*% t(a))) / (n^2) / (n+1))
}

dat2est <- function(dat){
  colMeans(dat)
}

dat2cov <- function(dat){
  cov(dat)/nrow(dat)
}

dat2var <- function(dat){
  apply(dat, 2, var)/nrow(dat)
}

dat2se <- function(dat){
  sqrt(dat2var(dat))
}

cov2var <- function(cov){
  diag(cov)
}

cov2se <- function(cov){
  sqrt(cov2var(cov))
}

# data_moments2 <- function(d){
#   M <- matrix(NA, ncol(d), ncol(d))
#   for(j1 in 1:ncol(d)){
#     for(j2 in j1:ncol(d)){
#       M[j1, j2] <- M[j2, j1] <- sum(d[, j1] + d[, j2]==2)
#     }
#   }
#   return(M)
# }