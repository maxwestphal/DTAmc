## data to moments calculations
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


## calculations based on moments 
mom2est <- function(mom){
  return(diag(mom[[2]]) / mom[[1]])
}

mom2cov <- function(mom){
  n <- mom[[1]]
  A <- mom[[2]]
  a <- diag(A)
  return((n*A - (a %*% t(a))) / (n^2) / (n+1)) ## TODO: check correctness (PUB3)
}

cov2var <- function(cov){
  diag(cov)
}

cov2se <- function(cov){
  sqrt(cov2var(cov))
}


## calculations based on raw data
dat2est <- function(dat, regu=c(0,0,0)){
  n <- nrow(dat)
  (colMeans(dat)*n + regu[2])/(n+regu[1]) 
}

data2est <- function(data, regu=c(0,0,0)){
  lapply(data, function(dat) dat2est(dat, regu=regu))
}

dat2var <- function(dat){
  apply(dat, 2, var)/nrow(dat)
}

dat2se <- function(dat){
  sqrt(dat2var(dat))
}

dat2cov <- function(dat){
  cov(dat)/nrow(dat)
}


## calculations based on stats
stats2est <- function(stats){
  lapply(stats, function(s) s$est)
}

stats2tstat <- function(stats, mu0, alternative="greater"){
  stopifnot(is.list(stats))
  G <- length(stats)
  m <- length(stats[[1]]$est)
  stopifnot(length(mu0) == G)
  if(is.numeric(mu0)){
    mu0 <- lapply(1:G, function(g) rep(mu0[g], m))
  }
  stopifnot(is.list(mu0))
  
  lapply(1:G, function(g){
    tstat_tr((stats[[g]]$est - mu0[[g]])/stats[[g]]$se, alternative)
  })
}

tstat_tr <- function(x, alternative="greater"){
  switch(alternative,
         greater = x,
         two.sided = abs(x),
         less = -x)
}

## TODO: needed at all?
b2bb <- function(b, G=1:2){
  lapply(G, function(x) as.numeric(b==x))
}

tstat_cpe <- function(est, tstat){
  # tstat <- lapply(tstat, function(x){
  #   y <- x; y[!is.finite(x)] <- 1e+6; y
  # })
  b <- pargmin(args=est, rdm=TRUE) 
  #bb <- do.call(cbind, b2bb(b, 1:length(est)))
  do.call(cbind, tstat)[cbind(1:length(b), b)]
  #apply(mapply('*', bb, tstat), 1, sum)
}


