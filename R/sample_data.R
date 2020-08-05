#' sample_data
#'
#' @param n 
#' @param m 
#' @param prev 
#' @param random 
#' @param method 
#' @param pars 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
sample_data <- function(n=200,
                        m=10,
                        prev= c(0.5, 0.5),
                        random = FALSE,
                        method="roc",
                        pars = list(),
                        ...){
  pars <- c(pars, list(...))
  ## sample subgroup sizes
  ng <- sample_ng(n=n, prev=prev, random=random)
  ## sample binary data
  do.call(paste0("sample_data_", method), pars)
}


sample_ng <- function(n, prev, random=FALSE){
  stopifnot(all(prev >= 0))
  prev <- prev/sum(prev)
  if(random){
    ng <- as.numeric(extraDistr::rmnom(n=1, size=n, prob=prev))
  }else{
    ng <- round(n*prev, 0)
    ng[1] <- n-sum(ng[-1])
  }
  return(ng)
}
#sample_ng(200, c(0.5, 0.5), TRUE)


sample_data_mvn <- function(theta){
  return(NULL)
}

sample_data_roc <- function(n = c(100, 300),
                            m = 10,
                            r = 3,
                            auc = seq(0.85, 0.95, length.out = r),
                            rho = c(0.25, 0.25),
                            delta = 0,
                            e = 1,
                            k = 100,
                            corrplot = FALSE
                            ){
  ## argument checks:
  stopifnot(length(n) == 2)
  stopifnot(length(auc) == r)
  
  G <- length(n)
  ## mean vector (diseased: 1, healthy: 2):
  mu1 <- sqrt(2)* qnorm(auc) 
  mu2 <- rep(0, r) 
  ## covariance Matrix:
  C1 <- matrix(rho[1], r, r) + diag(rep(1-rho[1], r))
  C2 <- matrix(rho[2], r, r) + diag(rep(1-rho[2], r))
  ## subgroup samples:
  S1 <- mvtnorm::rmvnorm(n[1], mu1, C1) 
  S2 <- mvtnorm::rmvnorm(n[2], mu2, C2)
  ## cutoffs:
  q <- quantile(rbind(S1, S2), c(0.05, 0.95))
  cu <- seq(q[1], q[2], length.out = k)
  
  ## derive binary data:
  Y1 <- apply(S1, 2, function(x) lapply(cu, function(t) as.numeric(x>t)))
  Y1 <- do.call(cbind, lapply(Y1, function(x) do.call(cbind, x))) 
  
  Y2 <- apply(S2, 2, function(x) lapply(cu, function(t) as.numeric(x<=t)))
  Y2 <- do.call(cbind, lapply(Y2, function(x) do.call(cbind, x))) 
  
  ## True parameters
  Se <- colMeans(Y1)
  Sp <- colMeans(Y2)
  tau <- pmin(Se, Sp + delta)
  
  ## selection of 
  s <- sort(sample(1:(r*k), m, prob=tau^e))
  
  
  ## true parameters values
  info <- data.frame(
    cutoff = rep(cu, r)[s],
    theta1 = sapply(s, function(x) 1-pnorm(rep(cu, r)[x], mean=rep(mu1, each=k)[x])),
    theta2 = sapply(s, function(x) pnorm(rep(cu, r)[x], mean=rep(mu2, each=k)[x]))
  )
  info <- dplyr::mutate(info, tau=pmin(theta1, theta2 + delta) )
  
  
  #########################
  ## inspect parameter values:
  if(corrplot){

  R1 <- cov2cor(cov(Y1[,s]))
  R2 <- cov2cor(cov(Y2[,s]))
  corrplot::corrplot(R1)
  corrplot::corrplot(R2)
  
  return(info)
  }
  #########################
  
  out <- list(Y1[, s], Y2[, s])
  attr(out, "info") <- info
  ## return
  return(out)
}

#sample_data_roc()
#sample_data(pars=list(m=5))
