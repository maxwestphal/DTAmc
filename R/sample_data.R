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
sample_data <- function(n = 200,
                        m = 10,
                        prev = c(0.5, 0.5),
                        random = FALSE,
                        method = "roc",
                        pars = list(),
                        ...) {
  # if(is.character(prev)){
  #   prev <- as.numeric(strsplit(prev, "_"))
  # }
  ## sample subgroup sizes
  ng <- sample_ng(n = n, prev = prev, random = random)
  ## sample binary data
  args <- c(list(ng=ng, m=m), pars, list(...))
  do.call(paste0("sample_data_", method), args)
}

sample_ng <- function(n, prev, random = FALSE) {
  stopifnot(all(prev >= 0))
  prev <- prev / sum(prev)
  if (random) {
    ng <- as.numeric(extraDistr::rmnom(n = 1, size = n, prob = prev))
  } else{
    ng <- round(n * prev, 0)
    ng[1] <- n - sum(ng[-1])
  }
  return(ng)
}

sample_data_lfc <- function(ng = c(100, 300),
                            m = 10,
                            theta) {
  # TODO: use SEPM.LFC function
  return(NULL)
}

sample_data_roc <- function(ng = c(100, 300),
                            m = 10,
                            auc = seq(0.85, 0.95, length.out = 5),
                            rho = c(0.25, 0.25),
                            delta = 0,
                            e = 10,
                            k = 100,
                            corrplot = FALSE,
                            ...) {
  ## argument checks:
  stopifnot(length(ng) == 2)
  #stopifnot(length(auc) == r)
  
  r <- length(auc)
  G <- length(ng)
  ## mean vector (diseased: 1, healthy: 2):
  mu1 <- sqrt(2) * qnorm(auc)
  mu2 <- rep(0, r)
  ## covariance Matrix:
  C1 <- matrix(rho[1], r, r) + diag(1-rho[1], nrow=r, ncol=r)
  C2 <- matrix(rho[2], r, r) + diag(1-rho[2], nrow=r, ncol=r)
  ## subgroup samples:
  S1 <- mvtnorm::rmvnorm(ng[1], mu1, C1)
  S2 <- mvtnorm::rmvnorm(ng[2], mu2, C2)
  ## cutoffs:
  q <- quantile(rbind(S1, S2), c(0.05, 0.95))
  cu <- seq(q[1], q[2], length.out = k)
  
  ## derive binary data:
  Y1 <-
    apply(S1, 2, function(x)
      lapply(cu, function(t)
        as.numeric(x > t)))
  Y1 <- do.call(cbind, lapply(Y1, function(x)
    do.call(cbind, x)))
  
  Y2 <-
    apply(S2, 2, function(x)
      lapply(cu, function(t)
        as.numeric(x <= t)))
  Y2 <- do.call(cbind, lapply(Y2, function(x)
    do.call(cbind, x)))
  
  ## True parameters
  Se <- colMeans(Y1)
  Sp <- colMeans(Y2)
  tau <- pmin(Se, Sp + delta)
  
  ## selection of
  s <- sort(sample(1:(r * k), m, prob = tau ^ e))
  
  ## true parameters values
  info <- data.frame(
    auc = rep(auc, each = k)[s],
    cutoff = rep(cu, r)[s],
    theta1 = 1 - pnorm(rep(cu, r)[s], mean = rep(mu1, each = k)[s]),
    theta2 = pnorm(rep(cu, r)[s], mean = rep(mu2, each = k)[s])
  )
  info$tau = with(info, pmin(theta1, theta2 + delta))
  
  #########################
  ## inspect parameter values:
  if (corrplot) {
    R1 <- cov2cor(cov(Y1[, s]))
    R2 <- cov2cor(cov(Y2[, s]))
    corrplot::corrplot(R1)
    corrplot::corrplot(R2)
    
    return(info)
  }
  #########################
  
  out <- list(Y1[, s, drop=FALSE], Y2[, s, drop=FALSE])
  attr(out, "info") <- info
  ## return
  return(out)
}

#sample_data_roc()
#sample_data(pars=list(m=5))
