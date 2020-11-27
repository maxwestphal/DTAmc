#' Sample binary data (ROC model)
#'
#' @param m 
#' @param auc 
#' @param rho 
#' @param delta 
#' @param e 
#' @param k 
#' @param corrplot 
#' @param ... 
#' @param n 
#' @param prev 
#' @param random 
#' @param modnames 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom corrplot corrplot
sample_data_roc <- function(n = 100,
                            prev = c(0.5, 0.5),
                            random = FALSE,
                            m = 10,
                            auc = seq(0.85, 0.95, length.out = 5),
                            rho = c(0.25, 0.25),
                            e = 10,
                            k = 100,
                            delta = 0,
                            modnames = paste0("model", 1:m),
                            seed = NULL,
                            corrplot = FALSE,
                            ...) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  ng <- sample_ng(n, prev, random)
  stopifnot(length(ng) == 2)

  
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
  
  
  ## selection of models
  theta1 <- 1 - pnorm(rep(cu, r), mean = rep(mu1, each = k))
  theta2 <- pnorm(rep(cu, r), mean = rep(mu2, each = k))
  tau <- pmin(theta1, theta2 + delta)
  s <- sort(sample(1:(r * k), m, prob = tau ^ e))
  
  #########################
  ## true parameters values
  info <- data.frame(
    model = modnames,
    auc = rep(auc, each = k)[s],
    cutoff = rep(cu, r)[s],
    se = theta1[s],
    sp = theta2[s]
  )
  
  Y1s <- Y1[, s, drop=FALSE]
  Y2s <- Y2[, s, drop=FALSE]
  
  #########################
  ## inspect parameter values:
  if (corrplot) {
    R1 <- cov2cor(cov(Y1s))
    R2 <- cov2cor(cov(Y2s))
    corrplot::corrplot(R1)
    corrplot::corrplot(R2)
    return(info)
  }
  
  colnames(Y1s) <- colnames(Y2s) <- modnames
  out <- list(Y1s, Y2s)
  
  names(out) <- c("sensitivity", "specificity")
  attr(out, "info") <- info

  return(out)
}

#sample_data_roc()
#sample_data(pars=list(m=5))
