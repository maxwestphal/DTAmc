#' dta
#'
#' @param data
#' @param method
#' @param pars
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
dta <- function(data = sample_data(),
                test = c("univariate", "bonferroni", "maxt", "bootstrap"),
                alpha = 0.05,
                alternative = c("two.sided", "less", "greater"),
                regu = FALSE,
                pars = list(),
                ...) {
  if (is.data.frame(data) | is.matrix(data)) {
    # TODO: transform to list form
    stop("transform to list format not yet implemented")
  }
  cov_needed <- c("maxt")
  pars <-
    c(
      list(
        data = data,
        stats = data2stats(
          data,
          regu = regu,
          covariance = match.arg(test) %in% cov_needed
        ),
        alpha = alpha,
        alternative = match.arg(alternative)
      ),
      list(...),
      pars
    )
  do.call(paste0("dta_", match.arg(test)), pars)
}

dta_univariate <- function(stats,
                           alpha = 0.05,
                           alternative = c("two.sided", "less", "greater"),
                           ...) {
  cv <- cv_uni(alpha, alternative)
  return(stats2results(stats, cv))
}

dta_bonferroni <- function(stats,
                           alpha = 0.05,
                           alternative = c("two.sided", "less", "greater"),
                           ...) {
  cv <- cv_uni(alpha / length(stats$est[[1]]), alternative)
  return(stats2results(stats, cv))
}

dta_maxt <- function(stats,
                     alpha = 0.05,
                     alternative = c("two.sided", "less", "greater"),
                     useSEPM=F,
                     ...) {
  
  if(useSEPM){
    #TODO: remove old version and ARGUMENT
    require(SEPM)
    define_hypothesis("accuracy.cp", threshold = c(0.5, 0.5)) %>%
      compare(comparison = data) %>%
      estimate(method = "beta.approx") %>%
      infer(method = "maxT") %>%
      summary() %>%
      return()
  }

  cv <- cv_maxt(stats, alpha, alternative)
  return(stats2results(stats, cv))
  
  
}

dta_boostrap <- function(stats,
                         alpha = 0.05,
                         alternative = c("two.sided", "less", "greater"),
                         ...) {
  ## prob: data needed as arg
}

cv_uni <- function(alpha = 0.05,
                   alternative = "two.sided") {
  cv <- numeric(2)
  cv[1] <- switch(
    alternative,
    two.sided = qnorm(alpha / 2),
    less = -Inf,
    greater = qnorm(alpha)
  )
  cv[2] <- switch(
    alternative,
    two.sided = qnorm(1 - alpha / 2),
    less = qnorm(1 - alpha),
    greater = Inf
  )
  return(cv)
}

#' @importFrom mvtnorm qmvnorm
cv_maxt <- function(stats, alpha = 0.05, alternative = "two.sided") {
  G <- length(stats$est)
  ## argmin vector am (where is minimum entry of estimates?)
  E <- do.call(cbind, stats$est)
  am <- apply(E, 1, argmin)
  
  ## binary coefficient matrices for all g (list)
  Bg <- lapply(1:G, function(g)
    diag(as.integer(am == g)))
  
  ## correlation matrices for all g (list)
  Rg <- lapply(stats$sigma, cov2cor)
  
  ## product BRG for all g (list)
  BRBg <-
    mapply(function(B, R) {
      B %*% R %*% B
    }, Bg, Rg, SIMPLIFY = FALSE)
  
  ## final overall correlation matrix of relevant estimates
  R <- do.call("+", BRBg)
  
  p <- switch (
    alternative,
    two.sided = 1 - alpha / 2,
    less = 1 - alpha,
    greater = 1 - alpha,
  )
  q <- mvtnorm::qmvnorm(p = p, corr = R)$quantile
  
  cv <- c(ifelse(alternative == "less",-Inf,-q),
          ifelse(alternative == "greater", Inf, q))
  
  return(cv)
}

stats2results <- function(stats, cv = 0) {
  G <- length(stats$est)
  do.call(cbind,
          lapply(1:G, function(g) {
            stat2result(
              est = stats$est[[g]],
              se = stats$se[[g]],
              cv = cv,
              g = g
            )
          }))
}

stat2result <- function(est, se, cv, g = "") {
  result <-
    cbind(
      estimate = est,
      lower = est + cv[1] * se,
      upper = est + cv[2] * se
    )
  colnames(result) <- paste0(colnames(result), g)
  return(result)
}
