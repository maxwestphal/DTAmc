#' Sample binary data
#'
#' @param n overall sample size
#' @param m number of models
#' @param prev prevalence
#' @param random random sampling (TRUE) or fixed group sample sizes
#' @param method 
#' @param seed
#' @param pars
#' @param ...
#' @return
#' @export
#'
#' @examples
sample_data <- function(n = 200,
                        prev = c(0.5, 0.5),
                        random = FALSE,
                        m = 10,
                        method = c("roc", "lfc", "prob"),
                        seed = NULL,
                        pars = list(),
                        ...) {
  method <- match.arg(method)
  
  ng <- sample_ng(n = n, prev = prev, random = random)
  
  args <- c(list(ng = ng, m = m), pars, list(...))
  do.call(paste0("sample_data_", method), args)
}

#' Sample binary data (one sample)
#'
#' @param n
#'
#' @param prob
#' @param R
#'
#' @importFrom bindata rmvbin
sample_data_prob <- function(n = 100,
                             prob = c(0.8, 0.8),
                             R = diag(length(prob))) {
  if (length(prob) == 1)
    return(bindata::rmvbin(n = n, margprob = prob))
  return(bindata::rmvbin(
    n = n,
    margprob = prob,
    bincorr = R
  ))
}

sample_ng <- function(n = 100,
                      prev = c(0.5, 0.5),
                      random = FALSE) {
  stopifnot(all(prev >= 0))
  prev <- prev / sum(prev)
  ng <- rep(0, length(prev))
  while (any(ng == 0)) {
    if (random) {
      ng <- as.numeric(extraDistr::rmnom(
        n = 1,
        size = n,
        prob = prev
      ))
    } else{
      ng <- round(n * prev, 0)
      ng[1] <- n - sum(ng[-1])
    }
  }
  return(ng)
}
