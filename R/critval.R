cv_uni <- function(alpha = 0.05,
                   alternative = "two.sided",
                   ...) {
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
cv_maxt <-
  function(stats,
           alpha = 0.05,
           alternative = "two.sided",
           ...) {
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
    
    cv <- c(ifelse(alternative == "less", -Inf, -q),
            ifelse(alternative == "greater", Inf, q))
    
    return(cv)
  }