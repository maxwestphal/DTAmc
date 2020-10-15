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
    
    cv <- c(ifelse(alternative == "less",-Inf,-q),
            ifelse(alternative == "greater", Inf, q))
    
    return(cv)
  }

teststat <- function(d, p0, regu){
  st <- DTAmc:::data2stats(list(d), regu=regu, covariance = FALSE)
  t <- (st$est[[1]]-p0)/sqrt(st$est[[1]]*(1-st$est[[1]])/st$n[1])
  return(t)
}

bootMaxMinT <- function(data, s, group, p0l, regu,
                        alternative = c("two.sided", "less", "greater")){
  a <- match.arg(alternative)
  f <- switch(a,
              two.sided = abs,
              less = function(x){-x},
              greater = function(x){x})
  s <- sort(s)
  tt <- lapply(unique(group), function(g) {
    f(teststat(d = data[s, ][group==g, ], p0=p0l[[g]], regu=regu))
  })
  max(do.call(pmin, tt))
}

# coverage <- function(b, v, q, 
#                      alternative = c("two.sided", "less", "greater")){
#   a <- match.arg(alternative)
#   f <- switch(a,
#               two.sided = abs,
#               less = function(x){-x},
#               greater = function(x){x})
#   # TODO: extend to two sided
#   mean(apply(f(v) <= b, 1, all)) - q 
# }

#' @importFrom boot boot
cv_bootstrap <- function(stats,
                         alpha = 0.05,
                         alternative = "two.sided",
                         data = NA,
                         nboot = nboot,
                         regu = list(),
                         ...){
  stopifnot(!is.na(data))
  
  p0l <- stats$est # data2stats(data, regu=regu, covariance=FALSE)$est
  
  group <- do.call(c, lapply(1:length(data), function(g) rep(g, nrow(data[[g]]))))
  data <- do.call(rbind, data)
  
  bs <- boot::boot(data, statistic = bootMaxMinT,
                   R = nboot, strata=group, group=group, p0l=p0l, regu=regu,
                   alternative=alternative) 
  
  cv <- switch(alternative,
               two.sided = quantile(bs$t, 1-alpha/2),
               greater = quantile(bs$t, 1-alpha),
               less = quantile(bs$t, 1-alpha))
  cv <- switch(alternative,
               two.sided = c(-cv, cv),
               greater = c(-cv, Inf),
               less = c(-Inf, cv))
  return(cv)
}

#DTAmc:::cv_uni()
#DTAmc:::cv_maxt(DTAmc:::data2stats(dat), alpha=0.05, alternative="two.sided")
