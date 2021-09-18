study_dta_maxt <- function(data = generate_data(seed=1337),
                     contrast = define_contrast("raw"),
                     benchmark = 0.5,
                     alpha = 0.05,
                     alternative = "greater",
                     transformation = "none",
                     regu = FALSE,
                     pars = list()
){
  
  stats <- data2stats(data, contrast=contrast, regu = regu)
  R <- stats::cov2cor(active_cov(stats))
  cv <- cv_maxt(R, alpha, alternative)
  
  ## output
  stats %>% 
    stats2results(
      cv=cv, pval_fun=pval_maxt(R), benchmark,  
      alternative, transformation
    ) %>% 
    setattr(
      n = sapply(stats, function(x) x$n), m=ncol(R), 
      alpha=alpha, alpha_adj=alpha_maxt(alpha, alternative, R), cv=cv
    ) %>% 
    return()
}


# Helper functions ----------------------------------------------------------------------------

active_cov <- function(stats){
  G <- length(stats)
  ## argmin vector am (where is minimum entry of estimates?)
  E <- sapply(1:G, function(g) stats[[g]]$est)
  am <- apply(E, 1, argmin)
  
  ## binary coefficient matrices for all g (list)
  Bg <- lapply(1:G, function(g)
    diag(as.integer(am == g)))
  
  ## covariance matrices for all g (list)
  Sg <- lapply(stats, function(x) x$cov)
  
  ## product BSB for all g (list)
  BSBg <-
    mapply(function(B, S) {
      B %*% S %*% B
    }, Bg, Sg, SIMPLIFY = FALSE)
  
  ## final overall covariance matrix of 'active' estimates
  do.call("+", BSBg)
}

#' @importFrom mvtnorm qmvnorm pmvnorm
cv_maxt <-
  function(R,
           alpha = 0.05,
           alternative = "greater",
           ...) {
    
    
    tail <- switch (
      alternative,
      two.sided = "both.tails",
      less = "upper.tail",
      greater = "lower.tail",
    )
    
    q <- mvtnorm::qmvnorm(p = 1-alpha, sigma = R, tail = tail)$quantile
    
    cv <- switch(alternative,
                 greater = c(-q, Inf),
                 two.sided = c(-q, q),
                 less = c(-Inf, -q)
    )
    
    return(cv)
  }

pval_maxt <- function(R){
  function(tstat, alternative){
    m <- length(tstat)
    switch(alternative,
           greater = sapply(tstat, function(x) 
             mvtnorm::pmvnorm(lower = rep(x, m), upper = rep(Inf, m), corr=R)[1]),
           two.sided = sapply(tstat, function(x) 
             mvtnorm::pmvnorm(lower = rep(-abs(x), m), upper = rep(abs(x), m), corr=R)[1]),
           less = sapply(tstat, function(x) 
             mvtnorm::pmvnorm(lower = rep(-Inf, m), upper = rep(x, m), corr=R)[1]))
  } 
}
# TODO: pval_maxt - correct pvals? NO, compare: study_dta(adj="maxt", regu=T)

alpha_maxt <- function(alpha, alternative, R){ 
  NA
}
# TODO: alpha_maxt
