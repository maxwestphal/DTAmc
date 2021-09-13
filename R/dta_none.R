dta_none <- function(data = sample_data(seed=1337),
                     contrast = define_contrast("raw"),
                     benchmark = 0.5,
                     alpha = 0.05,
                     alternative = "greater",
                     transformation = "none",
                     regu = FALSE,
                     pars = list()
){  
  stats <- data2stats(data, contrast=contrast, regu = regu)
  cv <- cv_uni(alpha, alternative)
  
  ## output
  stats %>% 
    stats2results(
      cv=cv, pval_fun=pval_uni, benchmark,  
      alternative, transformation
    ) %>% 
    setattr(
      n = sapply(stats, function(x) x$n), m=ncol(data[[1]]), 
      alpha=alpha, alpha_adj=alpha, cv=cv
    ) %>% 
    return()
}



# Helper functions ----------------------------------------------------------------------------

cv_uni <- function(alpha = 0.05,
                   alternative = "two.sided",
                   ...) {
  c(switch(
    alternative,
    two.sided = qnorm(alpha / 2),
    less = -Inf,
    greater = qnorm(alpha)
  ),
  switch(
    alternative,
    two.sided = qnorm(1 - alpha / 2),
    less = qnorm(1 - alpha),
    greater = Inf
  )
  )
}

pval_uni <- function(tstat, alternative = "two.sided"){
  switch(
    alternative,
    two.sided = pnorm(abs(tstat), lower.tail = FALSE), 
    less = pnorm(-tstat, lower.tail = FALSE),
    greater = pnorm(tstat, lower.tail = FALSE)
  )
}

# TODO: pval_uni - check correctness