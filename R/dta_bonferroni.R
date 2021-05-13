dta_bonferroni <- function(data = sample_data(seed=1337),
                           contrast = define_contrast("raw"),
                           benchmark = 0.5,
                           alpha = 0.05,
                           alternative = "greater",
                           transformation = "none",
                           regu = FALSE,
                           pars = list()
){
  
  K <- contrast(data)
  m <- nrow(K)
  
  stats <- data2stats(data, contrast=contrast, regu = regu)
  cv <- cv_bonferroni(m, alpha, alternative)
  
  ## output
  stats %>% 
    stats2results(
      cv=cv, pval_fun=pval_bonferroni, benchmark,  
      alternative, transformation
    ) %>% 
    setattr(
      n = sapply(stats, function(x) x$n), m=m, 
      alpha=alpha, alpha_adj=alpha/m, cv=cv
    ) %>% 
    return()
}

# Helper functions ----------------------------------------------------------------------------

cv_bonferroni <- function(m, alpha, alternative){
  cv_uni(alpha/m, alternative)
}

pval_bonferroni <- function(tstat, alternative){
  p <- pval_uni(tstat, alternative)*length(tstat)
  return(pmin(p, 1))
}