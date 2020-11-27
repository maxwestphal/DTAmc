cv_bonferroni <- function(m, alpha, alternative){
  cv_uni(alpha/m, alternative)
}

pval_bonferroni <- function(tstat, alternative){
  p <- pval_uni(tstat, alternative)*length(tstat)
  return(pmin(p, 1))
}