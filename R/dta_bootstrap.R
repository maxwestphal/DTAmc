#' @importFrom boot boot
study_dta_bootstrap <- function(data = sample_data(seed=1337),
                          contrast = define_contrast("raw"),
                          benchmark = 0.5,
                          alpha = 0.05,
                          alternative = "greater",
                          transformation = "none",
                          regu = FALSE,
                          pars = list()
){

  stats <- data2stats(data, contrast=contrast, regu = regu)
  bst <- bootstrap_sample(data, contrast, regu, alternative, pars) 
  cv <- cv_bootstrap(alpha, alternative, bst)
  
  ## output
  out <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_bootstrap(bst), benchmark,  
    alternative, transformation
  )
  
  out %>%
    setattr(
      n = sapply(stats, function(x) x$n), m=nrow(out[[1]]), 
      alpha=alpha, alpha_adj=alpha_bootstrap(alpha, alternative, bst), cv=cv
    ) %>% 
    return()
}
