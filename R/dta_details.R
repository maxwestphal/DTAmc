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
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_uni, benchmark,  
    alternative, transformation
  )
  
  out$info <- list(n = sapply(stats, function(x) x$n), 
                   m = ncol(data[[1]]),
                   alpha = alpha,
                   alpha_adj = alpha,
                   cv = cv)
  
  return(out)
}

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
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_bonferroni, benchmark,  
    alternative, transformation
  )
  
  out$info <- list(n = sapply(stats, function(x) x$n), 
                   m = m,
                   alpha = alpha,
                   alpha_adj = alpha/ncol(data[[1]]),
                   cv = cv)
  return(out)
}

dta_maxt <- function(data = sample_data(seed=1337),
                     contrast = define_contrast("raw"),
                     benchmark = 0.5,
                     alpha = 0.05,
                     alternative = "greater",
                     transformation = "none",
                     regu = FALSE,
                     pars = list()
){
  
  stats <- data2stats(data, contrast=contrast, regu = regu)
  R <- cov2cor(active_cov(stats))
  cv <- cv_maxt(R, alpha, alternative)
  
  ## output
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_maxt(R), benchmark,  
    alternative, transformation
  )
  
  out$info <- list(n = sapply(stats, function(x) x$n), 
                   m = ncol(data[[1]]),
                   alpha = alpha,
                   alpha_adj = alpha_maxt(alpha, alternative, R),
                   cv = cv)
  return(out)
}

#' @importFrom boot boot
dta_bootstrap <- function(data = sample_data(seed=1337),
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
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_bootstrap(bst), benchmark,  
    alternative, transformation
  )
  
  out$info <- list(n = sapply(stats, function(x) x$n), 
                   m = ncol(data[[1]]),
                   alpha = alpha,
                   alpha_adj = alpha_bootstrap(alpha, alternative, bst),
                   cv = cv)
  
  return(out)
}

dta_mbeta <- function(data = sample_data(seed=1337),
                      contrast = define_contrast("raw"),
                      benchmark = 0.5,
                      alpha = 0.05,
                      alternative = "greater",
                      transformation = "none",
                      regu = FALSE,
                      pars = list()
) {
  
  ## output
  out <- list()
  out$results <- results_mbeta(data, contrast, benchmark, ## TODO: contrast needed here?
                               alpha, alternative, transformation, regu, pars)
  
  
  
  out$info <- list(n = sapply(data, nrow), 
                   m = ncol(data[[1]]),
                   alpha = alpha,
                   alpha_adj = alpha_bootstrap(alpha, alternative, bst), ## TODO: bst here makes not sense?
                   cv = NA)
  
  return(out)
}


## TODO:
## one sided versus two sided tests/intervals - interpretation?




## TODO: prob not needed?!?
# dta_wildbs <- function(data = sample_data(seed=1337),
#                        comparator = NULL,
#                        benchmark = 0.5,
#                        alpha = 0.05,
#                        alternative = "greater",
#                        transformation = "none",
#                        regu = FALSE,
#                        pars = list(),
#                        ...){
#   ## TODO
#   stop("adjustment wildbs not yet implemented")
# }

# dta_generic <- function(stats,
#                         alpha_adj,
#                         cv,
#                         data,
#                         comparator,
#                         alpha,
#                         alternative,
#                         transformation,
#                         regu,
#                         pars = list(),
#                         ...) {
#   ## output
#   out <- list()
#   out$results <- stats2results(
#     stats = stats,
#     cv = cv,
#     comp = comparator,
#     alt = alternative,
#     tf = get_tf(transformation)
#   )
#   
#   out$info <- list(n = stats$n, 
#                    m = length(stats$est[[1]]),
#                    alpha = alpha,
#                    alpha_adj = alpha,
#                    cv = cv)
#   
#   return(out)
# } 
# 
