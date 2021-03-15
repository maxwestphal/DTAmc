#' Analysis for Diagnostic Test Accuracy studies with multiple index tests.
#' 
#' Assess classification accuracy of multiple classifcation rules stratified 
#' by subpopulation, in the simplest case diseased and healthy individuals.
#'
#' @param data nxm binary matrix or data.frame (n observations of m binary decisions). 
#' Data should have values 0 (incorrect prediction) or 1 (correct prediction). 
#' \link{compare} provides a simple way to match predictions against true labels.
#' Alternatively, data can be a list of such arrays (all with m columns) defining different subsamples.
#' @param contrast (NOT YET IMPLEMENTED)
#' @param benchmark value to compare against (RHS), should have same length as data.
#' @param alpha numeric, significance level (default: 0.05)
#' @param alternative character, specify alternative hypothesis
#' @param adjustment character, specify type of statistical adjustment taken to address multiplicity
#' @param transformation character, define transformation to ensure results 
#' (e.g. point estimates, confidence limits) lie in unit interval ("none" (default) or "logit")
#' @param regu numeric vector of length 3, specify type of shrinkage.
#' Alternatively, logical of length one (TRUE := c(2, 1, 1/2), FALSE := c(0, 0, 0))
#' @param pars further parameters given as named list
#' @param ... additional named parameters
#' @param comparator LEGACY ARGUMENT, has no effect!
#'
#' @return DTAmcResults object, which is a list of analysis results
#' @details 
#' Adjustment methods:
#' - "none" (default): no adjustment for multiplicity
#' - "bonferroni": Bonferroni adjustment
#' - "maxt": maxT adjustment
#' ' - "mbeta": Bayesian approach
#' - "bootstrap": Bootstrap approach, whereby several details can be set via the 
#' 'pars' argument (see below)
#' 
#' Additional parameters to be adjusted via 'pars=list(par1=val1, par2=val2, ...)', e.g. list(type="pairs", nboot=10000), alternatively via '...'
#' - type: "pairs" (default) or "wild" = type of bootstrap
#' - nboot = number of bootstrap draws (default 5000)
#' - res_tra = 0,1,2 or 3 = type of residual transformation of wild boostrap (default = 0: no transformation)
#' (see https://www.math.kth.se/matstat/gru/sf2930/papers/wild.bootstrap.pdf)
#'
#' @export
#'
#' @examples
dta <- function(data = sample_data(seed=1337),
                contrast = NULL,
                benchmark = 0.5, 
                comparator = 1,
                alpha = 0.05,
                alternative = c("greater", "two.sided", "less"), 
                adjustment = c("none", "bonferroni", "maxt", "bootstrap", "mbeta"),
                transformation = c("none", "logit"),
                regu = FALSE,
                pars = list(),
                ...) {
  ## preprocess & check arguments:
  stopifnot(is.list(data))
  stopifnot(all(sapply(data, function(x) 
    any(class(x) %in% c("data.frame", "matrix")))))
  if(any(sapply(data, function(x) any(class(x) == "data.frame")))){
    data <- lapply(data, as.matrix)
  }
  stopifnot(all(diff(sapply(data, ncol))==0))
  if(!all(apply(sapply(data, colnames), 1, function(x) length(unique(x))==1 ))){
    stop("Expecting identical column names!")
  }
  if(length(benchmark) == 1){
    benchmark <- rep(benchmark, length(data))
  }
  stopifnot(all(abs(benchmark) < 1))
  stopifnot(is.numeric(alpha))
  stopifnot(alpha > 0 & alpha < 1)
  stopifnot(is.list(pars))
  
  ## prepare arguments for specific dta_XYZ function:
  args <-
    list(
      data = data,
      comparator = preproc_comp(comparator, data), # TODO; remove
      benchmark = benchmark,
      alpha = alpha,
      alternative = match.arg(alternative),
      transformation = match.arg(transformation),
      regu = preproc_regu(regu),
      pars = c(list(...), pars)
    )
  
  ## calculate & label result
  out <- do.call(paste0("dta_", match.arg(adjustment)), args)
  names(out$results) <- names(data)
  class(out) <- append(class(out), "DTAmcResults")
  
  return(out)
}

dta_none <- function(data = sample_data(seed=1337),
                     comparator = NULL,
                     benchmark = 0.5,
                     alpha = 0.05,
                     alternative = "greater",
                     transformation = "none",
                     regu = FALSE,
                     pars = list()
                     ){  
  
  stats <- data2stats(data, regu = regu)
  cv <- cv_uni(alpha, alternative)
  
  ## output
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_uni,
    comparator, benchmark,  
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
                           comparator = NULL,
                           benchmark = 0.5,
                           alpha = 0.05,
                           alternative = "greater",
                           transformation = "none",
                           regu = FALSE,
                           pars = list()
                           ){
  
  m <- ncol(data[[1]])
  
  stats <- data2stats(data, regu = regu)
  cv <- cv_bonferroni(m, alpha, alternative)
  
  ## output
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_bonferroni,
    comparator, benchmark,  
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
                     comparator = NULL,
                     benchmark = 0.5,
                     alpha = 0.05,
                     alternative = "greater",
                     transformation = "none",
                     regu = FALSE,
                     pars = list()
                     ){
  
  stats <- data2stats(data, regu = regu)
  R <- cov2cor(active_cov(stats))
  cv <- cv_maxt(R, alpha, alternative)
  
  ## output
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_maxt(R),
    comparator, benchmark,  
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
                          comparator = NULL,
                          benchmark = 0.5,
                          alpha = 0.05,
                          alternative = "greater",
                          transformation = "none",
                          regu = FALSE,
                          pars = list()
                          ){
  
  stats <- data2stats(data, regu = regu)
  bst <- bootstrap_sample(data, regu, pars, alternative) 
  cv <- cv_bootstrap(alpha, alternative, bst)
  
  ## output
  ## TODO: write function that preps output for all dta_XYZ, reduce redundancy
  out <- list()
  out$results <- stats2results(
    stats=stats, cv=cv, pval_fun=pval_bootstrap(bst),
    comparator, benchmark,  
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
                      comparator = NULL,
                      benchmark = 0.5,
                      alpha = 0.05,
                      alternative = "greater",
                      transformation = "none",
                      regu = FALSE,
                      pars = list()
                      ) {
  
  ## output
  out <- list()
  out$results <- results_mbeta(data, comparator, benchmark,
                               alpha, alternative, transformation, regu, pars)
  
  
  
  out$info <- list(n = sapply(data, nrow), 
                   m = ncol(data[[1]]),
                   alpha = alpha,
                   alpha_adj = alpha_bootstrap(alpha, alternative, bst),
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
