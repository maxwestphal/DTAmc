#' Analysis for Diagnostic Test Accuracy studies with multiple index tests.
#' 
#' Assess classification accuracy of multiple classifcation rules stratified 
#' by subpopulation, in the simplest case diseased and healthy individuals.
#'
#' @param data nxm binary matrix or data.frame (n observations of m binary decisions). 
#' Data should have values 0 (incorrect prediction) or 1 (correct prediction). 
#' \link{compare} provides a simple way to match predictions against true labels.
#' Alternatively, data can be a list of such arrays (all with m columns) defining different subsamples.
#' @param contrast specify contrast via \link{define_contrast}
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
#'
#' @return DTAmcResults object, which is a list of analysis results
#' @details 
#' Adjustment methods:
#' - "none" (default): no adjustment for multiplicity
#' - "bonferroni": Bonferroni adjustment
#' - "maxt": maxT adjustment
#' - "bootstrap": Bootstrap approach, whereby several details can be set via the 
#' 'pars' argument (see below)
#' - "mbeta": Bayesian approach
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
                contrast = define_contrast("raw"),
                benchmark = 0.75, 
                alpha = 0.025,
                alternative = c("greater", "two.sided", "less"), 
                adjustment = c("none", "bonferroni", "maxt", "bootstrap", "mbeta"),
                transformation = c("none", "logit"),
                regu = FALSE,
                pars = list(),
                ...) {
  
  ## check 'data' argument:
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
  
  ## check 'contrast' argument:
  stopifnot("contrastFun" %in% class(contrast))
  
  ## check 'benchmark' argument:
  if(! (length(benchmark) %in% c(1, length(data)))){
    stop("benchmark argument needs to be numeric of length one or of same length as data argument!")
  }
  if(length(benchmark) == 1){
    benchmark <- rep(benchmark, length(data))
  }
  stopifnot(all(abs(benchmark) < 1))
  
  ## check 'alpha' argument:
  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(alpha > 0 & alpha < 1)
  
  ## check 'pars' argument
  stopifnot(is.list(pars))
  
  ## prepare arguments for specific dta_xyz function:
  args <-
    list(
      data = data,
      contrast = contrast,
      benchmark = benchmark,
      alpha = alpha,
      alternative = match.arg(alternative),
      transformation = match.arg(transformation),
      regu = preproc_regu(regu),
      pars = c(list(...), pars)
    )
  
  ## calculate & label result:
  out <- do.call(paste0("dta_", match.arg(adjustment)), args)
  names(out) <- names(data)
  class(out) <- append(class(out), "DTAmcResults")
  
  return(out)
}
