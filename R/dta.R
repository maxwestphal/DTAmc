#' Diagnostic test accuracy analysis
#' 
#' Assess classification accuracy stratified by subpopulation, in the simplest case 
#' diseased and healthy individuals.
#'
#' @param data nxm binary matrix or data.frame, i.e. with values 0 and 1 (n observations of m binary decisions). 
#' Alternatively, data can be a list of such arrays (all with n columns) defining different subsamples.
#' @param group integer vector (default: NULL) of length n with G distinct values, usually 1:G.
#' Defines subsamples of data. Not needed if data is given as a list.
#' @param comparator numeric vector 
#' @param adjustment character, specify type of statistical adjustment taken to adress multiplicity
#' @param alpha numeric (default: 0.05), in unit interval 
#' @param alternative character, specify alternative hypothesis
#' @param transformation character, define transformation to ensure results 
#' (e.g. point estimates, confidence limits) lie in unit interval
#' @param regu numeric vector of length 3, specify type of shrinkage. TODO: DETAILS
#' Alternatively, logical of length one (TRUE := c(2, 1, 1/2), FALSE := c(0, 0, 0))
#' @param pars further parameters given as named list
#' @param ... 
#'
#' @return data.frame with analysis results
#' @export
#'
#' @examples
dta <- function(data = sample_data(),
                group = NULL,
                comparator = 0.5,
                adjustment = c("none", "bonferroni", "maxt", "bootstrap"),
                alpha = 0.05,
                alternative = c("two.sided", "less", "greater"),
                transformation = c("none", "logit"),
                regu = FALSE,
                pars = list(),
                ...) {
  
  if (is.data.frame(data) | is.matrix(data)) {
    stopifnot(!is.null(group) & length(group) == nrow(data))
    data <- split(data, group)
  }
  
  adjustment <- match.arg(adjustment)
  alternative <- match.arg(alternative)
  
  cv <- switch(
    adjustment,
    none = cv_uni(alpha, alternative),
    bonferroni = cv_uni(alpha / length(stats$est[[1]]), alternative),
    maxt = cv_maxt(stats, alpha, alternative),
    bootstrap = stop("bootstap not yet implemented")
  )
  
  cov_needed <- c("maxt")
  args <-
    c(
      list(
        data = data,
        stats = data2stats(
          data,
          regu = regu,
          covariance = match.arg(adjustment) %in% cov_needed
        ),
        cv = cv, 
        comp = comp_preproc(comparator, data),
        alt = match.arg(alternative),
        tf = get_tf(match.arg(transformation))
      ),
      list(...),
      pars
    )
  #do.call(paste0("dta_", match.arg(adjustment)), args)
  return(do.call(stats2results, args))
}

stats2results <- function(stats, cv, comp, alt, tf, ...) {
  do.call(cbind,
          lapply(1:length(stats$est), function(g) {
            stat2result(
              est = stats$est[[g]],
              se = stats$se[[g]],
              n = stats$n[g],
              cv = cv,
              comp = comp[g], # TODO:
              alt = alt,
              tf = tf,
              g = g
            )
          }))
}

stat2result <- function(est, se, n, cv, comp, alt, tf, g = "") {
  est.t <- tf$est_link(est)
  comp.t <- tf$est_link(comp)
  se.t <- tf$se_link(se, n, est)
  tstat <- (est.t - comp.t)/se.t
  result <-
    cbind(
      estimate = est,
      lower = tf$inv(est.t + cv[1] * se.t), 
      upper = tf$inv(est.t + cv[2] * se.t),
      comparator = comp,
      pvalue = ifelse(alt=="two.sided", 2, 1) * pnorm(tstat, lower.tail=FALSE) 
    )
  colnames(result) <- paste0(colnames(result), g)
  return(result)
}


# TODO: remove stuff
# est.t <- do.call(t$link, list(est))
# se.t <- do.call(t$se.t, list(se=se, n=e$n, estimate=est))
# 
# 
# 
# sigma.t <- diag(se.t, length(se.t)) %*% R %*% diag(se.t, length(se.t))
# 
# CI.t <- data.frame(estimate = est.t,
#                    lower = est.t - cv*se.t,
#                    upper = est.t + cv*se.t)
# CI <- as.data.frame(do.call(t$inv, list(CI.t)))



# dta_none <- function(stats,
#                      alpha = 0.05,
#                      alternative = c("two.sided", "less", "greater"),
#                      comp, alt, tf, ...) {
#   cv <- cv_uni(alpha, alternative)
#   return(stats2results(stats, cv, comp, alt, tf))
# }
# 
# dta_bonferroni <- function(stats,
#                            alpha = 0.05,
#                            alternative = c("two.sided", "less", "greater"),
#                            ...) {
#   cv <- cv_uni(alpha / length(stats$est[[1]]), alternative)
#   return(stats2results(stats, cv, comp, alt, tf))
# }
# 
# dta_maxt <- function(stats,
#                      alpha = 0.05,
#                      alternative = c("two.sided", "less", "greater"),
#                      useSEPM=F,
#                      ...) {
#   
#   if(useSEPM){
#     #TODO: remove old version and ARGUMENT
#     require(SEPM)
#     define_hypothesis("accuracy.cp", threshold = c(0.5, 0.5)) %>%
#       compare(comparison = data) %>%
#       estimate(method = "beta.approx") %>%
#       infer(method = "maxT") %>%
#       summary() %>%
#       return()
#   }
#   
#   cv <- cv_maxt(stats, alpha, alternative)
#   return(stats2results(stats, cv, comp, alt, tf))
#   
#   
# }
# 
# dta_boostrap <- function(stats,
#                          alpha = 0.05,
#                          alternative = c("two.sided", "less", "greater"),
#                          ...) {
#   ## prob: data needed as arg
#   stop("Bootstrap not yet implemented")
# }