stats2results <- function(stats, cv=c(-1.96, 1.96), pval_fun=pval_uni, 
                          benchmark=rep(0.5, length(stats)), 
                          alternative="greater", transformation="none", ...) {
  lapply(1:length(stats), function(g) stat2results(stats[[g]], cv, pval_fun, benchmark[g],
                                                   alternative, transformation))
}

#' @importFrom dplyr mutate_if
stat2results <- function(stat, cv=c(-1.96, 1.96), pval_fun=pval_uni,
                         benchmark=0.5, 
                         alternative="greater", transformation="none") {
  tf <- get_tf(transformation)
  est.t <- tf$est_link(stat$est)
  bm.t <- tf$est_link(benchmark)
  se.t <- tf$se_link(stat$se, stat$n, stat$est)
  tstat <- (est.t - bm.t)/se.t

  result <-
    data.frame(
      parameter = stat$names,
      alternative = altstr(alternative, benchmark),
      estimate = stat$est,
      lower = tf$inv(est.t + cv[1] * se.t), 
      upper = tf$inv(est.t + cv[2] * se.t),
      pvalue = pval_fun(tstat, alternative) 
    ) %>% dplyr::mutate_if(is.numeric, round, 4)
  rownames(result) <- NULL
  return(result)
}

altstr <- function(alternative, benchmark){
  paste0(switch(alternative,
                greater = " >= ",
                two.sided = " = ",
                less = " <= "),
         benchmark)
}