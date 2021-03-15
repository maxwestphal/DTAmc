stats2results <- function(stats, cv=c(-1.96, 1.96), pval_fun=pval_uni,
                          comparator=NULL, benchmark=rep(0.5, length(stats)), 
                          alternative="greater", transformation="none", ...) {
  lapply(1:length(stats), function(g) stat2results(stats[[g]], cv, pval_fun,
                                                   comparator, benchmark[g],
                                                   alternative, transformation))
}

stat2results <- function(stat, cv=c(-1.96, 1.96), pval_fun=pval_uni,
                         comparator=NULL, benchmark=0.5, 
                         alternative="greater", transformation="none") {
  tf <- get_tf(transformation)
  est.t <- tf$est_link(stat$est)
  bm.t <- tf$est_link(benchmark)
  se.t <- tf$se_link(stat$se, stat$n, stat$est)
  tstat <- (est.t - bm.t)/se.t

  result <-
    data.frame(
      hypothesis = hypotheses(names(stat$est), comparator, benchmark, alternative),
      estimate = stat$est,
      lower = tf$inv(est.t + cv[1] * se.t), 
      upper = tf$inv(est.t + cv[2] * se.t),
      pvalue = pval_fun(tstat, alternative) 
    )
  rownames(result) <- NULL
  return(result)
}

# TODO: needed after all?
hypotheses <- function(modnames, comparator, benchmark, alternative){
  paste0(modnames, 
         ifelse(is.null(comparator), "",
                paste0(" - ", modnames[comparator])),
         switch(alternative, 
                greater = " <= ",
                two.sided = " = ",
                less = " >= "),
         benchmark)
}