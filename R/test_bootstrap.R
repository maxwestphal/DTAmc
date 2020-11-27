b2bb <- function(b, G=1:2){
  lapply(G, function(x) as.numeric(b==x))
}

minStats <- function(data, s, group, regu){
  dd <- lapply(unique(group), function(g) data[s[group==g], , drop=FALSE])
  #dd <- lapply(unique(group), function(g) data[s, , drop=FALSE][group==g, , drop=FALSE])
  #pp <- lapply(dd, colMeans) ## TODO: replace with data2stats? SUBSTRACT cutoff?
  pp <- lapply(dd, function(dat) dat2est(dat, regu))
  b <- pargmin(args=pp)
  bb <- b2bb(b, sort(unique(group)))
  p <- apply(mapply('*', bb, pp), 1, sum)
  nn <- sapply(dd, nrow)
  data.frame(p=p, n=nn[b])
}

maxMinT <- function(data, s, group, regu){
  s <- sort(s)
  st <- minStats(data, s, group, regu)
  st0 <- minStats(data, 1:nrow(data), group, regu)
  tt <- (st$p-max(st0$p)) / sqrt(st$p*(1-st$p)/st$n)
  return(max(tt))
}


cv_bootstrap <- function(alpha, alternative, bst){
  
  cv <- switch(alternative,
               greater = c(quantile(bst, alpha), Inf),
               two.sided = quantile(bst, c(alpha/2, 1-alpha/2)),
               less = c(-Inf, quantile(bst, 1-alpha))
               )
  
  return(unname(cv))
}

pval_bootstrap <- function(bst){
  function(tstat, alternative){
    switch(alternative,
           greater = sapply(tstat, function(x) mean(bst > x)),
           two.sided = sapply(tstat, function(x) mean(abs(bst) > abs(x))),
           less = sapply(tstat, function(x) mean(bst > x))
           
    )
  }
}

alpha_bootstrap <- function(alpha, alternative, bst){
  ## TODO
  NA
}



### TODO: REMOVE BELOW

# mstar <- function(d, comp=0){
#   do.call(pmin, list(1:5, 5:1))
# }

#' teststat <- function(d, p0, regu){
#'   st <- DTAmc:::data2stats(list(d), regu=regu, covariance = FALSE)
#'   t <- (st$est[[1]]-p0)/sqrt(st$est[[1]]*(1-st$est[[1]])/st$n[1])
#'   return(t)
#' }
#' 
#' bootMaxMinT <- function(data, s, group, p0l, regu,
#'                         alternative = c("two.sided", "less", "greater")){
#'   ## TODO: p0l corecct?
#'   ## TODO: p0l name?
#'   a <- match.arg(alternative)
#'   f <- switch(a,
#'               two.sided = abs,
#'               less = function(x){-x},
#'               greater = function(x){x})
#'   s <- sort(s)
#'   tt <- lapply(unique(group), function(g) {
#'     f(teststat(d = data[s, ][group==g, ], p0=p0l[[g]], regu=regu))
#'   })
#'   ## TODO: correct min 
#'   #apply(do.call(rbind, p0l), 2, which.min)
#'   max(do.call(pmin, tt))
#' }
#' 
#' # coverage <- function(b, v, q, 
#' #                      alternative = c("two.sided", "less", "greater")){
#' #   a <- match.arg(alternative)
#' #   f <- switch(a,
#' #               two.sided = abs,
#' #               less = function(x){-x},
#' #               greater = function(x){x})
#' #   # TODO: extend to two sided
#' #   mean(apply(f(v) <= b, 1, all)) - q 
#' # }
#' 

#' cv_bootstrap <- function(stats,
#'                          alpha = 0.05,
#'                          alternative = "two.sided",
#'                          data = NA,
#'                          nboot = nboot,
#'                          regu = list(),
#'                          ...){
#'   stopifnot(!is.na(data))
#'   
#'   p0l <- stats$est # data2stats(data, regu=regu, covariance=FALSE)$est
#'   
#'   group <- do.call(c, lapply(1:length(data), function(g) rep(g, nrow(data[[g]]))))
#'   data <- do.call(rbind, data)
#'   
#'   bs <- boot::boot(data, statistic = bootMaxMinT,
#'                    R = nboot, strata=group, group=group, p0l=p0l, regu=regu,
#'                    alternative=alternative) 
#'   
#'   cv <- switch(alternative,
#'                two.sided = quantile(bs$t, 1-alpha/2),
#'                greater = quantile(bs$t, 1-alpha),
#'                less = quantile(bs$t, 1-alpha))
#'   cv <- switch(alternative,
#'                two.sided = c(-cv, cv),
#'                greater = c(-cv, Inf),
#'                less = c(-Inf, cv))
#'   return(cv)
#' }

#DTAmc:::cv_uni()
#DTAmc:::cv_maxt(DTAmc:::data2stats(dat), alpha=0.05, alternative="two.sided")


