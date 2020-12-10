b2bb <- function(b, G=1:2){
  lapply(G, function(x) as.numeric(b==x))
}

minStats <- function(data, s, group, regu, se_type=c("est", "data")){
  dd <- lapply(unique(group), function(g) data[s[group==g], , drop=FALSE])
  
  pp <- lapply(dd, function(dat) dat2est(dat, regu))
  b <- pargmin(args=pp)
  bb <- b2bb(b, sort(unique(group)))
  p <- apply(mapply('*', bb, pp), 1, sum)

  nn <- sapply(dd, nrow)
  n <- nn[b]
  
  se <- switch(match.arg(se_type), 
               est = sqrt((p*(1-p))/n),
               data = apply(mapply('*', bb, lapply(dd, dat2se)), 1, sum) )
  
  data.frame(p=p, se=se, n=b)
}

maxMinT <- function(data, s, group, regu, st0, 
                    se_type = c("est", "data"), res=0, rgen=rv, ...){
  if(any(res>0)){
    rb <- do.call(rgen, c(list(n=prod(dim(data))), list(...)))
    rb <- matrix(rb, nrow=nrow(data), ncol=ncol(data))
    data <- data + rb
  }
  ## TODO: REMOVE
  #sss <<- data 
  s <- sort(s)
  st <- minStats(data, s, group, regu, se_type = match.arg(se_type))
  tt <- (st$p-max(st0$p)) / st$se
  return(max(tt))
}


cv_bootstrap <- function(alpha, alternative, bst){
  if(any(!is.finite(bst))){
    message(paste0("Attention: ", sum(!is.finite(bst)), " bootstrap sample statistics (" ,
                   100*mean(!is.finite(bst)) , "%) are not finite!"))
  }
  
  cv <- switch(alternative,
               greater = c(quantile(bst, alpha, na.rm=TRUE), Inf),
               two.sided = quantile(bst, c(alpha/2, 1-alpha/2), na.rm=TRUE),
               less = c(-Inf, quantile(bst, 1-alpha, na.rm=TRUE))
               )
  
  return(unname(cv))
}

pval_bootstrap <- function(bst){
  
  function(tstat, alternative){
    switch(alternative,
           greater = sapply(tstat, function(x) mean(bst > x, na.rm=TRUE)), # TODO: correct?
           two.sided = sapply(tstat, function(x) mean(abs(bst) > abs(x), na.rm=TRUE)), ##TODO: FALSE, min(p1, p2)???
           less = sapply(tstat, function(x) mean(bst > x, na.rm=TRUE))# TODO: correct???
           
    )
  }
}

alpha_bootstrap <- function(alpha, alternative, bst){
  ## TODO
  NA
}

bootstrap_sample <- function(data, regu, pars){
  type <- ifelse(is.null(pars$type), "pairs", pars$type)
  nboot <- ifelse(is.null(pars$nboot), 2000, pars$nboot)
  stopifnot(type %in% c("pairs", "wild"))
  stopifnot(nboot %% 1 == 0)
  message(paste0("DTAmc: Drawing ", nboot, " '", type, "' bootstrap samples..."))
  group <- unlist(lapply(1:length(data), function(g) rep(g, nrow(data[[g]]))))
  db <- do.call(rbind, data)
  
  if(type == "pairs"){
    st0 <- minStats(db, 1:nrow(db), group, regu, se_type="est") 
    return(boot::boot(db, statistic = maxMinT,
                      R = nboot, strata=group, group=group,
                      regu=regu, st0 = st0, se_type="est")$t)
  }
  if(type == "wild"){
    dist <- ifelse(is.null(pars$dist), "Rademacher", pars$dist)
    ## TODO: colMeans --> stats with regu respected
    pp <- lapply(data, dat2est, regu=regu)
    n <- nrow(db); m <- ncol(db)
    pb <- do.call(rbind, lapply(group, function(g)pp[[g]]))
    rb <- db-pb
    #st0 <- maxMinT(db, 1:nrow(rb), group, regu, res=rb, rgen=rv)
    st0 <- minStats(db, 1:nrow(db), group, regu, se_type = "data") 
    
    ## TODO: rescaling of residuals??? http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1127.pdf
    #vb <- matrix(rv(n*m), nrow=n, ncol=m)
    #dbn <- db + rb*vb 
    #lapply(unique(group), function(g) dbn[group==g, , drop=FALSE]) %>% lapply(colMeans)

    return(boot::boot(pb, statistic = maxMinT,
                      R = nboot, strata=group, group=group,
                      regu=regu, st0 = st0, se_type="data",
                      res=rb, rgen=rv, dist=dist)$t)
    ## TODO: hand over dist arg 
  }
}

rv <- function(n=100, dist="Rademacher"){
  if(dist=="Rademacher"){
    return((rbinom(n, 1, 1/2)-0.5)*2)
  }
  if(dist=="Normal"){
    return(rnorm(n))
  }else{
    stop("argument 'dist' not recognized")
  }
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


