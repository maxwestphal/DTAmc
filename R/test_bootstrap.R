cv_bootstrap <- function(alpha, alternative, bst){
  if(any(!is.finite(bst))){
    message(paste0("Attention: ", sum(!is.finite(bst)), " bootstrap sample statistics (" ,
                   100*mean(!is.finite(bst)) , "%) are not finite!"))
  } # TODO: remove this
  
  cv <- quantile(bst, 1-alpha, na.rm=TRUE)
  
  cv <- switch(alternative,
               greater = c(-cv, Inf), 
               two.sided = c(-cv, cv),
               less = c(-Inf, cv)
  )
  
  # cv <- switch(alternative,
  #              greater = c(quantile(bst, alpha, na.rm=TRUE), Inf), ### TODO: correct?
  #              two.sided = quantile(bst, c(alpha/2, 1-alpha/2), na.rm=TRUE),
  #              less = c(-Inf, quantile(bst, 1-alpha, na.rm=TRUE))
  #              )
  
  return(unname(cv))
}

pval_bootstrap <- function(bst){
  function(tstat, alternative){
    # switch(alternative,
    #        greater = sapply(tstat, function(x) mean(bst > x, na.rm=TRUE)), # TODO: correct?
    #        two.sided = sapply(tstat, function(x) mean(abs(bst) > abs(x), na.rm=TRUE)), ##TODO: FALSE, min(p1, p2)???
    #        less = sapply(tstat, function(x) mean(bst > x, na.rm=TRUE))# TODO: correct??
    #)
    sapply(tstat, function(x) mean(bst > x, na.rm=TRUE))  
  }
}

alpha_bootstrap <- function(alpha, alternative, bst){
  ## TODO
  NA
}


## create bootstrap sample of max test statistic
bootstrap_sample <- function(data, regu, pars, alternative="greater"){
  pars$type <- ifelse(is.null(pars$type), "pairs", pars$type)
  pars$nboot <- ifelse(is.null(pars$nboot), 5000, pars$nboot)
  
  stopifnot(pars$type %in% c("pairs", "wild"))
  stopifnot(pars$nboot %% 1 == 0)
  
  message(paste0("DTAmc: Drawing ", pars$nboot, " '", pars$type, "' bootstrap samples..."))
  
  do.call(paste0("bootstrap_sample_", pars$type), 
          args=list(data=data, regu=regu, pars=pars, alternative=alternative))
}


## pairs bootstrap
bs_draw_pairs <- function(data, G=length(data), ng=sapply(data, nrow)){
  lapply(1:G, function(g) data[[g]][sample(ng[g], replace=TRUE), ] )
}

bootstrap_sample_pairs <- function(data, regu=c(0,0,0), pars=list(nboot=5000),
                                   alternative="greater"){
  G <- length(data); ng=sapply(data, nrow)
  mu0 <- stats2est(data2stats(data, regu=regu))
  sapply(1:pars$nboot, function(b){
    st <- data2stats(bs_draw_pairs(data, G=G, ng=ng), regu=regu);
    tstat_cpe(stats2est(st), stats2tstat(st, mu0, alternative)) %>% max()
  })
}


## wild bootstrap
bs_draw_wild <- function(M, D, 
                         G=length(M), ng=sapply(M, nrow), m=ncol(M[[1]]),
                         dist = "Normal", res_tra=0){
  R <- rm(ng, m, dist); 
  lapply(1:G, function(g){
    M[[g]] + ( res_transform(D[[g]], res_tra = res_tra) *R[[g]] )
  })
}

bootstrap_sample_wild <- function(data, regu=c(0,0,0), pars=list(nboot=5000),
                                  alternative="greater"){
  pars$dist <- ifelse(is.null(pars$dist), "Normal", pars$dist)
  pars$res_tra <- ifelse(is.null(pars$res_tra), 0, pars$res_tra)
  
  ## insert pseudo obs
  ## TODO: wild BS should only work with regu=c(2,1,?) or c(0,0,?)
  if(!all(regu==0)){
    data <- lapply(data, function(d){
      rbind(d, 0, 1)
    })
  }
  G <- length(data); ng=sapply(data, nrow); m <- ncol(data[[1]])
  
  mu0 <- stats2est(data2stats(data, regu=regu))
  
  M <- lapply(1:G, function(g){
    matrix(mu0[[g]], nrow=ng[g], ncol=length(mu0[[g]]), byrow=TRUE)
  })
  D <- lapply(1:G, function(g){
    data[[g]] - M[[g]] 
  })
  sapply(1:pars$nboot, function(b){
    st <- data2stats(bs_draw_wild(M, D, G, ng, m,
                                  dist=pars$dist, res_tra=pars$res_tra),
                     regu=regu, raw=TRUE);
    tstat_cpe(stats2est(st), stats2tstat(st, mu0, alternative)) %>% max()
  })
}


## random vector for wild bootstrap
rv <- function(n=100, dist="Normal"){
  if(dist=="Rademacher"){
    return((rbinom(n, 1, 1/2)-0.5)*2)
  }
  if(dist=="Normal"){
    return(rnorm(n))
  }else{
    stop("argument 'dist' not recognized")
  }
}

## random matrix of correct dimensions
rm <- function(ng=c(5, 10), m=4, dist="Normal"){
  lapply(ng, function(n){
    matrix(rv(n, dist=dist), nrow=n, ncol=m, byrow=FALSE)
  })
}

res_transform <- function(x, h=rep(1/nrow(x), nrow(x)), res_tra=0){
  ## reference:
  # https://www.math.kth.se/matstat/gru/sf2930/papers/wild.bootstrap.pdf
  if(res_tra == 0){
    return(x)
  }
  if(res_tra == 1){
    sqrt(nrow(x)/(nrow(x) - ncol(x))) * x
  }
  if(res_tra == 2){
    matrix( (1-h)^(-0.5), nrow=nrow(x), ncol=ncol(x), byrow=FALSE) * x
  }
  if(res_tra == 3){
    matrix( (1-h)^(-1), nrow=nrow(x), ncol=ncol(x), byrow=FALSE) * x
  }
  stop("pars$res_tra needs to be 0,1,2 or 3.")
}

## LITERATURE:
# https://www.math.kth.se/matstat/gru/sf2930/papers/wild.bootstrap.pdf
# https://halshs.archives-ouvertes.fr/halshs-00175910/document
# http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1127.pdf






## TODO: remove old stuff below
# minStats <- function(data, s, group, regu, se_type=c("est", "data")){
#   dd <- lapply(unique(group), function(g) data[s[group==g], , drop=FALSE])
#   
#   pp <- lapply(dd, function(dat) dat2est(dat, regu))
#   b <- pargmin(args=pp)
#   bb <- b2bb(b, sort(unique(group)))
#   p <- apply(mapply('*', bb, pp), 1, sum)
#   
#   nn <- sapply(dd, nrow)
#   n <- nn[b]
#   
#   se <- switch(match.arg(se_type), 
#                est = sqrt((p*(1-p))/n),
#                data = apply(mapply('*', bb, lapply(dd, dat2se)), 1, sum) )
#   
#   data.frame(p=p, se=se, n=b)
# }
# 
# maxMinT <- function(data, s, group, regu, st0, 
#                     se_type = c("est", "data"), res=0, rgen=rv, ...){
#   if(any(res>0)){
#     rb <- do.call(rgen, c(list(n=prod(dim(data))), list(...)))
#     rb <- matrix(rb, nrow=nrow(data), ncol=ncol(data))
#     data <- data + rb
#   }
#   ## TODO: REMOVE
#   #sss <<- data 
#   s <- sort(s)
#   st <- minStats(data, s, group, regu, se_type = match.arg(se_type))
#   tt <- (st$p-max(st0$p)) / st$se
#   return(max(tt))
# }



# data <- sample_data_lfc()
# 
# 
# db <- bs_draw_pairs(data)
# mu0 <- stats2est(data2stats(data, regu=c(2,1,0.5)))
# 
# #lapply(data, colMeans)
# #lapply(db, colMeans)
# 
# #sapply(1:1000, function(i)lapply(bs_draw_pairs(data), colMeans)[[1]]) %>% apply(1, function(x) all(diff(x) == 0))
# 
# est <- data2stats(db, regu=c(2,1,0.5)) %>% stats2est()
# tstat <- data2stats(db, regu=c(2,1,0.5)) %>% stats2tstat(mu0)
# tstat_cpe(est, tstat) %>% max()






## TODO: REMOVE OLD VERSION BELOW
# bootstrap_sample <- function(data, regu, pars){
#   type <- ifelse(is.null(pars$type), "pairs", pars$type)
#   nboot <- ifelse(is.null(pars$nboot), 2000, pars$nboot)
#   stopifnot(type %in% c("pairs", "wild", "wild21")) # TODO: fix
#   stopifnot(nboot %% 1 == 0)
#   message(paste0("DTAmc: Drawing ", nboot, " '", type, "' bootstrap samples..."))
#   group <- unlist(lapply(1:length(data), function(g) rep(g, nrow(data[[g]]))))
#   db <- do.call(rbind, data)
#   
#   if(type == "pairs"){
#     st0 <- minStats(db, 1:nrow(db), group, regu, se_type="est") 
#     return(boot::boot(db, statistic = maxMinT,
#                       R = nboot, strata=group, group=group,
#                       regu=regu, st0 = st0, se_type="est")$t)
#   }
#   if(type == "wild"){
#     dist <- ifelse(is.null(pars$dist), "Rademacher", pars$dist)
#     ## TODO: colMeans --> stats with regu respected
#     pp <- lapply(data, dat2est, regu=regu)
#     n <- nrow(db); m <- ncol(db)
#     pb <- do.call(rbind, lapply(group, function(g)pp[[g]]))
#     rb <- db-pb
#     #st0 <- maxMinT(db, 1:nrow(rb), group, regu, res=rb, rgen=rv)
#     st0 <- minStats(db, 1:nrow(db), group, regu, se_type = "data") 
#     
#     ## TODO: rescaling of residuals??? http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1127.pdf
#     #vb <- matrix(rv(n*m), nrow=n, ncol=m)
#     #dbn <- db + rb*vb 
#     #lapply(unique(group), function(g) dbn[group==g, , drop=FALSE]) %>% lapply(colMeans)
# 
#     return(boot::boot(pb, statistic = maxMinT,
#                       R = nboot, strata=group, group=group,
#                       regu=regu, st0 = st0, se_type="data",
#                       res=rb, rgen=rv, dist=dist)$t)
#     ## TODO: hand over dist arg 
#   }
#   if(type == "wild21"){
#     # TODO: allow other than normal weights
#     return(wildbs(data, nboot=nboot))
#   }
# }



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


### OLD
# if(type == "wild"){
#   dist <- ifelse(is.null(pars$dist), "Rademacher", pars$dist)
#   ## TODO: colMeans --> stats with regu respected
#   pp <- lapply(data, dat2est, regu=regu)
#   n <- nrow(db); m <- ncol(db)
#   pb <- do.call(rbind, lapply(group, function(g)pp[[g]]))
#   rb <- db-pb
#   #st0 <- maxMinT(db, 1:nrow(rb), group, regu, res=rb, rgen=rv)
#   st0 <- minStats(db, 1:nrow(db), group, regu, se_type = "data")
#   
#   ## TODO: rescaling of residuals??? http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1127.pdf
#   #vb <- matrix(rv(n*m), nrow=n, ncol=m)
#   #dbn <- db + rb*vb
#   #lapply(unique(group), function(g) dbn[group==g, , drop=FALSE]) %>% lapply(colMeans)
#   
#   return(boot::boot(pb, statistic = maxMinT,
#                     R = nboot, strata=group, group=group,
#                     regu=regu, st0 = st0, se_type="data",
#                     res=rb, rgen=rv, dist=dist)$t)
#   ## TODO: hand over dist arg
# }

