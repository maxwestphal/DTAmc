cv_bootstrap <- function(alpha, alternative, bst){
  if(any(!is.finite(bst))){
    message(paste0("Attention: ", sum(!is.finite(bst)), " bootstrap sample statistics (" ,
                   100*mean(!is.finite(bst)) , "%) are not finite!"))
  } # TODO
  
  cv <- quantile(bst, 1-alpha, na.rm=TRUE)
  
  cv <- switch(alternative,
               greater = c(-cv, Inf), 
               two.sided = c(-cv, cv),
               less = c(-Inf, cv)
  )
  
  return(unname(cv))
}

pval_bootstrap <- function(bst){
  function(tstat, alternative){
    sapply(tstat, function(x) mean(bst > x, na.rm=TRUE))  
  }
}

alpha_bootstrap <- function(alpha, alternative, bst){ # TODO
  NA
}


## create bootstrap sample of max test statistic
bootstrap_sample <- function(data, contrast, regu, alternative, pars){
  pars$type <- ifelse(is.null(pars$type), "pairs", pars$type)
  pars$nboot <- ifelse(is.null(pars$nboot), 2000, pars$nboot)
  pars$variant <- ifelse(is.null(pars$variant ), "lfc", pars$variant ) # TODO: ?!?
  
  stopifnot(pars$type %in% c("pairs", "wild"))
  stopifnot(pars$nboot %% 1 == 0)
  
  message(paste0("DTAmc: Drawing ", pars$nboot, " '", pars$type, "' bootstrap samples..."))
  
  do.call(paste0("bootstrap_sample_", pars$type), 
          args=list(data, contrast, regu, alternative, pars))
}


## pairs bootstrap
bs_draw_pairs <- function(data, G=length(data), ng=sapply(data, nrow)){
  lapply(1:G, function(g) data[[g]][sample(ng[g], replace=TRUE), ] )
}

bootstrap_sample_pairs <- function(data, contrast, regu=c(0,0,0), 
                                   alternative="greater", pars=list(nboot=2000)){
  G <- length(data); ng=sapply(data, nrow)
  mu0 <- stats2est(data2stats(data, contrast, regu))
  # pars=list(type="pairs", nboot=5000)
  # TODO: remove variant completely?!?
  # if(pars$variant == "max"){
  #   y <- max(do.call(pmin, args=mu0))
  #   mu0 <- lapply(1:G, function(g) rep(y, m))
  # }
  sapply(1:pars$nboot, function(b){
    st <- data2stats(bs_draw_pairs(data, G=G, ng=ng), contrast, regu);
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

bootstrap_sample_wild <- function(data, contrast, regu=c(0,0,0), 
                                  alternative="greater", pars=list(nboot=2000)){
  pars$dist <- ifelse(is.null(pars$dist), "Normal", pars$dist)
  pars$res_tra <- ifelse(is.null(pars$res_tra), 0, pars$res_tra)
  
  ## TODO: wild BS should only work with regu=c(2,1,?) or c(0,0,?)
  ## insert pseudo obs
  if(!all(regu==0)){
    data <- lapply(data, function(d){
      rbind(d, 0, 1)
    })
  }
  G <- length(data); ng=sapply(data, nrow); m <- ncol(data[[1]])
  
  mu0_raw <- stats2est(data2stats(data, contrast=define_contrast("raw"), regu))
  mu0 <- lapply(mu0_raw, function(x) as.numeric(contrast(data) %*% x))
  ## TODO: needed?
  # if(pars$variant == "max"){
  #   y <- max(do.call(pmin, args=mu0))
  #   mu0 <- lapply(1:G, function(g) rep(y, m))
  # }
  
  M <- lapply(1:G, function(g){
    matrix(mu0_raw[[g]], nrow=ng[g], ncol=length(mu0_raw[[g]]), byrow=TRUE)
  })
  D <- lapply(1:G, function(g){
    data[[g]] - M[[g]] 
  })
  sapply(1:pars$nboot, function(b){
    st <- data2stats(bs_draw_wild(M, D, G, ng, m,
                                  dist=pars$dist, res_tra=pars$res_tra),
                     contrast, regu, raw=TRUE);
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
  if(res_tra == 0){
    return(x)
  }
  if(res_tra == 1){
    return(sqrt(nrow(x)/(nrow(x) - ncol(x))) * x)
  }
  if(res_tra == 2){
    return(matrix( (1-h)^(-0.5), nrow=nrow(x), ncol=ncol(x), byrow=FALSE) * x)
  }
  if(res_tra == 3){
    return(matrix( (1-h)^(-1), nrow=nrow(x), ncol=ncol(x), byrow=FALSE) * x)
  }
  stop("pars$res_tra needs to be 0,1,2 or 3.")
} # TODO: 1 and 2 equivalent?

## LITERATURE:
# https://www.math.kth.se/matstat/gru/sf2930/papers/wild.bootstrap.pdf
# https://halshs.archives-ouvertes.fr/halshs-00175910/document
# http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1127.pdf




