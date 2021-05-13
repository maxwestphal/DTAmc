dta_mbeta <- function(data = sample_data(seed=1337),
                      contrast = define_contrast("raw"),
                      benchmark = 0.5,
                      alpha = 0.05,
                      alternative = "greater",
                      transformation = "none",
                      regu = FALSE,
                      pars = list()
) {
  
  stopifnot(transformation == "none")
  nrep <- ifelse(is.null(pars$nrep), 5000, pars$nrep)
  
  G <- length(data)
  m <- ncol(data[[1]])
  type <- attr(contrast(data), "type")
  
  moms <- lapply(data, function(d){
    add_moments(prior_moments(m, regu), data_moments(d))})
  stats <- data2stats(data, contrast, regu, raw=TRUE)
  
  ## posterior sample:
  pss <- sample_mbeta(moms) %>% transform_sample(tr=function(s){s %*% contrast(data)})
  
  ## credible region:
  qstar <- uniroot(f=eval_cr, interval=c(0, 0.5),
                   moms=moms, pss=pss, type=type,
                   alpha=alpha, alternative=alternative)$root
  crs <- get_cr(moms, pss, type=type, q=qstar, alternative=alternative) 
  
  ## output:
  lapply(1:G, function(g) {
    data.frame(
      parameter = stats[[g]]$names,
      alternative = altstr(alternative, benchmark),
      estimate = stats[[g]]$est,
      lower = crs[[g]]$lower,
      upper = crs[[g]]$upper,
      pvalue = NA
    )
  }) %>%
    setattr(
      n = sapply(data, nrow), m=m, 
      alpha=alpha, alpha_adj=qstar, cv=NA,
      class = c("list", "DTAmcResults")
    ) %>% 
    return()
}



# Helper functions ----------------------------------------------------------------------------

eval_cr <- function(q, moms, pss, type, alpha, alternative){
  crs <- get_cr(moms, pss, type, q, alternative)
  coverage(crs, pss) - (1-alpha)
}

coverage <- function(crs, pss){
  nrep <- nrow(pss[[1]]); m <- ncol(pss[[1]]); G <- length(pss)

  L <- lapply(crs, function(x) matrix(x$lower, nrow=nrep, ncol=m, byrow=TRUE))
  U <- lapply(crs, function(x) matrix(x$upper, nrow=nrep, ncol=m, byrow=TRUE))

  C <- lapply(1:G, function(g) covered(pss[[g]], L[[g]], U[[g]]))
  W <- matrix(sample(1:G, nrep*m, replace=TRUE), nrow=nrep, ncol=m, byrow=TRUE)

  mean(apply(Reduce("+", lapply(1:G, function(g) (W==g)*C[[g]])), 1, min) == 1)
}




covered <- function(S, L, U){
  S >= L & S <= U
}

sample_mbeta <- function(moms, nrep=5000, proj.pd=FALSE){
  lapply(moms, function(mom) sample_mbeta1(mom, nrep, proj.pd))
}

#' @importFrom  MCMCpack rdirichlet
#' @importFrom  copula normalCopula
#' @importFrom  copula P2p
#' @importFrom  copula mvdc
#' @importFrom  copula rMvdc
#' @importFrom  Matrix nearPD
sample_mbeta1 <- function(mom, nrep=5000, proj.pd = FALSE){
  m <- nrow(mom$moments)
  R <- cov2cor(mom2cov(mom))
  if(proj.pd){R <- as.matrix(Matrix::nearPD(R)$mat)}
  cop <- copula::normalCopula(param=copula::P2p(R), dim=m, dispstr = "un")
  mp <- margin_params(mom, c("shape1", "shape2"))
  
  ## output:
  copula::mvdc(cop, margins = rep("beta", m), paramMargins = mp) %>%
    copula::rMvdc(n=nrep) %>%
    return()
}

margin_params <- function(mom, n=c("alpha", "beta")){
  lapply(1:nrow(mom$moments), function(j){
    y <- list(mom$moments[j,j], mom$n - mom$moments[j,j])
    names(y) <- n
    return(y)
  })
}

transform_sample <- function(sample, tr=function(s){s %*% contrast(data)}){
  lapply(sample, tr)
}

get_cr1 <- function(mom, ps, type="raw", q=0.05, alternative="two.sided"){
  mp <- margin_params(mom); m <- length(mp); s <- c(0,1)
  p  <- switch(alternative,
               less = c(min(s), 1-q),
               greater = c(q, max(s)) ,
               two.sided = c(q/2, 1-q/2))
  
  if(type=="raw"){
    cr <- data.frame(t(sapply(1:m, function(j){ qbeta(p, mp[[j]]$alpha, mp[[j]]$beta)})))
  }else{
    cr <- data.frame(t(apply(ps, 2, function(x) quantile(x, probs=p)))) 
  }
  colnames(cr) <- c("lower", "upper")
  
  return(cr)
}

get_cr <- function(moms, pss, type="raw", q=0.05, alternative="two.sided"){
  lapply(1:length(moms), function(j) get_cr1(moms[[j]], pss[[j]], type, q, alternative))
}




