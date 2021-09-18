study_dta_mbeta <- function(data = generate_data(seed=1337),
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
  lfp <- ifelse(is.null(pars$lfp), 1, pars$lfp)
  
  G <- length(data)
  m <- ncol(data[[1]])
  type <- attr(contrast(data), "type")
  
  moms <- lapply(data, function(d){
    add_moments(prior_moments(m, regu), data_moments(d))})
  stats <- data2stats(data, contrast, regu, raw=TRUE)
  
  ## posterior sample:
  pss <- sample_mbeta(moms) %>%
    lapply(function(s){s %*% contrast(data)})

  ## credible region:
  qstar <- stats::uniroot(f=eval_cr, interval=c(0, 0.5),
                   moms=moms, pss=pss, type=type,
                   alpha=alpha, alternative=alternative, lfp=lfp)$root
  crs <- get_cr(moms, pss, type=type, q=qstar, alternative=alternative) 
  
  ## output:
  lapply(1:G, function(g) {
    data.frame(
      parameter = stats[[g]]$names,
      alternative = altstr(alternative, benchmark[g]),
      estimate = stats[[g]]$est,
      lower = crs[[g]]$lower,
      upper = crs[[g]]$upper,
      pvalue = NA
    )
  }) %>%
    setattr(
      n = sapply(data, nrow), m=m, 
      alpha=alpha, alpha_adj=qstar, cv=NA,
      class = c("list", "DTAmc_result")
    ) %>% 
    return()
}



# Helper functions ----------------------------------------------------------------------------
eval_cr <- function(q, moms, pss, type, alpha, alternative, lfp=1){
  crs <- get_cr(moms, pss, type, q, alternative)
  coverage(crs, pss, lfp) - (1-alpha)
}

coverage <- function(crs, pss, lfp=1){
  nrep <- nrow(pss[[1]]); m <- ncol(pss[[1]]); G <- length(pss)

  L <- lapply(crs, function(x) matrix(x$lower, nrow=nrep, ncol=m, byrow=TRUE))
  U <- lapply(crs, function(x) matrix(x$upper, nrow=nrep, ncol=m, byrow=TRUE))
  
  C <- lapply(1:G, function(g) covered(pss[[g]], L[[g]], U[[g]]))
  
  mean(apply(Reduce("+", postproc(C, nrep, m, G, lfp)), 1, min) >= 1)
}

postproc <- function(C, nrep, m, G, lfp=1){
  ## 'split prior' approach (1-lfp prior mass normal, lfpr prior mass on 'LFC'):
  # matrix(sample(0:1, nrep*m, replace=TRUE, prob=c(1-lfp, lfp)) * 
  #               sample(1:G, nrep*m, replace=TRUE),
  #             nrow=nrep, ncol=m, byrow=FALSE) %>% 
  #   {lapply(1:G, function(g) {(. %in% c(0, g)) * C[[g]]} )}
  M <- matrix(sample(0:1, nrep*m, replace=TRUE, prob=c(1-lfp, lfp)) * 
           sample(1:G, nrep*m, replace=TRUE),
         nrow=nrep, ncol=m, byrow=FALSE) 
  lapply(1:G, function(g){(M %in% c(0, g)) * C[[g]]} )
}
# TODO: document this in study_dta_mbeta
# TODO: remove old version

covered <- function(S, L, U){
  S >= L & S <= U
}

sample_mbeta <- function(moms, nrep=5000, proj.pd=FALSE){
  lapply(moms, function(mom) sample_mbeta1(mom, nrep, proj.pd))
}

#' @importFrom  copula normalCopula
#' @importFrom  copula P2p
#' @importFrom  copula mvdc
#' @importFrom  copula rMvdc
#' @importFrom  Matrix nearPD
sample_mbeta1 <- function(mom, nrep=5000, proj.pd = FALSE){
  m <- nrow(mom$moments)
  R <- stats::cov2cor(mom2cov(mom))
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

# TODO: remove?!
# transform_sample <- function(sample, tr){
#   lapply(sample, tr)
# }

get_cr1 <- function(mom, ps, type="raw", q=0.05, alternative="two.sided"){
  mp <- margin_params(mom); m <- length(mp); s <- c(0,1)
  p  <- switch(alternative,
               less = c(min(s), 1-q),
               greater = c(q, max(s)) ,
               two.sided = c(q/2, 1-q/2))
  
  if(type=="raw"){
    cr <- data.frame(t(sapply(1:m, function(j){ stats::qbeta(p, mp[[j]]$alpha, mp[[j]]$beta)})))
  }else{
    cr <- data.frame(t(apply(ps, 2, function(x) stats::quantile(x, probs=p)))) 
  }
  colnames(cr) <- c("lower", "upper")
  
  return(cr)
}

get_cr <- function(moms, pss, type="raw", q=0.05, alternative="two.sided"){
  lapply(1:length(moms), function(j) get_cr1(moms[[j]], pss[[j]], type, q, alternative))
}




