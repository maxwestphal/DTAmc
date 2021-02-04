data2stats <- function(data, regu=c(0,0,0), raw=FALSE){
  lapply(data, dat2stats, regu=regu, raw=raw)
}

dat2stats <- function(dat, regu=c(0,0,0), raw=FALSE){
  
  if(raw){
    est <- dat2est(dat)
    sig <- NA
    se <- dat2se(dat)
  }else{
    mom <- add_moments(prior_moments(ncol(dat), regu), data_moments(dat))
    est <- mom2est(mom)
    sig <- mom2cov(mom)
    se <- cov2se(sig)
  }

  out <- list(est = est, 
              cov = sig,
              se = se,
              n = nrow(dat),
              npseudo = regu[1],
              ntotal = nrow(dat) + regu[1])
}