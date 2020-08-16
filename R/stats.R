# in: data in list format and further args
# out: stats = list(est, cov)

data2stats <- function(data, regu=FALSE, covariance=FALSE){
  regu <- preproc_regu(regu)
  ## case 1: no regularisation
  if(regu[1] == 0){
    est <- lapply(data, dat2est)
    if(covariance){
      sigma <- lapply(data, dat2cov)
      se <- lapply(sigma, cov2se)
    }
    if(!covariance){
      sigma <- NULL
      se <- lapply(data, dat2se)
    }
  }
  ## case 2: enabled regularisation
  ## TODO: simpify covariance = FALSE case
  if(regu[1] > 0){
    mom <- lapply(data, function(d){
      add_moments(prior_moments(ncol(d), regu),
                  data_moments(d))
    })
    est <- lapply(mom, mom2est)
    sigma <- lapply(mom, mom2cov)
    se <- lapply(sigma, cov2se)
  }
  return(list(est=est, se=se, sigma=sigma))
}

preproc_regu <- function(regu="2_1_0.5"){
  if(is.logical(regu)){
    if(!regu){
      return(c(0,0,0))
    }
    if(regu){
      return(c(2,1,1/2))
    }
  }
  if(is.character(regu)){
    regu <- as.numeric(strsplit(regu, "_")[[1]])
  }
  stopifnot(is.numeric(regu))
  stopifnot(length(regu)==3)
  stopifnot(all(regu >= 0))
  return(regu)
}

#preproc_regu(regu=c(1,1,1, 1))
