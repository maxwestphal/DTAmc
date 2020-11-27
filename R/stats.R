data2stats <- function(data, regu=c(0,0,0)){
  lapply(data, dat2stats, regu=regu)
}

dat2stats <- function(dat, regu=c(0,0,0)){
  mom <- add_moments(prior_moments(ncol(dat), regu),
                     data_moments(dat))
  sigma <- mom2cov(mom)
  list(est = mom2est(mom), 
       cov = sigma,
       se = cov2se(sigma),
       n = mom$n,
       nraw = nrow(dat),
       npseudo = regu[1])
}

## TODO: allow to NOT calc covariance matrix if not needed
## TODO: REMOVE OLD VERSION
# data2stats <- function(data, regu=c(0,0,0), covariance=TRUE){
#   
#   n <- sapply(data, nrow)
#   
#   ## case 1: no regularisation
#   if(regu[1] == 0){
#     est <- lapply(data, dat2est)
#     if(covariance){
#       sigma <- lapply(data, dat2cov)
#       se <- lapply(sigma, cov2se)
#     }
#     if(!covariance){
#       sigma <- NULL
#       se <- lapply(data, dat2se)
#     }
#   }
#   ## case 2: enabled regularisation
#   ## TODO: simpify covariance = FALSE case
#   if(regu[1] > 0){
#     mom <- lapply(data, function(d){
#       add_moments(prior_moments(ncol(d), regu),
#                   data_moments(d))
#     })
#     est <- lapply(mom, mom2est)
#     sigma <- lapply(mom, mom2cov)
#     se <- lapply(sigma, cov2se)
#   }
#   
#   
#   out <- lapply(1:length(data), function(g){
#     list(est = est[[g]],
#          se = se[[g]],
#          cov = sigma[[g]],
#          n = n[g])
#   })
#   names(out) <- names(data)
#   return(out)
# }