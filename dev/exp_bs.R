# 
# 
# ?bootstrap::jackknife
# 
# n <- 100
# p <- 0.8
# #y <- rbinom(n, 1, p)
# y <- c(rep(1, round(p*n)), rep(0, n-round(p*n)))
# y <- sample(y)
# 
# set.seed(1337)
# ?DTAmc::sample_data
# data <- DTAmc::sample_data(n=100, m=50, prev=c(0.5, 0.5))
# data
# lapply(data, colMeans)
# data <- do.call(rbind, data)
# colMeans(data)
# 
# phat <- mean(y)
# phat
# 
# sqrt( phat*(1-phat)/n)
# 
# bootstrap::jackknife(y, mean)
# 
# ?bootstrap::bcanon
# bootstrap::bcanon(y, 10000, mean)
# 
# 
# 
# #?bootstrap::abcnon
# #bootstrap::abcnon(y, 10000, mean)
# 
# s <- sample(1:n)
# bootstrap::bcanon(1:n, 10000, function(s){
#   max(colMeans(data[s, ]))
# })
# 
# data
# 
# DTAmc::dta(data = data.frame(y), group=rep(1,n))
# 
# ?boot::boot.ci
# ?boot::boot.return
# ?boot::boot
# 
# ?DTAmc::dta
# 
# s <- sample(1:nrow(data), replace = TRUE)
# S <- numeric(100)
# group <- rep(1:2, each=50)
# maxMin <- function(data, s, group){
#   s <- sort(s)
#   #S <<- s
#   phat <- lapply(unique(group), function(g) colMeans(data[s, ][group==g, ]) )
#   #do.call(rbind, phat)
#   #do.call(pmin, phat)
#   max(do.call(pmin, phat))
# }
# 
# 




######################################
# pmin <- function(bb, pp){
#   bb <- apply(unique(bb))
#   stopifnot(length(bb))
# }

# data2dd <- function(data, s, group){
#   lapply(unique(group), function(g) data[s, ][group==g, ])
# }

# teststat <- function(phat, p0, n){
#   (phat-p0)/sqrt(phat*(1-phat)/n)
# }



argmin <- function(x, det=TRUE){
  am <- which(x == min(x))
  # deterministic output, if required, otherwise randomize in case of ties
  if(det){ 
    return(min(am))
  }
  if(length(am) == 1){
    return(am)
  }else{
    return(sample(am, 1))
  }
}

pargmin <- function(..., a=NULL, det=TRUE){
  if(is.null(a)){
    a <- list(...)
  }
  stopifnot(do.call(all.equal, lapply(a, length)))
  d <- as.data.frame(a)
  apply(d, 1, argmin, det=det)
}

mstar <- function(d, comp=0){
  do.call(pmin, list(1:5, 5:1))
}

b2bb <- function(b){
  lapply(sort(unique(b)), function(x) as.numeric(b==x))
}

minStats <- function(data, s, group){
  dd <- lapply(unique(group), function(g) data[s, ][group==g, ])
  pp <- lapply(dd, colMeans) ## TODO: replace with data2stats? SUBSTRACT cutoff?
  b <- pargmin(args=pp)
  bb <- b2bb(b)
  p <- apply(mapply('*', bb, pp), 1, sum)
  nn <- sapply(dd, nrow)
  n <- nn[b]
  data.frame(p=p, n=n)
}

maxMinT <- function(data, s, group){
  #                  alternative = c("two.sided", "less", "greater")){
  #a <- match.arg(alternative)
  # f <- switch(a,
  #             two.sided = abs,
  #             less = function(x){-x},
  #             greater = function(x){x}) ## TODO: ??????????????
  
  s <- sort(s)

  st <- minStats(data, s, group)
  st0 <- minStats(data, 1:nrow(data), group)

  max((st$p-max(st0$p))/sqrt(st$p*(1-st$p)/st$n))
}


dat <- DTAmc::sample_data_lfc()

group <- as.numeric(sapply(1:length(data), function(g) rep(g, nrow(data[[g]]))))
data <- do.call(rbind, dat)

#s <- 1:nrow(data)
s <- as.numeric(sapply(sort(unique(group)), function(g){
  sample( (1:nrow(data))[group==g], sum(group==g), replace=TRUE )
}))


maxMinT(data, s, group)

bsT <- boot::boot(data, statistic = maxMinT,
                  R = 5000, strata=group, group=group)


hist(bsT$t)

quantile(bsT$t, c(0.95, 0.975))

boot::boot.ci(bsT, type="perc")
boot::boot.ci(bsT, type="bca")



DTAmc::dta(dat, adjustment = "maxt", regu = T)
dat


devtools::install_github("fdhidalgo/multitestr")


mintest <- function(n0 = 50, n1=100, p0=0.8, p1=0.9, c0=0.8, c1=0.8){
  e0 <- rbinom(1, n0, p0)/n0
  e1 <- rbinom(1, n1, p1)/n1
  t0 <- (e0-c0)/sqrt(e0*(1-e0)/n0)
  t1 <- (e1-c1)/sqrt(e1*(1-e1)/n1)
  c(t1, t0)[which.min(c(e1, e0))]
}

ttt <- replicate(100000, mintest(10000, 20000))

plot(density(ttt))
mean(ttt); sd(ttt)


## TODO 0: sample data LFC

## TODO 1: MULTIPLE MCNEMAR
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2902578/pdf/nihms187308.pdf

## TODO 2: BOOSTRAP (BCA/WILD)

## TODO 3: mBeta


################################################################################
# maxMinT1 <- function(data, s, group, 
#                      alternative = c("two.sided", "less", "greater")){
#   a <- match.arg(alternative)
#   f <- switch(a,
#               two.sided = abs,
#               less = function(x){-x},
#               greater = function(x){x})
#   s <- sort(s)
#   tt <- lapply(unique(group), function(g) {
#     (teststat(d = data[s, ][group==g, ], p0=colMeans(data[group==g, ])))
#   })
#   max(do.call(pmin, tt))
# }

# coverage <- function(b, v, q, 
#                      alternative = c("two.sided", "less", "greater")){
#   a <- match.arg(alternative)
#   f <- switch(a,
#               two.sided = abs,
#               less = function(x){-x},
#               greater = function(x){x})
#   # TODO: extend to two sided
#   mean(apply((v) <= b, 1, all)) - q 
# }
# 
# 
# a = "g"
# ## (1) bsT
# 
# 
# uniroot(coverage, c(1, 10), v=bsT$t, q=0.95, alternative=a)$root
# hist(bsT$t)
# summary(bsT$t)
# coverage(2, bsT$t, 0)
# 
# ## (2) bsT1
# bsT1 <- boot::boot(data, statistic = maxMinT1,
#                    R = 10000, strata=group, group=group, alternative=a)
# 
# uniroot(coverage, c(1, 10), v=bsT1$t, q=0.95, alternative=a)$root
# hist(bsT1$t)
# summary(bsT1$t)
# coverage(2, bsT1$t, 0)
# 
# quantile(1-alpha)
# 
# 
# boot::boot.ci(bsT1, type="perc")
# 
# boot::boot.ci(bsT1, type="bca")
# summary()
# 
# hist(bs$t[,4])
# 
# hist(bs$t)
# str(bs, 1)
# sum(S > 50)
# 
# hist(bs$t)
# boot::boot.ci
# 

# 
# data <- DTAmc::sample_data()
# DTAmc::dta(data=data, adjustment="bootstrap")
# 
# 
# 
# set.seed(123)
# M <- as.data.frame(mvtnorm::rmvnorm(20, mean=rep(0, 3), sigma=2*diag(3)))
# M
# threshold(M)
# C <- matrix(rep(c(-1, 0, 1, -2, 0, 2), 3), ncol=3, byrow = TRUE)
# C
# w <- c(1, 1, 2, 2, 3, 3)
# threshold(M, C, w)
# 
# M
# cu <- c(0, 1, 0, 2, 1, 2)
# cu
# w <- c(1, 1, 2, 2, 3, 3)
# threshold(M, cu, w)




