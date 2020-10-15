

?bootstrap::jackknife

n <- 100
p <- 0.8
#y <- rbinom(n, 1, p)
y <- c(rep(1, round(p*n)), rep(0, n-round(p*n)))
y <- sample(y)

set.seed(1337)
?DTAmc::sample_data
data <- DTAmc::sample_data(n=100, m=50, prev=c(0.5, 0.5))
data
lapply(data, colMeans)
data <- do.call(rbind, data)
colMeans(data)

phat <- mean(y)
phat

sqrt( phat*(1-phat)/n)

bootstrap::jackknife(y, mean)

?bootstrap::bcanon
bootstrap::bcanon(y, 10000, mean)



#?bootstrap::abcnon
#bootstrap::abcnon(y, 10000, mean)

s <- sample(1:n)
bootstrap::bcanon(1:n, 10000, function(s){
  max(colMeans(data[s, ]))
})

data

DTAmc::dta(data = data.frame(y), group=rep(1,n))

?boot::boot.ci
?boot::boot.return
?boot::boot

?DTAmc::dta

s <- sample(1:nrow(data), replace = TRUE)
S <- numeric(100)
group <- rep(1:2, each=50)
maxMin <- function(data, s, group){
  s <- sort(s)
  #S <<- s
  phat <- lapply(unique(group), function(g) colMeans(data[s, ][group==g, ]) )
  #do.call(rbind, phat)
  #do.call(pmin, phat)
  max(do.call(pmin, phat))
}

teststat <- function(d, p0){
  phat <- colMeans(d)
  stopifnot(length(phat) == length(p0))
  t <- (phat-p0)/sqrt(phat*(1-phat)/nrow(d))
  return(t)
}

#teststat(data[[2]], colMeans(data[[2]]))

maxMinT <- function(data, s, group, 
                    alternative = c("two.sided", "less", "greater")){
  a <- match.arg(alternative)
  f <- switch(a,
              two.sided = abs,
              less = function(x){-x},
              greater = function(x){x})
  s <- sort(s)
  tt <- lapply(unique(group), function(g) {
    (teststat(d = data[s, ][group==g, ], p0=colMeans(data[group==g, ])))
  })
  do.call(pmin, tt)
}

maxMinT1 <- function(data, s, group, 
                     alternative = c("two.sided", "less", "greater")){
  a <- match.arg(alternative)
  f <- switch(a,
              two.sided = abs,
              less = function(x){-x},
              greater = function(x){x})
  s <- sort(s)
  tt <- lapply(unique(group), function(g) {
    (teststat(d = data[s, ][group==g, ], p0=colMeans(data[group==g, ])))
  })
  max(do.call(pmin, tt))
}

coverage <- function(b, v, q, 
                     alternative = c("two.sided", "less", "greater")){
  a <- match.arg(alternative)
  f <- switch(a,
              two.sided = abs,
              less = function(x){-x},
              greater = function(x){x})
  # TODO: extend to two sided
  mean(apply((v) <= b, 1, all)) - q 
}


a = "g"
## (1) bsT
bsT <- boot::boot(data, statistic = maxMinT,
                  R = 10000, strata=group, group=group, alternative=a)

uniroot(coverage, c(1, 10), v=bsT$t, q=0.95, alternative=a)$root
hist(bsT$t)
summary(bsT$t)
coverage(2, bsT$t, 0)

## (2) bsT1
bsT1 <- boot::boot(data, statistic = maxMinT1,
                   R = 10000, strata=group, group=group, alternative=a)

uniroot(coverage, c(1, 10), v=bsT1$t, q=0.95, alternative=a)$root
hist(bsT1$t)
summary(bsT1$t)
coverage(2, bsT1$t, 0)

quantile(1-alpha)

boot::boot.ci(bsT, type="perc")
boot::boot.ci(bsT1, type="perc")

boot::boot.ci(bsT1, type="bca")
summary()

hist(bs$t[,4])

hist(bs$t)
str(bs, 1)
sum(S > 50)

hist(bs$t)
boot::boot.ci

## TODO 0: sample data LFC

## TODO 1: MULTIPLE MCNEMAR
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2902578/pdf/nihms187308.pdf

## TODO 2: BOOSTRAP (BCA/WILD)

## TODO 3: mBeta

data <- DTAmc::sample_data()
DTAmc::dta(data=data, adjustment="bootstrap")



set.seed(123)
M <- as.data.frame(mvtnorm::rmvnorm(20, mean=rep(0, 3), sigma=2*diag(3)))
M
threshold(M)
C <- matrix(rep(c(-1, 0, 1, -2, 0, 2), 3), ncol=3, byrow = TRUE)
C
w <- c(1, 1, 2, 2, 3, 3)
threshold(M, C, w)

M
cu <- c(0, 1, 0, 2, 1, 2)
cu
w <- c(1, 1, 2, 2, 3, 3)
threshold(M, cu, w)




