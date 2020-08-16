# DEV/TESTING -------------------------------------------------------------
prior_moments(5, pars=c(5, 3, 2))

#mom <- dm

set.seed(1337)
data = sample_data(m=5)
d <- data[[1]]
d
n <- nrow(d)
0.75^2
pm <- prior_moments(ncol(d), 
                    pars= c(0,0,0))
#pars=c(2,1,0.5 ))
#pars=c(4,3,9/4 ))
dm <- data_moments(d)
dm
#dm[[2]] == data_moments2(d)
pm
dm
mapply("+", pm, dm, SIMPLIFY = FALSE)

C <- mom2cov(add_moments(pm, dm))
C
Q <- cov(d)/(n) *(n-1)/(n+1)
Q
cov2cor(C)
cov2cor(Q)

(C-Q)/C

mvtnorm::qmvnorm(0.95, sigma=cov2cor(C))
mvtnorm::qmvnorm(0.95, sigma=cov2cor(Q))

## VARIANCE 1: (N-1) vs N denominator
## VARIANCE 2: raw var vs sample mean

mom2cov(data_moments(d))
C / cov(d)
mom2cov(data_moments(d)) == cov(d)

?dta
dta(data, method="maxt")
