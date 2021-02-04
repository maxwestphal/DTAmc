
n1 <- 100
y1 <- 80
a1 <- 1
nu1 <- 2
e1 <- (y1+a1)/(n1+nu1)
s1 <- sqrt((e1*(1-e1))/n1)

n0 <- 100
y0 <- 100
a0 <- 1
nu0 <- 2
e0 <- (y0+a0)/(n0+nu0)
s0 <- sqrt((e0*(1-e0))/n0)

e1
e0


xx <- seq(0, 1, 0.01)


### SINGLE PROP
# single prop wald interval (lower)
qnorm(0.05, e1, s1)
e1 + qnorm(0.05)*s1

dn <- dnorm(xx, e1, s1)
plot(xx, dn, type="l")


# single prop credible interval
qbeta(0.05, y1+a1, n1+nu1-y1-a1)

db <- dbeta(xx, y1+a1, n1+nu1-y1-a1)
plot(xx, db, type="l")


### co-primary endpoint
# CPE normal 
qnorm(0.05, e1, s1)

# CPE mbeta
nrep=100000

r1 <- rbeta(nrep,  y1+a1, n1+nu1-y1-a1)
r0 <- rbeta(nrep,  y0+a0, n0+nu0-y0-a0)

rr <- pmin(r1, r0)

quantile(rr, 0.05)

plot(density(rr))


#####################
set.seed(1337)
dd1 <- sample_data_roc()
dd2 <- sample_data_lfc()

g <- 1
alt <- "two"

## ROC
dta(dd1, adj="none", alt=alt, regu=T)[[1]][[g]]
dta(dd1, adj="bonf", alt=alt, regu=T)[[1]][[g]]
dta(dd1, adj="maxt", alt=alt, regu=T)[[1]][[g]]
dta(dd1, adj="mbeta", alt=alt, regu=T, pars=list(nrep=10000, loss="fwer"))[[1]][[g]]
dta(dd1, adj="mbeta", alt=alt, regu=T, pars=list(nrep=10000, loss="fwer2"))[[1]][[g]]


dta(dd1[1], adj="mbeta", regu=T, loss="fwer")



## LFC


dta(dd2, adj="maxt", regu=T, pars=list(nrep=10000, loss="fwer"))[[1]][[1]]







## CHECK
dd <- dd2

d1 <- dd[[1]]
d2 <- dd[[2]]

j <- 1

n1 <- nrow(d1)
y1 <- colSums(d1)[j]
a1 <- 1
nu1 <- 2
e1 <- (y1+a1)/(n1+nu1)
s1 <- sqrt((e1*(1-e1))/n1)

n0 <- nrow(d2)
y0 <- colSums(d2)[j]
a0 <- 1
nu0 <- 2
e0 <- (y0+a0)/(n0+nu0)
s0 <- sqrt((e0*(1-e0))/n0)

r1 <- rbeta(nrep,  y1+a1, n1+nu1-y1-a1)
r0 <- rbeta(nrep,  y0+a0, n0+nu0-y0-a0)

rr <- pmin(r1, r0)

quantile(rr, 0.05)

p <- sqrt(0.05)

mean(r1 < qbeta(p, a1+y1, n1+nu1-y1-a1))
mean(r0 < qbeta(p, a0+y0, n0+nu0-y0-a0))

mean(r1 < qbeta(p, a1+y1, n1+nu1-y1-a1) & r0 < qbeta(p, a0+y0, n0+nu0-y0-a0))










