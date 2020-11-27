

set.seed(NULL)
n1 <- 10
n0 <- 90
p <- 0.8
L <- 0.9


nsim <- 100000

n1t <- n1+2; n0t <- n0+2
x1 <- (rbinom(nsim, n1, p)+1)/n1t
x0 <- (rbinom(nsim, n0, L)+1)/n0t

t1 <- (x1-p)/sqrt(x1*(1-x1)/n1t)
t0 <- (x0-p)/sqrt(x0*(1-x0)/n0t)

p1 <- 1-pnorm(t1)
p0 <- 1-pnorm(t0)

b <- pargmin(x1, x0)

tm <- pmin(t1, t0)
tb <- sapply(1:length(b), function(i) c(t1[i], t0[i])[b[i]])

pm <- pmax(p1, p0)
pb <- sapply(1:length(b), function(i) c(p1[i], p0[i])[b[i]])

qnorm(0.95)
#quantile(t1, 0.95)
mean(tm>qnorm(0.95))
mean(tb>qnorm(0.95))
mean(pm < 0.05)
mean(pb < 0.05)

hist(tm)
hist(tb)

table(b)

all(tm-tb==0)

hist(tm)
hist(tm-tb)

plot(tm, tb)

any(!is.finite(t1) & !is.finite(t0))
sum(!is.finite(t0))

sum(x1==1)
