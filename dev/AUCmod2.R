require(dplyr)
require(corrplot)
require(mvtnorm)

n <- 1000
L <- 10
AUC <- c(0.86, 0.88, 0.9, 0.92)
M <- length(AUC)
mu1 <- sqrt(2)* qnorm(AUC)
mu1
mu0 <- rep(0, M)
#mu1 <- c(1, 1.25, 1.5, 1.75)
rho = 0.5
C <- matrix(rho, M, M)+diag(rep(1-rho,M))
H <- rmvnorm(n, rep(0, M), C) 
D <- rmvnorm(n, mu1, C)


# AUC
pnorm(mu1/sqrt(2))

cu <- seq(0, 2, length.out=100)

sapply(cu, function(t) colMeans(D>t))
sapply(cu, function(t) colMeans(H<=t))




G <- expand.grid(cu, cu)

Y1 <- apply(D, 2, function(col) lapply(cu, function(t) col>t))
Y1 <- do.call(cbind, lapply(Y1, function(x) do.call(cbind, x))) 

Y0 <- apply(H, 2, function(col) lapply(cu, function(t) col<=t))
Y0 <- do.call(cbind, lapply(Y0, function(x) do.call(cbind, x))) 

Se <- colMeans(Y1)
Sp <- colMeans(Y0)
tau <- pmin(Se, Sp+0.0)

plot(1:400, Se)
plot(1:400, Sp)
plot(1:400, tau)

R1 <- cov2cor(cov(Y1))
R0 <- cov2cor(cov(Y0))

set.seed(1337)
S <- 20
sel <- sample(1:400, S, prob=tau^4) %>% sort()
plot(1:S, tau[sel])
R1s <- R1[sel, sel]; R0s <- R0[sel, sel]
corrplot(R1s)
corrplot(R0s)

#replicate_sim




