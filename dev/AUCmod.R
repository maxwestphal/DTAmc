require(dplyr)
require(corrplot)

n <- 10000
L <- 10
mu1 <- 1.5
D <- rnorm(n, mu1, 1)
H <- rnorm(n, 0, 1)

# AUC
pnorm(mu1/sqrt(2))

cu <- seq(0.5, 1.5, length.out=L)

sapply(cu, function(t) mean(D>t))
sapply(cu, function(t) mean(H>t))



G <- expand.grid(cu, cu)

U1 <- matrix(apply(G, 1, function(t) mean(D>t[1] & D>t[2])), L, L)
u1 <- diag(U1)

U0 <- matrix(apply(G, 1, function(t) mean(H<=t[1] & H<=t[2])), L, L)
u0 <- diag(U0)

S1 <- (n*U1 - u1 %*%t(u1))/n^3
S0 <- (n*U0 - u0 %*%t(u0))/n^3

R1 <- cov2cor(S1) 
R0 <- cov2cor(S0)

R0
cbind(cu=cu, Se=u1, Sp=u0)
corrplot(R1)
corrplot(R0)

