# single classifier, m = 1

se <- 0.8
sp <- 0.8
se0 <- 0.8
sp0 <- 0.8
n1 <- 50
n0 <- 150

set.seed(1337)
y1 <- rbinom(n1, 1, se)
y0 <- rbinom(n0, 1, sp)

a <- 0; b <- 0

e1 <- (sum(y1)+a)/(n1+a+b)
e0 <- (sum(y0)+a)/(n0+a+b)

s1 <- sqrt(e1*(1-e1)/n1)
s0 <- sqrt(e0*(1-e0)/n0)

t1 <- (e1-se0)/s1
t0 <- (e0-sp0)/s0

#b <- e1 < e0
tm <- min(t1, t0)

## wild B
?DTAmc::sample_data_lfc
data <- DTAmc::sample_data_lfc(m=5)
data

wildbs <- function(data, nboot=5000){
  m <- ncol(data[[1]])
  G <- length(data)
  n <- sapply(data, nrow)
  W <- lapply(1:G, function(g){
    # TODO: replace with generic weight function
    do.call(abind::abind, c(
      lapply(1:nboot, function(b){
        matrix(rnorm(n), nrow=n[g], ncol=m, byrow=FALSE)
      }),
      along=3))
  })
  
  M <- lapply(1:G, function(g) matrix(colMeans(data[[g]]), nrow=n[g], ncol=m, byrow=T))
  
  Z <- lapply(1:G, function(g){
    do.call(abind::abind, c(
      lapply(1:nboot, function(b){
        W[[g]][, , b] * (data[[g]] - M[[g]])
      }),
      along=3))
  })
  
  D <- lapply(1:G, function(g){
    apply(Z[[g]], 2:3, function(x){mean(x)}) 
  })
  
  E <- lapply(1:G, function(g){
    apply(Z[[g]], 2:3, function(x){sqrt(var(x)/n[g])}) 
  })
  
  S <- lapply(1:G, function(g){D[[g]]/E[[g]]})
  
  R <- Reduce(function(a, b) pmin(a, b, na.rm=T), S)
  
  A <- apply(R, 2, max) ## TODO: extend to stepwise??
  #hist(R[4,])
  #hist(A)
  
  return(A)
  # Z[[g]] %>% str()
  # Z[[g]][1:4, 1:5, 1:3]
  # W[[g]][1:3, 1:5, 1:3]
  
}
wildbs(data)


z1 <- w1*(y1-e1)
z0 <- w0*(y0-e0)


# 
t1 <- (mean(z1)-sd(z1)/n1)
t0 <- (mean(z0)-sd(z1)/n0)

bm <- pmin(b1, b0)

b0

### HAT matrix
m = 5
n = 50

#predictor matrix (wrong: wrong n and false iid)
X <- matrix(0, n*m, m)
for(j in 1:m){
  I <- ((j-1)*n+1):(j*n)
  X[I, j] <- 1
}
X <- matrix(1, nrow=n, ncol=1)
#X <- cbind(X, c(1,1,0))
X


H <- X %*% solve(t(X) %*% X) %*% t(X)
diag(H)

y <- rbinom(n*m, 1, 0.75)

Q <- matrix(y, ncol=m)
colMeans(Q)

solve(t(X) %*% X) %*% t(X) %*% Q

## cov of beta hat
solve(t(X) %*% X)

## other formulation
y <- matrix(rep(1, n), ncol=1)
X <- Q

solve(t(X) %*% X) %*% t(X) %*% y

solve(t(X) %*% X)


## experiments
qnorm(0.95)


qnorm(0.025)
rnorm(1000000) %>% abs %>% quantile(1-0.05)

qnorm(0.025)


d <- c(1000, 10, 100000); array(1:prod(d), dim=d)




