

e <- rnorm(10000, 10, 2)
w <- rnorm(10000, 0, 1)
w <- runif(10000, -1, 1)
w <- rbinom(10000, 1, 0.5)-0.5

z <- (e - 10)/2 #rnorm(10000, 2, 0.01) 


quantile(z, 0.975)
hist(z)

sd(z)
mean(z)

quantile(w*z, 0.975)
hist(z*w)

sd(w*z)
mean(w*z)




set.seed(123456)
n <- 100
x <- rbinom(n, 1, 0.75)
mean(x)

?boot::boot
bs <- boot::boot(x, function(x, s) {c(mean(x[s]), sqrt(var(x[s])/n))}, 10000)
bs$t
tt <- (bs$t[,1]-mean(x))/bs$t[,2]
tt <- (bs$t[,1]-mean(x))/sqrt(bs$t[,1]*(1-bs$t[,1])/n)

quantile(tt, 0.975) 
qnorm(0.975)

head(bs$t)

sqrt(0.75*0.25/n)










