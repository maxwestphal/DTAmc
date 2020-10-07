

r <- c(0.4, 0.6, -0.2)
R <- matrix(c(1, r[1:2], r[1], 1, r[3], r[2:3], 1), 3, ,3)



d <- mvtnorm::rmvnorm(100, sigma=R)

corrplot::corrplot(cov2cor(cov(d)))
