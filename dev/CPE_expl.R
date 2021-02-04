

gen_lfc <- function(se0=0.8, sp0=0.8, nrep=10000){
  bb <- rbinom(nrep, 1, 0.5)
  se <- runif(nrep, se0, 1)
  sp <- runif(nrep, sp0, 1)
  data.frame(se = bb*se0 + (1-bb)*se, sp= (1-bb)*sp0 + bb*sp)
}

gen_lfc(0.2, 0.9) %>% plot()
gen_lfc(0.2, 0.1) %>% cor()

require(SIMPle)

define_mBeta(2, 5, c(0.7, 0.8), corr=-0.5) %>% visualize()


## MIN DIFF
gen_md <- function(am=4, bm=2, ad=1, bd=1, nrep=10000){
  m <- rbeta(nrep, am, bm)
  w <- rbeta(nrep, ad, bd)
  # d = se-sp
  dl <- -(1-m); du = 1-m
  d <- dl + w*(du-dl)
  se <- (d > 0)*(m+d) + (d<=0)*m
  sp <- (d < 0)*(m-d) + (d>=0)*m
  data.frame(se=se, sp=sp)
}

xx <- seq(0, 1, 0.01)
plot(xx, dbeta(xx, 4, 2), type="l")

pbeta(0.8, 4, 2)

gen_md(30, 10, 5, 5) %>% plot()

require(ggplot2)

gen_md(40, 20, 0.25, 0.25) %>% 
  ggplot(aes(x=se, y=sp)) +
  geom_point(alpha=0.25) +
  geom_density2d()


library(dplyr)
library(ggplot2)

util <- function(se, sp, w=c(1,1)){
  se^w[1] * sp^w[2]
}

df <- expand.grid(se = seq(0, 1, 1/100),
                  sp = seq(0, 1, 1/100)) %>%
  mutate(u = util(se, sp, -10*c(1,1)/2)) 

df %>%
  ggplot(aes(se, sp, z=u)) +
  geom_contour(lwd=1.1)


















