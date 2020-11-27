cv_uni <- function(alpha = 0.05,
                   alternative = "two.sided",
                   ...) {
  c(switch(
    alternative,
    two.sided = qnorm(alpha / 2),
    less = -Inf,
    greater = qnorm(alpha)
  ),
  switch(
    alternative,
    two.sided = qnorm(1 - alpha / 2),
    less = qnorm(1 - alpha),
    greater = Inf
  )
  )
}

pval_uni <- function(tstat, alternative = "two.sided"){
  switch(
    alternative,
    two.sided = 2*pnorm(abs(tstat), lower.tail = FALSE),
    less = pnorm(abs(tstat), lower.tail = FALSE),
    greater = pnorm(abs(tstat), lower.tail = FALSE)
  )
}




