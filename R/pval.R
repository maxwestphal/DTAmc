pval <- function(tstat, alt, adjustment="none"){
  ifelse(alt=="two.sided", 2, 1) * pnorm(tstat, lower.tail=FALSE)
}

## TODO: not adjusted
## TODO: functions for different adjustment