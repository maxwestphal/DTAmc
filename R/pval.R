pval <- function(tstat, alt, adjustment="none"){
  p <- ifelse(alt=="two.sided", 2, 1) * pnorm(tstat, lower.tail=FALSE)
  return(round(p, 4))
}

## TODO: round
## TODO: not adjusted
## TODO: functions for different adjustment