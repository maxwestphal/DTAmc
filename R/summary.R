#' @export
summary.DTAmc_results <- function(object, ...){ 
  cat("\n")
  cat("<> Diagnostic Test Accuracy <> \n")
  
  cat("\n")
  cat("++ General Information:\n")
  cat("+  Sample sizes:", object$info$n, "\n")
  cat("+  Dimension: ", object$info$m, "\n")
  
  cat("\n")
  cat("++ Analysis: \n")
  cat(paste0("+ Hypotheses: candidate[j] - ", "comparator", " <= ", "benchmark" ,"\n"))
  cat("+  Number of comparisons:", ".....", "\n")
  cat("+  Adjustment:", ".....", "\n")
  
  cat("\n")
  cat("++ Global Result: \n")
  cat("+  ...\n")
  cat("+  ...\n")
  cat("+  ...\n")
  
}

## TODO: summary.DTAmc_results
