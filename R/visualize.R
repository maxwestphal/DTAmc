visualize <- function(object, args=list()){
  stopifnot(inherits(object, "DTAmc_result"))
  stopifnot(length(object) == 2)
  plot(object[[1]]$estimate, object[[2]]$estimate)
  # TODO: see preprint figures
}

visualise <- visualize
