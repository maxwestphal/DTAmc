#' @importFrom magrittr %>%

argmin <- function(x) {
  am <- which(x == min(x))
  if (length(am) > 1) {
    return(sample(am, 1))
  }
  return(am)
}

argmax <- function(x) {
  am <- which(x == max(x))
  if (length(am) > 1) {
    return(sample(am, 1))
  }
  return(am)
}

# matrix_product <- function(...){
#   L <- list(...)
#   stopifnot(length(L) > 0)
#   R <- L[[1]]
#   if(length(L) > 1){
#     for(g in 1:length(L)){
#       R <- R %*% L[[g]]
#     }
#   }
#   return(R)
# }

