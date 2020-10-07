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

comp_preproc <- function(comp, data){
  # TODO checks
  if(length(comp) == 1){
    comp <- rep(comp, ncol(data[[1]]))
  }
  return(comp) # TODO: implemented comp integer kess
  # TODO: return data? stats? adoptoted?
}


## TODO: REMOVE:
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

# df2list <- function(df, group){
#   lapply(unique(group), function(g) as.data.frame(df)[group==g, ])
# }