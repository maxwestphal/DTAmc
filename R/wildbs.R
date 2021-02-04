# TODO: OLD/NOT NEEDED?

# wildbs <- function(data, nboot=5000){
#   m <- ncol(data[[1]])
#   G <- length(data)
#   n <- sapply(data, nrow)
#   W <- lapply(1:G, function(g){
#     # TODO: replace with generic weight function
#     do.call(abind::abind, c(
#       lapply(1:nboot, function(b){
#         matrix(rnorm(n), nrow=n[g], ncol=m, byrow=FALSE)
#       }),
#       along=3))
#   })
#   
#   M <- lapply(1:G, function(g) matrix(colMeans(data[[g]]), nrow=n[g], ncol=m, byrow=T))
#   
#   Z <- lapply(1:G, function(g){
#     do.call(abind::abind, c(
#       lapply(1:nboot, function(b){
#         W[[g]][, , b] * (data[[g]] - M[[g]]) # TODO: or substract max
#       }),
#       along=3))
#   })
#   
#   D <- lapply(1:G, function(g){
#     apply(Z[[g]], 2:3, function(x){mean(x)}) 
#   })
#   
#   E <- lapply(1:G, function(g){
#     apply(Z[[g]], 2:3, function(x){sqrt(var(x)/n[g])}) 
#   })
#   
#   S <- lapply(1:G, function(g){ (D[[g]]-0) /E[[g]]}) # what to place instead of 0?
#   
#   R <- Reduce(function(a, b) pmin(a,b, na.rm=T), S)
#   
#   A <- apply(R, 2, max) ## TODO: extend to stepwise??
#   #hist(R[4,])
#   hist(A)
#   
#   return(A)
#   # Z[[g]] %>% str()
#   # Z[[g]][1:4, 1:5, 1:3]
#   # W[[g]][1:3, 1:5, 1:3]
#   
# }
