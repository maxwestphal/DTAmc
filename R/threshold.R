#' Threshold continuous (bio)markers 
#' 
#' This function allows to split continuous (bio)markers into two or more categories 
#' by specifying one or more cutoffs.
#' 
#' @param markers numeric matrix of continuous (bio)markers to be thresholded.
#' Assume an (n x r) matrix with n observations (subjects) of m continuous markers.
#' @param cutoffs numeric matrix of dimension L x K. Each row of cutoffs defines a split
#' into K+1 distinct categories. Each row must contain distinct values. In the simplest case,
#' cutoffs is a single column matrix whereby is row defines a binary split (<=t vs. >t).
#' In this case (K=1), cutoffs can also be a numeric vector.
#' @param map integer vector of length L with values in 1:r, whereby r = ncol(markers).
#' map_l gives the value which column of markers should be tresholded by ...
#'
#' @return numeric (n x L) matrix with categorical outcomes after thresholding
#' @export
#'
#' @examples
#' set.seed(123)
#' M <- as.data.frame(mvtnorm::rmvnorm(20, mean=rep(0, 3), sigma=2*diag(3)))
#' M
#' threshold(M)
#' C <- matrix(rep(c(-1, 0, 1, -2, 0, 2), 3), ncol=3, byrow = TRUE)
#' C
#' w <- c(1, 1, 2, 2, 3, 3)
#' threshold(M, C, w)
threshold <- function(markers, cutoffs=rep(0, ncol(markers)), map=1:ncol(markers)){
  stopifnot(is.matrix(cutoffs) | is.numeric(cutoffs))
  cutoffs <- as.matrix(cutoffs)
  ## TODO: why needed??
  # if(ncol(cutoffs) == 1){
  #   cutoffs <- do.call(cbind, lapply(1:ncol(markers), function(k) cutoffs))
  # }
  if(!all(apply(cutoffs, 1, function(x) length(x) == length(unique(x))))){
    stop("Marker split cannot be based on duplicate cutoffs.")
  }
  stopifnot(all(map %in% 1:ncol(markers)))
  stopifnot(length(map) == nrow(cutoffs))
  if(is.null(names(markers))){names(markers) <- paste0("marker", 1:ncol(markers))}
  
  C <- as.data.frame(matrix(NA, nrow=nrow(markers), ncol=nrow(cutoffs)))
  for(k in 1:nrow(cutoffs)){
    C[, k] <- threshold1(markers[, map[k]], cutoffs[k, ])
    n <- sum(map[1:k] == map[k])
    a <- ifelse(sum(map == map[k])>1, paste0("_", letters[n]), "")
    names(C)[k] <- paste0(names(markers)[map[k]], "_cat", a)
  }
  
  return(C)
}

threshold1 <- function(x, cutoffs){
  sapply(x, function(xi) sum(xi>cutoffs))
}
