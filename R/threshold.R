#' Categorize continuous (bio)markers 
#' 
#' This function allows to split continuous (bio)markers into two or more categories 
#' by specifying one or more cutoff values.
#' 
#' @param markers numeric matrix of continuous (bio)markers to be categorized.
#' Assume an (n x r) matrix with n observations (subjects) of r continuous markers.
#' @param cutoffs numeric matrix of dimension r x k. Each row of cutoffs defines a split
#' into k+1 distinct categories. Each row must contain distinct values. In the simplest case,
#' cutoffs is a single column matrix whereby is row defines a binary split (<=t vs. >t).
#' In this case (k=1), cutoffs can also be a numeric vector.
#' @param map integer vector of length k with values in 1:r, whereby r = ncol(markers).
#' map_l gives the value which column of markers should be categorized by ...
#'
#' @return numeric (n x k) matrix with categorical outcomes after categorizing.
#' @export
#'
#' @examples
#' set.seed(123)
#' M <- as.data.frame(mvtnorm::rmvnorm(20, mean=rep(0, 3), sigma=2*diag(3)))
#' M
#' categorize(M)
#' C <- matrix(rep(c(-1, 0, 1, -2, 0, 2), 3), ncol=3, byrow = TRUE)
#' C
#' w <- c(1, 1, 2, 2, 3, 3)
#' categorize(M, C, w)
categorize <- function(markers,
                      cutoffs=rep(0, ncol(markers)),
                      map=1:ncol(markers)){
  stopifnot(is.matrix(cutoffs) | is.numeric(cutoffs))
  cutoffs <- as.matrix(cutoffs)
  if(!all(apply(cutoffs, 1, function(x) length(x) == length(unique(x))))){
    stop("Marker split cannot be based on duplicate cutoffs.")
  }
  stopifnot(all(map %in% 1:ncol(markers)))
  stopifnot(length(map) == nrow(cutoffs))
  markers <- as.data.frame(markers)
  if(is.null(names(markers))){names(markers) <- paste0("marker", 1:ncol(markers))}
  
  C <- as.data.frame(matrix(NA, nrow=nrow(markers), ncol=nrow(cutoffs)))
  for(k in 1:nrow(cutoffs)){
    C[, k] <- categorize1(markers[, map[k]], cutoffs[k, ])
    if(ncol(cutoffs)==1){
      a <- as.character(cutoffs[k, ])
    }else{
      a <- letters[sum(map[1:k] == map[k])]
    }
    names(C)[k] <- paste0(names(markers)[map[k]], "_", a)
  }
  return(C)
}

categorize1 <- function(x, cutoffs){
  sapply(x, function(xi) sum(xi>cutoffs))
}
