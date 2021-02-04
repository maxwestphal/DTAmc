preproc_comp <- function(comparator, data){
  if(is.null(comparator)){
    return(NULL)
  }
  if(comparator == 0){
    return(NULL)
  }
  
  stopifnot(is.numeric(comparator) | is.character(comparator))
  stopifnot(length(comparator) == 1)
  
  m <- ncol(data[[1]])
  modnames <- colnames(data[[1]])
  
  if(is.numeric(comparator)){
    stopifnot(comparator %in% 1:m)
  }
  if(is.character(comparator)){
    stopifnot(comparator %in% modnames)
    comparator <- which(comparator == modnames)
  }
  
  return(comparator) 
}

preproc_regu <- function(regu="2_1_0.5"){
  if(is.logical(regu)){
    if(!regu){
      return(c(0,0,0))
    }
    if(regu){
      return(c(2,1,1/2))
    }
  }
  if(is.character(regu)){
    regu <- as.numeric(strsplit(regu, "_")[[1]])
  }
  stopifnot(is.numeric(regu))
  stopifnot(length(regu)==3)
  stopifnot(all(regu >= 0))
  return(regu)
}