
#' dta
#'
#' @param data 
#' @param method 
#' @param pars 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dta <- function(data = sample_data(),
                method = "maxt",
                pars=list(),
                ...){
  pars <- c(list(data=data), pars, list(...))
  do.call(paste0("dta_", method), pars)
}

dta_uni <- function(data,
                    correction = c("none", "bonf"),
                    ...){
  return(NULL)
}

dta_maxt <- function(data, ...){
  require(SEPM)
  define_hypothesis("accuracy.cp", threshold=c(0.5, 0.5)) %>% 
    compare(comparison=data) %>% 
    estimate(method="beta.approx") %>% 
    infer(method="maxT") %>% 
    summary()
}

dta_mBeta <- function(data, ...){
  
}

# DEV/TESTING -------------------------------------------------------------


#dta()
