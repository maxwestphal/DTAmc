results_mbeta <- function(data, comparator, benchmark,
                          alpha, alternative, transformation, regu, pars){
  message("mbeta adjustment is experimental")
  if(transformation != "none"){
    warning("transformation argument ignored for mbeta adjustment")
    transformation <- "none"
  }
  ## TODO: remove message
  ## TODO: own prior specification via pars$prior = 
  ## TODO: TODO adjust prior params
  prior <- ifelse(is.null(pars$prior), "regu", pars$prior)
  
  G <- length(data)
  m <- ncol(data[[1]])

  
  if(prior == "regu"){
    mom <- prior_moments(m, regu)
    pr <- SIMPle::define_dlist(dlist = lapply(1:G, function(g) 
      SIMPle::define_mBeta(m, nu=mom$n, moments=mom$moments, mode="reduced", msg=FALSE)))
  }
  if(prior == "empirical"){
    pr <- SIMPle::define_mBeta_empirical(data)
  }
  
  po <- pr %>% 
    SIMPle::update(data) %>% 
    SIMPle::inference(loss=ifelse(is.null(pars$loss), "fwer2", pars$loss),
                      target = alpha, 
                      alternative = alternative,
                      method = ifelse(is.null(pars$method), "mcmc", pars$method),
                      n.rep=ifelse(is.null(pars$nrep), 5000, pars$nrep))
  
  cr <- lapply(po, function(x) x$inference)
  est <- lapply(po, function(dist) {
    sapply(SIMPle::margins(dist), function(m) m$features$mean) 
  })
  
  modnames <- paste0("model", 1:m) # TODO
  comparator <- NULL # TODO
  
  lapply(1:G, function(g) {
    cbind(
      hyothesis = hypotheses(modnames, comparator, benchmark[g], alternative),
      estimate = est[[g]],
      cr[[g]]$result,
      pvalue = NA
    )
  })

}


pval_mbeta <- function(){
  NA
}

alpha_mbeta <- function(){
  NA
}