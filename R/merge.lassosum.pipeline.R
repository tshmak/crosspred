merge.lassosum.pipeline <- function(..., pthresh=F) {
  #' merge.lassosum.pipeline for lassosum.pipeline objects produced by pthresh.pipeline
  
  result <- lassosum:::merge.lassosum.pipeline(...)
  if(pthresh) {
    #### Fix merge (for pthresh beta) ####
    weird <- result$beta$pthresh
    ncomps <- length(weird) / 4
    for(i in 1:ncomps) {
      if(i == 1) {
        obj <- list(beta=weird[[1]], pbin=weird[[1+ncomps]], 
                    nbin=weird[[1+ncomps*2]], p.thresholds=weird[[1+ncomps*3]])
      } else {
        obj$beta <- c(obj$beta, weird[[i]])
        obj$pbin <- c(obj$pbin, weird[[i+ ncomps]])
      }
    }
    class(obj) <- "pthresh"
    result$beta$pthresh <- obj
  }
  return(result)
}