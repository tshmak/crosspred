cp.pthresh <- function(cp.plink.linear, pvals=NULL, ...) {
  #' @title p-value thresholding with cross-prediction
  #' @description A wrapper around \code{\link{cp.lassosum}}
  #' @param cp.plink.linear A \code{\link{cp.plink.linear}} object
  #' @param ... Other parameters to pass to lassosum.pipeline()
  #' @export

  cp <- cp.plink.linear # shorten 
  
  if(is.null(pvals)) {
    for(i in 1:length(cp$cor)) {
      pvals <- cor2p(cp$cor[[i]], n=cp$nonmiss[[i]])
      pvals[is.na(pvals)] <- 1
      attr(cp$cor[[i]], "pvals") <- pvals
    }
  } else {
    stopifnot(is.list(pvals)) 
    stopifnot(length(pvals) == length(cp$cor))
    stopifnot(all(sapply(cp$cor, length) == sapply(pvals, length)))
    for(i in 1:length(cp$cor)) {
      attr(cp$cor[[i]], "pvals") <- pvals[[i]]
    }
  }
  result <- cp.lassosum(cp, pipeline.FUN = pthresh.pipeline, ...)
  

  class(result) <- "cp.lassosum"
  
  return(result)
  #' @return A \code{cp.lassosum} object 

}

