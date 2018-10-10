pval.thresh <- function(pvals, p.thresholds, beta, bfile,
                        keep=NULL, remove=NULL, extract=NULL, exclude=NULL,
                        chr=NULL, cluster=NULL, trace=0) {
  #' Fast way to do p-value thresholding (without looping over the thresholds)
  #' @export
  beta <- as.vector(beta)
  stopifnot(is.vector(pvals) & is.vector(beta))
  stopifnot(length(pvals) == length(beta))
  Pvals <- sort(unique(c(0,p.thresholds,1)))
  cut <- cut(pvals, Pvals, include.lowest = TRUE)
  stopifnot(!any(is.na(cut)))
  nlevels <- nlevels(cut)

  parsed <- parseselect(bfile, extract, exclude, keep, remove, chr)
  pbin <- as.integer(cut) - 1

  obj <- list(beta=beta, pbin=pbin, nbin=nlevels, p.thresholds=Pvals[-1])
  class(obj) <- "pthresh"

  result <- pgs(bfile=bfile, weights = obj, keep=parsed$keep, extract=parsed$extract,
                cluster=cluster, trace=trace)

  attr(result, "beta") <- obj

  return(result)

}
