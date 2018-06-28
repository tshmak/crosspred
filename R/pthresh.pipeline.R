pthresh.pipeline <- function(beta=cor, pvals,
                             p.thresholds=c(seq(0.05, 0.99, by=0.01),
                                           10^(seq(-20, -2)),
                                           5*10^(seq(-20, -2))),
                             chr=NULL, pos=NULL, snp=NULL,
                             A1=NULL, A2=NULL,
                             test.bfile=NULL,
                             trace=1,
                             exclude.ambiguous=TRUE,
                             keep.test=NULL, remove.test=NULL,
                             cluster=NULL, 
                             destandardize=FALSE, 
                             cor=NULL, ...) {
  #' @title Run p-value thresholding to create a lassosum.pipeline object
  #' @description p-value thresholding for lassosum
  #'
  #' @export
  #'

  time.start <- proc.time()

  ### Checks ###
  if(missing(pvals)) pvals <- attr(beta, "pvals")
  if(is.null("pvals")) stop("pvals must be given.")
  if(destandardize) stop("destandardization not supported in pthresh.pipeline.")
  
  ### Parse ###
  results <- lassosum.pipeline(cor = beta, chr=chr, pos=pos, snp=snp,
                             A1=A1, A2=A2, test.bfile=test.bfile,
                             trace=trace, exclude.ambiguous = exclude.ambiguous,
                             keep.test=keep.test, remove.test=remove.test,
                             destandardize = destandardize,
                             s=numeric(0), ...)
  ss <- results$sumstats

  ### pval.thresh ###
  p.thresholds <- sort(unique(p.thresholds))
  Beta <- ss$cor
  Pvals <- pvals[ss$order]

  pt <- pval.thresh(beta=Beta, pvals = Pvals,
                    p.thresholds = p.thresholds,
                    bfile=test.bfile, keep = keep.test,
                    extract=results$test.extract,
                    cluster = cluster, trace=trace)

  ### Save results ###
  results$beta <- list(pthresh=attr(pt, "beta"))
  results$pgs <- list(pthresh=pt[,]) # [,] for de-attributing
  results$time <- (proc.time() - time.start)["elapsed"]
  results$lambda <- attr(pt, "beta")$p.thresholds
  results$s <- "pthresh"

  #' @return A \code{lassosum.pipeline} object
  class(results) <- "lassosum.pipeline"
  return(results)

}
