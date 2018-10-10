pthresh.pipeline <- function(beta=cor, pvals, 
                             p.thresholds=c(seq(0.01, 0.99, by=0.01),
                                           10^(seq(-20, -2)),
                                           5*10^(seq(-20, -2))),
                             clump=F, clump.options=NULL, 
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

  #### Checks ####
  if(missing(pvals)) pvals <- attr(beta, "pvals")
  if(is.null("pvals")) stop("pvals must be given.")
  # if(destandardize) stop("destandardization not supported in pthresh.pipeline.")
  
  #### Parse ####
  scale <- max(abs(beta)) * 2 # Trick to bypass correlation checking
  results <- lassosum.pipeline(cor = beta / scale, chr=chr, pos=pos, snp=snp,
                             A1=A1, A2=A2, test.bfile=test.bfile,
                             trace=trace, exclude.ambiguous = exclude.ambiguous,
                             keep.test=keep.test, remove.test=remove.test,
                             destandardize = destandardize,
                             s=numeric(0), ...) # This is just for parsing!!!
  results$sumstats$cor <- results$sumstats$cor * scale
  ss <- results$sumstats

  #### pval.thresh ####
  p.thresholds <- sort(unique(p.thresholds))
  Beta <- ss$cor
  Pvals <- pvals[ss$order]

  #### clumping ####
  if(clump) {
    if(trace > 0) {
      cat("Performing clumping...\n")
    }
    opts <- clump.options
    opts$bfile <- results$ref.bfile
    opts$keep <- results$keep.ref
    opts$cluster <- cluster
    ref.bim <- lassosum:::read.table2(paste0(opts$bfile, ".bim"))
    if(results$test.bfile != results$ref.bfile) {
      test.bim <- lassosum:::read.table2(paste0(results$test.bfile, ".bim"))
      ss.test.bim <- test.bim[results$test.extract,]
      ss.test.bim$Pvals <- Pvals
      m <- lassosum:::matchpos(tomatch = ss.test.bim, ref.df = ref.bim, 
                               auto.detect.tomatch = F,auto.detect.ref = F,
                               chr = 'V1',ref.chr = 'V1',
                               pos = "V4",ref.pos = "V4",snp = "V2", ref.snp = "V2",
                               ref = 'V6', ref.ref='V6', 
                               exclude.ambiguous = results$exclude.ambiguous)
      ss.test.bim$order <- Inf
      ss.test.bim$order[m$order] <- m$order
      opts$pvals <- Pvals[m$order]
      opts$extract <- m$ref.extract
    } else {
      test.bim <- ref.bim
      ss.test.bim <- test.bim[results$test.extract,]
      ss.test.bim$Pvals <- Pvals
      ss.test.bim$order <- 1:nrow(ss.test.bim)
      opts$pvals <- Pvals
      opts$extract <- results$test.extract
    }
    opts$trace <- trace - 1
    tab <- do.call(plink.clump, opts)
    ss.test.bim$toexclude <- ss.test.bim$order < Inf & !(ss.test.bim$V2 %in% tab$SNP)
    Pvals[ss.test.bim$toexclude] <- 1
    Beta[ss.test.bim$toexclude] <- 0
  }
  
  #### destandardize ####
  if(destandardize) {
    sd <- sd.bfile(bfile=test.bfile, keep=results$keep.test, 
                   extract = results$test.extract)
    Beta <- Beta / sd
    Beta[is.infinite(Beta)] <- 0
  }
  pt <- pval.thresh(beta=Beta, pvals = Pvals,
                    p.thresholds = p.thresholds,
                    bfile=test.bfile, keep = keep.test,
                    extract=results$test.extract,
                    cluster = cluster, trace=trace-1)

  #### Save results ####
  results$beta <- list(pthresh=attr(pt, "beta"))
  results$pgs <- list(pthresh=pt[,]) # [,] for de-attributing
  results$time <- (proc.time() - time.start)["elapsed"]
  results$lambda <- attr(pt, "beta")$p.thresholds
  results$s <- "pthresh"

  #' @return A \code{lassosum.pipeline} object
  class(results) <- "lassosum.pipeline"
  return(results)

}
