#' @title Internal engine for polygenic scores with p-value thresholding
#'
#' @param bfile A plink bfile stem
#' @param weights A pthresh object!!! (I use this name to trick validate.lassosum.pipeline)
#' @param extract SNPs to extract (see \code{\link{parseselect}})
#' @param exclude SNPs to exclude (see \code{\link{parseselect}})
#' @param keep samples to keep (see \code{\link{parseselect}})
#' @param remove samples to remove (see \code{\link{parseselect}})
#' @param chr a vector of chromosomes
#' @param cluster A \code{cluster} object from the \code{parallel} package.
#' For parallel processing.
#'

#' @export
.pgs.pthresh <- function(weights, bfile, keep=NULL, extract=NULL, exclude=NULL, remove=NULL,
                   chr=NULL, cluster=NULL, trace=0) {

  beta <- weights$beta
  pbin <- weights$pbin
  nbin <- weights$nbin

  if(length(beta) != length(pbin)) {
    stop("Length of pbin doesn't match length of beta")
  }
  stopifnot(max(pbin) < nbin)

  if(length(bfile) > 1) {
    pgs.vec(bfile=bfile, weights=weights, extract=extract, exclude=exclude, 
            keep=keep, remove=remove, chr=chr, cluster=cluster, trace=trace)
  }

  stopifnot(is.numeric(beta))
  stopifnot(!any(is.na(beta)))

  parsed <- parseselect(bfile, extract=extract, exclude = exclude,
                        keep=keep, remove=remove,
                        chr=chr)
  if(length(beta) != parsed$p) stop("Vector length of beta does not match number of selected columns in bfile")
  # stopifnot(length(cor) == parsed$p)

  if(!is.null(cluster)) {
    nclusters <- length(cluster)
    if(nclusters > 1) {
      split <- ceiling(seq(1/parsed$p, nclusters, length=parsed$p))
      t <- table(split)
      compute.size <- as.double(min(t)) * parsed$N * length(beta)
      if(compute.size < 1e8) {
        # Too many clusters
        f <- 1e8 / compute.size
        recommended <- min(ceiling(nclusters / f), nclusters - 1)
        return(pgs(bfile, weights, keep=parsed$keep, extract=parsed$extract,
                   cluster=cluster[1:recommended], trace=trace))
      }
      Bfile <- bfile # Define this within the function so that it is copied
                      # to the child processes
      if(trace > 0) cat("Divided into", nclusters, "chunks\n")
      l <- parallel::parLapply(cluster, 1:nclusters, function(i) {
        toextract <- if(!is.null(parsed$extract)) parsed$extract else
          rep(TRUE, parsed$P)
        touse <- split == i
        toextract[toextract] <- touse

        obj <- weights
        obj$beta <- beta[touse]
        obj$pbin <- pbin[touse]

        return(pgs(Bfile, obj, keep=parsed$keep, extract=toextract, trace=trace))
      })
      result <- l[[1]]
      if(nclusters > 1) for(i in 2:nclusters) result <- result + l[[i]]
      return(result)
    }
  }

  if(is.null(parsed$extract)) {
    extract2 <- list(integer(0), integer(0))
  } else {
    extract2 <- lassosum:::selectregion(!parsed$extract)
    extract2[[1]] <- extract2[[1]] - 1
  }

  if(is.null(parsed$keep)) {
    keepbytes <- integer(0)
    keepoffset <- integer(0)
  } else {
    pos <- which(parsed$keep) - 1
    keepbytes <- floor(pos/4)
    keepoffset <- pos %% 4 * 2
  }

  bfile <- paste0(bfile, ".bed")

  return(multiBed4(bfile, parsed$N, parsed$P,
                   beta, pbin, nbin,
                   extract2[[1]], extract2[[2]],
                   keepbytes, keepoffset, trace = trace))

  #' @return A matrix of Polygenic Scores

}
