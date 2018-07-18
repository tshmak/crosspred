cp.meta <- function(cp.plink.linear, 
                    cor, chr=NULL, pos=NULL, snp=NULL, 
                    A1=NULL, A2=NULL, 
                    nonmiss=NULL, weight.pl=NULL, 
                    exclude.ambiguous=TRUE) {
  #' Function to meta-analyse raw data summary statistics with external
  #' summary statistics
  #' @param cp.plink.linear An \code{cp.plink.linear} object
  #' @param cor,chr,pos,snp,A1,A2 see \code{\link{lassosum.pipeline}}
  #' @param nonmiss Number of non-missing observations (can be either a scalar or a vector)
  #' @param weight.pl An arbitrary weight between 0 and 1 given to the plink.linear correlation
  #' @note Either \code{n} or \code{nonmiss} or \code{weight.pl} must be specified
  #' @param exclude.ambiguous Should ambiguous SNPs be excluded? 
  #' @details This function performs a meta-analysis of the correlations 
  #' coefficients as calculated in \code{cp.plink.linear} and some external 
  #' correlations. If correlation coefficients are not available, these
  #' can be converted from p-values using \code{\link{p2cor}}. 
  
  pl <- cp.plink.linear

  ref.bim <- read.table2(paste0(pl$bfile, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
  if(!is.null(pl$extract)) ref.bim <- ref.bim[pl$extract, ]
  
  chrpos <- !is.null(chr) && !is.null(pos)
  if(is.null(snp) && !chrpos) {
    stop("Either snp or chr/pos must be specified.")
  } 
  
  ### Checks ###
  if(is.null(A1) && is.null(A2)) {
    stop("At least one of A1 (alternative allele) or A2 (reference allele) must be specified. Preferably both.")
  } else if(is.null(A1) || is.null(A2)) {
    # message("Matching on 1 allele only.")
  }
  
  stopifnot(!any(is.na(cor)))
  stopifnot(all(cor > -1 & cor < 1))
  
  if(chrpos) {
    stopifnot(length(chr) == length(pos))
    stopifnot(length(chr) == length(cor))
    chr <- as.character(sub("^chr", "", chr, ignore.case = T))
  }
  
  if(!is.null(snp)) {
    stopifnot(length(snp) == length(cor))
  }
  stopifnot(is.null(A1) || length(A1) == length(cor))
  stopifnot(is.null(A2) || length(A2) == length(cor))

  ### ss ###
  ss <- list(chr=chr, pos=pos, A1=A1, A2=A2, snp=snp, cor=cor)
  ss[sapply(ss, is.null)] <- NULL
  ss <- as.data.frame(ss)
  
  ### matchpos ###
  m <- matchpos(ss, ref.bim, auto.detect.ref = F, 
                    ref.chr = "V1", ref.snp="V2", 
                    ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                    rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                    silent=T)
  
  stopifnot(!is.null(nonmiss))
  if(length(nonmiss) == 1) {
    nonmiss <- rep(nonmiss, length(m$order))
  } else {
    stopifnot(length(nonmiss) == length(cor))
    nonmiss <- nonmiss[m$order]
  }
  
  ss <- ss[m$order,]
  ss$cor <- ss$cor * m$rev
  
  ### cor ###
  
  ext <- m$ref.extract
  
  for(i in 1:length(pl$cor)) {
    if(!is.null(weight.pl)) {
      stopifnot(length(weight.pl) == 1 && weight.pl >= 0 && weight.pl <= 1)
      weight1 <- weight.pl
      weight2 <- (1-weight.pl)
      N <- 1/(weight1^2 / pl$nonmiss[[i]][ext] + weight2^2 / nonmiss) # Supposing var(cor) = 1/N
    } else {
      N <- pl$nonmiss[[i]][ext]  + nonmiss
      weight1 <- pl$nonmiss[[i]][ext] /  N
      weight2 <- nonmiss / N
    }
    pl$cor[[i]] <- pl$cor[[i]][ext] * weight1 + 
      ss$cor * weight2
    pl$nonmiss[[i]] <- N
  }
  
  for(vars in c("chr", "pos", "A1", "A2", "snp")) {
    pl[[vars]] <- pl[[vars]][ext]
  }

  pl$extract[pl$extract] <- ext
  # We don't update pl$n, because this has to agree with pl$keep
  # pl$n <- pl$n + n
  pl$p <- sum(ext)
  return(pl)
  #' @return A \code{\link{cp.plink.linear}} object 
  

}