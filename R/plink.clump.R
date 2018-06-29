plink.clump <- function(pvals, r2=0.2, p1=1, p2=1, kb=250, allow.overlap=F, 
                        bfile, out=tempfile("clump"), 
                        keep=NULL, remove=NULL, 
                        extract=NULL, exclude=NULL, 
                        chr=NULL, 
                        extension=".clumped",
                        ignore.stdout=TRUE, 
                        plink = getOption("lassosum.plink"), 
                        ..., 
                        cluster=NULL) {
  #' FUnction to call plink --clump using user-supplied p-values
  #' A number of options not implemented here. Need to call plink() directly for those.
  #' @param pvals can be given as a vector (must equal number of SNPs), or
  #' a data.frame with two columns, the first column being SNP id and the second being the p-values.
  #' @param ... options passed to plink.parse()

  P <- plink.parse(bfile=bfile, out=out,
                   keep=keep, remove=remove, 
                   extract=extract, exclude=exclude, 
                   chr=chr, 
                   plink=plink, 
                   ..., details=TRUE)
  
  parsed <- attr(P, "parsed")
  options <- attr(P, "options")

  if(is.data.frame(pvals)) {
    stopifnot(ncol(pvals) == 2)
    table <- pvals
  } else {
    stopifnot(is.vector(pvals) & is.numeric(pvals))
    if(is.null(parsed$bim)) parsed$bim <- read.table2(options$bfile)
    stopifnot(parsed$p == length(pvals))
    bim <- parsed$bim
    if(!is.null(parsed$extract)) bim <- bim[parsed$extract, ]
    table <- cbind(bim$V2, pvals)
  }
  colnames(table) <- c("SNP", "P")
  write.table2(table, file=pvalsfile <- tempfile(pattern="pvals"), col.names=T)
  options <- paste("--clump", pvalsfile,
                   "--clump-p1", p1,
                   "--clump-p2", p2,
                   "--clump-r2", r2,
                   "--clump-kb", kb)
  if(allow.overlap) options <- paste(options, "--allow-overlap")
  options <- paste(P, options)
  
  #### run plink ####
  outfile <- do.call(plink, options)
  outtable  <- read.table2(paste0(outfile,extension), header=T)
  return(outtable)

}
