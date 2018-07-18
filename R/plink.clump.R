plink.clump <- function(bfile, pvals, 
                        r2=0.2, p1=1, p2=1, kb=250, allow.overlap=F, 
                        out=tempfile("clump"), 
                        keep=NULL, remove=NULL, 
                        extract=NULL, exclude=NULL, 
                        chr=NULL, 
                        extension=".clumped",
                        ignore.stdout=trace <= 1, 
                        trace=0, 
                        plink = getOption("lassosum.plink"), 
                        ..., 
                        parsed=NULL,
                        cluster=NULL, 
                        dryrun=FALSE) {
  #' FUnction to call plink --clump using user-supplied p-values
  #' A number of options not implemented here. Need to call plink() directly for those.
  #' @param pvals can be given as a vector (must equal number of SNPs), or
  #' a data.frame with two columns, the first column being SNP id and the second being the p-values.
  #' @param ... options passed to plink.parse()

  if(is.null(parsed)) {
    parsed <- lassosum:::parseselect(bfile=bfile, 
                                     keep=keep, remove=remove, 
                                     extract=extract, exclude=exclude, 
                                     chr=chr)
  }
  
  if(!is.null(cluster)) {
    ### parallel by chromosome ###
    nclusters <- length(cluster)
    if(is.null(parsed$bim)) parsed$bim <- lassosum:::read.table2(parsed$bimfile)
    if(is.null(parsed$extract)) bim <- parsed$bim else bim <- parsed$bim[parsed$extract,]
    chrs <- unique(bim$V1)
    chrs <- sub("^chr", "", chrs, ignore.case = TRUE)
    if(nclusters > 1 && length(chrs) > 1) {
      cmd <- plink.clump(bfile=bfile, pvals, 
                         r2=r2, p1=p1, p2=p2, kb=kb, 
                         allow.overlap=allow.overlap, 
                         out=sub <- "_TO_BE_SUBSTITUTED_", 
                         extension=extension,
                         plink = plink, 
                         ..., 
                         parsed=parsed, 
                         dryrun=TRUE)
      Ig <- ignore.stdout # Force copy for parallel
      Ex <- extension     # Force copy for parallel
      Tr <- trace     # Force copy for parallel
      l <- parallel::parLapplyLB(cl=cluster, chrs, function(i) {
        Cmd <- sub(sub, Out <- paste0(out, i), cmd)
        Cmd <- paste(Cmd, "--chr", i, "2>&1 | sed '/missing/d'")
        if(Tr > 0) cat("Processing chromosome ", i, "\n")
        system(Cmd, ignore.stdout = Ig)
        outtable  <- lassosum:::read.table2(paste0(Out,Ex), header=T, drop=12)
        return(outtable)
      })
      return(do.call(rbind, l))
    }
  }
  
  P <- plink.parse(parsed=parsed, out=out, plink=plink, ..., details=TRUE)

  parsed <- attr(P, "parsed")

  if(is.data.frame(pvals)) {
    stopifnot(ncol(pvals) == 2)
    table <- pvals
  } else {
    stopifnot(is.vector(pvals) & is.numeric(pvals))
    if(is.null(parsed$bim)) {
      parsed$bim <- lassosum:::read.table2(parsed$bimfile)
    }
    stopifnot(parsed$p == length(pvals))
    bim <- parsed$bim
    if(!is.null(parsed$extract)) bim <- bim[parsed$extract, ]
    table <- cbind(bim$V2, pvals)
  }
  colnames(table) <- c("SNP", "P")
  lassosum:::write.table2(table, file=pvalsfile <- tempfile(pattern="pvals"), col.names=T)
  options <- paste("--clump", pvalsfile,
                   "--clump-p1", p1,
                   "--clump-p2", p2,
                   "--clump-r2", r2,
                   "--clump-kb", kb)
  if(allow.overlap) options <- paste(options, "--allow-overlap")
  cmd <- paste(plink, P, options)
  
  #### run plink ####
  if(dryrun) {
    return(cmd)
  } else {
    system(cmd, ignore.stdout = ignore.stdout)
    outtable  <- lassosum:::read.table2(paste0(out,extension), header=T, drop=12)
    return(outtable)
  }

}
