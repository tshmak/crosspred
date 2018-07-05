plink.linear <- function(bfile, pheno, out=tempfile("lassosum.out"), 
                         keep=NULL, remove=NULL, 
                         extract=NULL, exclude=NULL, 
                         chr=NULL, 
                         covar=NULL, 
                         plink.cmd=NULL, 
                         extension=NULL,
                         ignore.stdout=trace<=0, 
                         trace=0, 
                         plink = getOption("lassosum.plink"), 
                         plink2 = grepl("plink2$", plink),
                         ..., 
                         parsed=NULL) {
  #' @rdname plink.parse
  #' @title Obtain standardized coefficients from linear regression in PLINK
  #' @param plink.cmd The command passed to plink
  #' @param keep,remove,extract,exclude,chr see parseselect() 
  #' @param extension The extension for the output file generated from the plink command
  #' @param ignore.stdout Ignore stdout. see \code{\link[base]{base::system}}
  #' @param plink plink executive 
  #' @param plink2 Whether we are using plink2
  #' @param ... Other parameters to pass to plink. See details
  #' @export
  #'           
  
  if(is.null(parsed)) {
    parsed <- parseselect(bfile=bfile, keep=keep, remove=remove, 
                          extract=extract, exclude=exclude, 
                          chr=chr)
  }

  P <- plink.parse(parsed=parsed, out=out,
                   pheno=pheno, 
                   covar=covar, 
                   plink=plink, 
                   ...)
  
  #### plink2 ? ####
  if(plink2) {
    if(is.null(plink.cmd)) plink.cmd <- "--linear omit-ref hide-covar"
    if(is.null(extension)) extension <- ".PHENO1.glm.linear"
  } else {
    if(is.null(plink.cmd)) plink.cmd <- "--linear standard-beta hide-covar"
    if(is.null(extension)) extension <- ".assoc.linear"
  }
  

  #### run plink ####
  cmd <- paste(plink, P, plink.cmd)
  system(cmd, ignore.stdout=ignore.stdout)
  
  tab <- lassosum:::read.table2(paste0(out, extension), header=TRUE)
  if(plink2) {
    colnames <- colnames(tab)
    colnames[colnames == "X.CHROM"] <- "CHR"
    colnames[colnames == "ID"] <- "SNP"
    colnames[colnames == "POS"] <- "BP"
    colnames[colnames == "OBS_CT"] <- "NMISS"
    colnames(tab) <- colnames
    tab$BETA <- p2cor(tab$P, n = tab$NMISS, sign=tab$BETA) 
      # I'm using this because standard-beta is not available in plink2 at the moment. 
  }
  attr(tab, "out") <- out

  return(tab)

}

