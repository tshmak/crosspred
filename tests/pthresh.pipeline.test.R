clear()
Tim.load(mysimtools)
# Tim.load(Rplink)
Tim.load(lassosum, export_all=F)
load_all(attachroot("~/WORK/Projects/validation/crosspred/."))

filename <- Rfilename("test", seed=1234)
setwd0(attachroot("~/WORK/Projects/validation/crosspred/tests/"))

bfile <- paste0(system.file("data", package="lassosum"), "/testsample")
# dim.bfile(bfile) 200 * 1800
pval <- runif(ncol.bfile(bfile))
cor <- p2cor(pval,n = nrow.bfile(bfile))
bim <- read.table.tim(paste0(bfile, '.bim'))

# options(lassosum.plink ="/home/tshmak/software/plink/v1.90b5.2/plink")
test <- pthresh.pipeline(test.bfile = bfile, beta = cor, pval=pval, snp=bim$V2,
                         A1=bim$V5, trace=2)

test <- pthresh.pipeline(test.bfile = bfile, beta = cor, pval=pval, snp=bim$V2,
                         A1=bim$V5, trace=2, clump=T)
refbfile <- paste0(system.file("data", package="lassosum"), "/refpanel")
test <- pthresh.pipeline(test.bfile = bfile, beta = cor, pval=pval, snp=bim$V2,
                         A1=bim$V5, trace=2, clump=T, ref.bfile=refbfile)
