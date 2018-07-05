clear()
Tim.load(mysimtools)
# Tim.load(Rplink)
Tim.load(lassosum, export_all=F)
load_all(attachroot("~/WORK/Projects/validation/crosspred/."))

filename <- Rfilename("test", seed=1234)
setwd0(attachroot("~/WORK/Projects/validation/crosspred/tests/"))

bfile <- paste0(system.file("data", package="lassosum"), "/testsample")
# dim.bfile(bfile) 200 * 1800
pheno <- rnorm(nrow.bfile(bfile))

# options(lassosum.plink ="/home/tshmak/software/plink/v1.90b5.2/plink")
test1 <- cp.plink.linear(bfile = bfile, nfolds=2, pheno = pheno, 
                        ignore.stdout=T)


options(lassosum.plink="/home/tshmak/software/plink/v2.00a2LM_AVX2_20180605/plink2")
test2 <- cp.plink.linear(bfile = bfile, nfolds=2, pheno = pheno, fold=test1$fold, 
                        ignore.stdout=T)
