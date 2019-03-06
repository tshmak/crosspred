Cross prediction for PGS in large cohorts
=======================
  
Cross prediction is a method for calculating polygenic scores in large cohorts without the use of summary statistics. For details, refer to this [paper](https://www.biorxiv.org/content/early/2018/01/23/252270). 

_Note: In the earlier version of the paper (prior to July, 2018), two methods of cross-prediction are proposed: Method 1 and Method 2. In the newer version, Method 1 is referred to as Stack and validate, while Method 2 is referred to as Split-validation. I recommend using Split-validation (Method 2) where possible._

# Installation

If **lassosum** is not yet installed, refer to the instruction [here](https://github.com/tshmak/lassosum#installation) for installation. **lassosum** v0.4.0 or above is required. 

Install `crosspred` using `devtools` (Note: Windows users need to have installed [Rtools](https://cran.r-project.org/bin/windows/Rtools/).): 
```r
# install.packages("devtools") # If devtools not yet installed. 
devtools::install_github("tshmak/crosspred")
```

Some functions in **crosspred** use [Plink](https://www.cog-genomics.org/plink2/) to calculate summary statistics. Before these functions can be used, the link to the plink executable needs to be specified by: 
```r
options(lassosum.plink='/path/to/plink')
```

# Tutorial

We assume that we have genotype in PLINK 1 [format](https://www.cog-genomics.org/plink/1.9/input#bed). For example, let's say our data files are: `mydata.bed`, `mydata.bim`, and `mydata.fam`. 

We illustrate `crosspred` using the toy example PLINK dataset that comes with the `lassosum` package: 
```r
library(lassosum) 
setwd(system.file("data", package="lassosum"))
random.pheno <- rnorm(nrow.bfile("testsample")) # A random phenotype
pl <- cp.plink.linear("testsample", pheno=random.pheno, nfold=2) 
  # This generates the cp.plink.linear object to be used for cross-prediction.
# plot(pl$cor[[1]], pl$cor[[2]]) # The correlation statistics for the two folds can be obtained thus. 
```
The default is 5-fold cross-prediction. The folds are randomly assigned to the samples. Use `nfolds` to specify the number of folds needed, or `fold` to allocate the fold yourself. Type `help(cp.plink.linear)` for more details. 

`sumstats` can then be fed into `cp.lassosum`. 
```r 
ldblocks <- data.table::fread("Berisa.EUR.hg19.bed")
cp <- cp.lassosum(pl, LDblocks = ldblocks)
# plot(cp$pheno, cp$best.pgs) # This is the best PGS using Method 1 (Stack and validate)
# plot(cp$pheno, cp$best.pgs.m2) # This is the best PGS using Method 2 (Split-validation)
```
We recommend you use one of the LD blocks given in `lassosum`, if you do not have your own LD blocks defined. These LD blocks are given by the paper [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium), and are based on the 1000 Genome data. Replace EUR with ASN or AFR for Asian or African LD regions, respectively. hg38 coordinates are also available by [liftOver](https://genome.sph.umich.edu/wiki/LiftOver). Simply replace `hg19` with `hg38` above. 

By default, if the sample size is > 5000, `cp.lassosum` uses a random subset of 5000 as the reference panel. This is to reduce the computation burden. Increase or decrease this using the `max.ref.bfile.n` option. Alternatively, specify the exact sample to use using the `keep.ref` option. 

For further options, please see the manual
```r
help(cp.lassosum)
```
or email me at <timmak@yahoo.com>. 

### Multi-threading 
Both `cp.plink.linear` and `cp.lassosum` can be run much quicker by multi-threading. For `cp.plink.linear`, multi-threading is performed by PLINK. For this to work, we need the experimental [PLINK 2](https://www.cog-genomics.org/plink/2.0/). For `cp.lassosum`, we can specify the `cluster` option. 

### Running cross-prediction using separated .bed files
In large datasets, data is often stored across chromosomes in separate `.bed` files. To run cross prediction across different chromosomes, simply run `cp.plink.linear` separately for the different `.bed` files, then run `cp.lassosum` separately for each output from `cp.plink.linear`, specifying `list.of.lpipe.output=TRUE`. This will generate a list of `lassosum.pipeline` objects for each fold. However, **remember to set the random number seed to the same number before you run `cp.plink.linear` so that the folds are defined consistently, and also before `cp.lassosum` if your sample size is > `max.ref.bfile.n`, so that the same reference sample is chosen!** Then, use `organize.by.fold` to reorganize the all the outputs for the different chromosomes. Finally, run `cp.lassosum` with the `list.of.lpipe.input` option. Below is an example for two chromosomes. 
```r
set.seed(42)
sumstats1 <- cp.plink.linear("chr1")
set.seed(42)
sumstats2 <- cp.plink.linear("chr2")
set.seed(1000)
lp1 <- cp.lassosum(sumstats1, LDblocks=ld, list.of.lpipe.output=TRUE)
set.seed(1000)
lp2 <- cp.lassosum(sumstats2, LDblocks=ld, list.of.lpipe.output=TRUE)
lp <- organize.by.fold(list(lp1, lp2))
cp <- cp.lassosum(list.of.lpipe.input=lp)
```
### Running cross-validation using results from cross-prediction
Method 1 in cross-prediction essentially identifies the best `lambda` and `s` for use with **lassosum**. We can then apply this to the entire dataset. This is the cross-validation procedure in our paper, which can be achieved by running: 
```r
cv <- cp.cv(sumstats, cp)
```
`sumstats` can be a `list` of `cp.plink.linear` objects if they were obtained separately for different chromosomes. Note that the order of the list should match that given to `cp.lassosum` earlier. Note that `cp.cv` does not calculate summary statistics afresh using the entire data, but simply averages those obtained across the different folds in `cp.lassosum`. This makes it very fast. 

### Incorporating external summary statistics into cross-prediction 
Suppose external summary statistics are loaded into a `data.frame` called `dat` with the columns `chr`, `pos`, `alt`, `rsid`, `beta`, `p` for chromosomes, position, alternative allele, rsID, beta coefficients, and p-values respectively. We first need to convert the p-values to correlation coefficients using the `p2cor` command: 
```r
dat$cor <- p2cor(dat$p, n = 60000, sign=dat$beta)
```
where `n = 60000` is the sample size. Assuming `cp.pl` is an object from `cp.plink.linear`, we can then run 
```r
merged <- cp.meta(cp.pl, chr=chr, pos=pos, snp=rsid, A1=alt, cor=cor, n=60000)
```
to merge it with `cp.pl`. We can then proceed with `merged` as if it were a `cp.plink.linear` object in running, e.g., 
```r 
cp <- cp.lassosum(merged, LDblocks=ld)
```
Note that not all of the variables `chr`, `pos`, `snp`, `A1`, `A2` need to be specified. At least one of `A1` and `A2` must be specified. Otherwise specify as many of `chr`, `pos`, and `snp` as you need. 

### Using clumping and thresholding instead of lassosum
There is an equivalent function `cp.pthresh` for running cross-prediction using clumping and p-value thresholding instead of lassosum. Refer to the documentation for details. 

Please email me <timmak@yahoo.com> for any bug reports, comments, or suggestions. Thanks for using **crosspred**!



