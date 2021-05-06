args <- commandArgs(TRUE)
file_height <- args[1]
file_cov <- args[2]
file_eigenvec <- args[3]
file_base <- args[4]
val_prefix_qc <- args[5]
val_prefix <- args[6]
file_ld <- args[7]

library(lassosum)
# Prefer to work with data.table as it speeds up file reading
library(data.table)
library(methods)
library(magrittr)
# For multi-threading, you can use the parallel package and 
# invoke cl which is then passed to lassosum.pipeline
library(parallel)
# This will invoke 2 threads. 
cl <- makeCluster(2)

sum.stat <- file_base
bfile <- val_prefix_qc
# Read in and process the covariates
covariate <- fread(file_cov)
pcs <- fread(file_eigenvec) %>%
    setnames(., colnames(.), c("FID","IID", paste0("PC",1:6)))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs)

# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <- file_ld
# output prefix
prefix <- val_prefix
# Read in the target phenotype file
target.pheno <- fread(file_height)[,c("FID", "IID", "Height")]
# Read in the summary statistics
ss <- fread(sum.stat)
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
        n = ss$N,
        sign = log(ss$OR)
        )
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]


# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file, 
    cluster=cl
)
# Store the R2 results
target.res <- validate(out, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov))
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2