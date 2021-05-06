args <- commandArgs(TRUE)
file_height <- args[1]
file_cov <- args[2]
file_eigenvec <- args[3]
file_bed <- args[4]
file_rds <- args[5]
file_base <- args[6]

# prepare workspace and load bigsnpr
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

# read in phenotype and covariates
library(data.table)
library(magrittr)
phenotype <- fread(file_height)
covariate <- fread(file_cov)
pcs <- fread(file_eigenvec)
# rename columns
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
# generate required table
pheno <- merge(phenotype, covariate) %>%
    merge(., pcs)

# load HapMap3 SNPs
info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))

# Read in the summary statistic file
sumstats <- bigreadr::fread2(file_base) 
# LDpred 2 require the header to follow the exact naming
names(sumstats) <-
    c("chr",
    "pos",
    "rsid",
    "a1",
    "a0",
    "n_eff",
    "beta_se",
    "p",
    "OR",
    "INFO",
    "MAF")
# Transform the OR into log(OR)
sumstats$beta <- log(sumstats$OR)
# Filter out hapmap SNPs
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

# Calculate the LD matrix (Genome Wide bed file option) (there is also Chromosome separated bed files option)
# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# preprocess the bed file (only need to do once for each data set)
snp_readBed(file_bed)
# now attach the genotype object
obj.bigSNP <- snp_attach(file_rds)
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map)
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
# calculate LD
for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
            genotype,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
        c("family.ID", "sample.ID"),
        c("FID", "IID"))

# Perform LD score regression
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

# Calculate the null R2
# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
y <- pheno[fam.order, on = c("FID", "IID")]
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    paste0("Height~Sex+", .) %>%
    as.formula %>%
    lm(., data = y) %>%
    summary
null.r2 <- null.model$r.squared

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

# Obtain model PRS
if(is.null(obj.bigSNP)){
    obj.bigSNP <- snp_attach(file_rds)
}
genotype <- obj.bigSNP$genotypes
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
pred_inf <- big_prodVec(    genotype,
                            beta_inf,
                            ind.row = ind.test,
                            ind.col = info_snp$`_NUM_ID_`)

# Get the final performance of the LDpred models
reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    paste0("Height~PRS+Sex+", .) %>%
    as.formula
reg.dat <- y
reg.dat$PRS <- pred_inf
inf.model <- lm(reg.formula, dat=reg.dat) %>%
    summary
(result <- data.table(
    infinitesimal = inf.model$r.squared - null.r2,
    null = null.r2
))
