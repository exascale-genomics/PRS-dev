args <- commandArgs(TRUE)
cov_file <- args[1]
eigenvec_file <- args[2]
covariate_file <- args[3]

covariate <- read.table(cov_file, header=T)
pcs <- read.table(eigenvec_file, header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,covariate_file, quote=F, row.names=F)