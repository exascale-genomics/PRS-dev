args <- commandArgs(TRUE)
valid_sample <- args[1]
sex_check <- args[2]
output_valid <- args[3]

# Read in file
valid <- read.table(valid_sample, header=T)
dat <- read.table(sex_check, header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], output_valid, row.names=F, col.names=F, sep="\t", quote=F) 