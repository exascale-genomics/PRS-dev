# HPC Pipeline For Calculating Polygenic Risk Scores (PRS) in Cancer

## Introduction
Polygenic risk scores (PRSs) are weighted genetic scores that are important in assessing and predicting risk for various clinical disorders. Multiple methods are being used to generate PRSs from summary genome-wide association study (GWAS) statistics. These methods are developed using multiple, complementary techniques like pruning and thresholding, bayesian approaches, and penalized regression and have relative strengths and weaknesses. Preparing the input files, configuring the analysis programs, and performing the post-analysis quality control and plots is a cumbersome process. We created a high performance, parallel analysis pipeline that enables rapid generation of scores using DOE computational resources. 

## Methods
We developed a high performance, robust workflow that uses the input summary GWAS data and performs standard quality control and generates Polygenic Risk Scores with multiple software tools such as PLINK, PRSice-2, LDpred-2, lassosum, PRS-CSx and SBayesR. This workflow was configured to leverage the high-performance computational resources available in the Oak Ridge National Laboratory KDI cluster and can be easily configured to run on other HPC resources and on cloud resources.

### PRS-dev

Use nextflow to build a PRS pipeline based on: 
```
https://choishingwan.github.io/PRS-Tutorial/ # PRS with PLINK, PRSice-2, LDpred-2, lassosum
https://github.com/getian107/PRScs           # PRS-CS and PRS-CSx
https://github.com/getian107/PRScsx
https://cnsgenomics.com/software/gctb        # SBayesR
```

### Requirements:
install nextflow:
```
require java
https://www.nextflow.io/docs/latest/getstarted.html
```

### R:
```
https://cran.r-project.org/bin/linux/ubuntu/
run as root:
apt update -qq
apt install --no-install-recommends software-properties-common dirmngr
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt install --no-install-recommends r-base
```

### plink:
```
https://www.cog-genomics.org/plink
export PATH="/home/ubuntu/plink:$PATH"
```

### PRSice:
```
https://choishingwan.github.io/PRS-Tutorial/prsice/
put it in /home/ubuntu/PRSice or somewhere accessible to the program
```

### LDpred-2:
```
run as root
apt install r-base r-base-dev
apt install make build-essential libcurl4-openssl-dev libssl-dev libxml2-dev
R
install.packages("bigsnpr")
```

### lassosum:
```
run as root
R
install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")
```

### PRS-CS and PRS-CSx:
```
git clone https://github.com/getian107/PRScs.git
Download the LD reference panels and extract files
sudo apt install python3-pip
sudo pip3 install scipy h5py

git clone https://github.com/getian107/PRScsx.git
Download the SNP information file and put it in the same folder containing the reference panels
```

### SBayesR:
```
https://cnsgenomics.com/software/gctb/#Tutorial
export PATH="/scratch/SBayesR/gctb_2.03beta_Linux:$PATH"
```

#### Base data applied QC:
File transfer
Standard GWAS QC
Duplicate SNPs
Ambiguous SNPs

#### target data applied QC:
File transfer
Standard GWAS QC: removing SNPs with low genotyping rate, low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing individuals with low genotyping rate; perform pruning to remove highly correlated SNPs; compute heterozygosity rates; remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean
Mismatching SNPs
Sex chromosomes
Relatedness


#### base data missing QC: 
Heritability check
Effect allele
Genome build
Mismatching SNPs
Sex chromosomes
Sample overlap
Relatedness

#### target data missing QC:
Sample size
Genome build
Sample overlap

