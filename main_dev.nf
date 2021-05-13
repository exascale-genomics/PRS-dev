// TODO: add Genome build QC
// TODO: PRS-CS and PRS-CSx: trim base date columes; redo the inputs part; more work for the output; PRS-CSx code is broken in testing, need investigation, so it is not added here
// TODO: SBayesR: input prepararion steps; look into outputs

// parameters
params.base_data = 'Height.gwas.txt.gz'
params.target_data = '/home/ubuntu/inputs_target_test'
params.outdir = '/home/ubuntu/outputs'
params.PRSice = '/home/ubuntu/PRSice'
params.workingDir = System.getProperty("user.dir")

// settings
//// PRSice
params.PRSice = '/home/ubuntu/PRSice'
//// lassosum
lassosum_prefix_to_ld_mapping = ["EUR":"EUR.hg19", "TMP":"EUR.hg19"]
//// PRS-CS
params.prs_cs_script_dir = '/scratch/PRScs'
params.prs_cs_ref_dir = '/scratch/ldblk/ldblk_1kg_eur' 
params.prs_cs_bim_prefix = '/scratch/PRScs/test_data/test'
params.prs_cs_sst_file = '/scratch/PRScs/test_data/sumstats.txt'
params.prs_cs_n_gwas = '200000'
params.prs_cs_chrom = '22'
//// SBayesR
params.sbayesr_ldm = '/scratch/SBayesR/gctb_2.0_tutorial/ldm/sparse/chr22/1000G_eur_chr22.ldm.sparse'
params.sbayesr_gwas_summary = '/scratch/SBayesR/gctb_2.0_tutorial/ma/sim_1.ma'


// helper: checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

bin = "${params.workingDir}/bin"
base_data = file(params.base_data)
Channel.fromFilePairs("$params.target_data/*.{bed,bim,cov,fam,height}", size: 5, flat: true)\
                     { file -> file.baseName }\
                     .ifEmpty { error "No matching target files" }\
                     .map { a -> [a[1].baseName, checker(a[1]), checker(a[2]), checker(a[3]), checker(a[4]), checker(a[5])] }\
                     .multiMap { it -> copy1: copy2: copy3: copy4: copy5: copy6: copy7: copy8: copy9: copy10: it }
                     .set { target_data }
/*
///////// base data QC /////////
// base data QC: File transfer
process base_qc_file_transfer {
    
    publishDir "${params.outdir}/md5sum_results", mode: 'copy'

    input:
    file base_data

    output:
    file("base_data_md5sum") into base_qc_file_transfer_output

    """
    md5sum $base_data > base_data_md5sum
    """
}

// base data QC: Standard GWAS QC
process base_qc_standard_gwas_qc {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    file base_data

    output:
    file('base_qc_standard_gwas_qc_output.gz') into base_qc_standard_gwas_qc_output

    shell:
    '''
    gunzip -c !{base_data} |\
    awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
    gzip  > base_qc_standard_gwas_qc_output.gz
    '''
}

// base data QC: Duplicate SNPs
process base_qc_duplicate_snps_get_dup {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    file base_qc_standard_gwas_qc_output

    output:
    file('base_qc_duplicate_snps.snp') into base_qc_duplicate_snps_get_dup_output

    shell:
    '''
    gunzip -c !{base_qc_standard_gwas_qc_output} |\
    awk '{ print $3}' |\
    sort |\
    uniq -d > base_qc_duplicate_snps.snp
    '''
}

process base_qc_duplicate_snps_remove_dup {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    file base_qc_standard_gwas_qc_output
    file base_qc_duplicate_snps_get_dup_output

    output:
    file('base_qc_remove_dup_output.nodup.gz') into base_qc_duplicate_snps_remove_dup_output

    shell:
    '''
    gunzip -c !{base_qc_standard_gwas_qc_output}  |\
    grep -vf !{base_qc_duplicate_snps_get_dup_output} |\
    gzip - > base_qc_remove_dup_output.nodup.gz
    '''
}

// base data QC: Ambiguous SNPs
process base_qc_ambiguous_snps {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    file base_qc_duplicate_snps_remove_dup_output

    output:
    file('base_qc_output.QC.gz') into base_qc_final_output

    shell:
    '''
    gunzip -c !{base_qc_duplicate_snps_remove_dup_output} |\
    awk '!( ($4=="A" && $5=="T") || \
            ($4=="T" && $5=="A") || \
            ($4=="G" && $5=="C") || \
            ($4=="C" && $5=="G")) {print}' |\
        gzip > base_qc_output.QC.gz
    '''
}

///////// target data QC /////////
// target data QC: File transfer
process target_qc_file_transfer {

    publishDir "${params.outdir}/md5sum_results", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy1

    output:
    file("${prefix}_target_data_md5sum") into target_qc_file_transfer_output

    """
    md5sum $bed >> ${prefix}_target_data_md5sum
    md5sum $bim >> ${prefix}_target_data_md5sum
    md5sum $cov >> ${prefix}_target_data_md5sum
    md5sum $fam >> ${prefix}_target_data_md5sum
    md5sum $height >> ${prefix}_target_data_md5sum
    """
}

// target data QC: Standard GWAS QC
//// removing SNPs with low genotyping rate, low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing individuals with low genotyping rate
//// perform pruning to remove highly correlated SNPs
//// compute heterozygosity rates
process target_qc_standard_gwas_qc_1 {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy2

    output:
    tuple val("${prefix}"), file("${prefix_qc}.irem"), file("${prefix_qc}.hh"), file("${prefix_qc}.snplist"), \
    file("${prefix_qc}.fam"), file("${prefix_qc}.prune.in"), file("${prefix_qc}.prune.out"), \
    file("${prefix_qc}.het") into target_qc_standard_gwas_qc_1_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    plink --bfile ${prefix} --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out ${prefix_qc}
    mkdir -p ${params.outdir}/logs
    if test -f "${prefix_qc}.log";
        then mv ${prefix_qc}.log ${params.outdir}/logs/${prefix_qc}.log.1
    fi

    plink --bfile ${prefix} --keep ${prefix_qc}.fam --extract ${prefix_qc}.snplist --indep-pairwise 200 50 0.25 --out ${prefix_qc}
    if test -f "${prefix_qc}.log";
        then mv ${prefix_qc}.log ${params.outdir}/logs/${prefix_qc}.log.2
    fi

    plink --bfile ${prefix} --extract ${prefix_qc}.prune.in --keep ${prefix_qc}.fam --het --out ${prefix_qc}
    if test -f "${prefix_qc}.log";
        then mv ${prefix_qc}.log ${params.outdir}/logs/${prefix_qc}.log.3
    fi
    """
}

target_qc_standard_gwas_qc_1_output.into { target_qc_standard_gwas_qc_1_output_copy1; target_qc_standard_gwas_qc_1_output_copy2; target_qc_standard_gwas_qc_1_output_copy3; target_qc_standard_gwas_qc_1_output_copy4; target_qc_standard_gwas_qc_1_output_copy5 }

//// remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean
process target_qc_standard_gwas_qc_2 {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het) from target_qc_standard_gwas_qc_1_output_copy1

    output:
    file("${prefix}.valid.sample") into target_qc_standard_gwas_qc_2_output

    script:
    """
    #!/usr/bin/env Rscript

    dat <- read.table("${het}", header=T) # Read in the EUR.het file, specify it has header
    m <- mean(dat\$F) # Calculate the mean  
    s <- sd(dat\$F) # Calculate the SD
    valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
    write.table(valid[,c(1,2)], "${prefix}.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
    """
}

// target data QC: Mismatching SNPs
process target_qc_mismatching_snps {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy3
    tuple _tmp, file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het) from target_qc_standard_gwas_qc_1_output_copy2
    file base_qc_final from base_qc_final_output

    output:
    tuple file("${prefix}.a1"), file("${prefix}.mismatch") into target_qc_mismatching_snps_output

    script:
    """
    R --file=${bin}/target_qc_mismatching_snps.R --args ${bim} ${snplist} ${base_qc_final} ${prefix}.a1 ${prefix}.mismatch
    """  
}

// target data QC: Sex chromosomes
process target_qc_sex_check {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy4
    tuple _tmp, file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het) from target_qc_standard_gwas_qc_1_output_copy3
    file valid_sample from target_qc_standard_gwas_qc_2_output

    output:
    tuple file("${prefix_qc}.sexcheck"), file("${prefix_qc}.valid"), file("${prefix_qc}.hh") into target_qc_sex_check_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    plink --bfile ${prefix} --extract ${prune_in} --keep ${valid_sample} --check-sex --out ${prefix_qc}
    if test -f "${prefix_qc}.log";
        then mv ${prefix_qc}.log ${params.outdir}/logs/${prefix_qc}.log.4
    fi

    R --file=${bin}/target_qc_sex_check.R --args ${valid_sample} ${prefix_qc}.sexcheck ${prefix_qc}.valid
    """  
}

// target data QC: Relatedness
process target_qc_relatedness {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy5
    tuple _tmp, file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het) from target_qc_standard_gwas_qc_1_output_copy4
    tuple file(sexcheck), file(valid), file("tmp_hh") from target_qc_sex_check_output

    output:
    file("${prefix_qc}.rel.id") into target_qc_relatedness_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    plink --bfile ${prefix} --extract ${prune_in} --keep ${valid} --rel-cutoff 0.125 --out ${prefix_qc}
    if test -f "${prefix_qc}.log";
        then mv ${prefix_qc}.log ${params.outdir}/logs/${prefix_qc}.log.5
    fi
    """  
}

// target data QC: Generate final QC'ed target data file
process target_qc_final {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy6
    tuple _tmp, file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het) from target_qc_standard_gwas_qc_1_output_copy5
    tuple file(a1), file(mismatch) from target_qc_mismatching_snps_output
    file rel_id from target_qc_relatedness_output

    output:
    tuple val("${prefix}"), file("${prefix_qc}.fam"), file("${prefix_qc}.bed"), file("${prefix_qc}.bim") into target_qc_final_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    plink --bfile ${prefix} --make-bed --keep ${rel_id} --out ${prefix_qc} --extract ${snplist} --exclude ${mismatch} --a1-allele ${a1}
    if test -f "${prefix_qc}.log";
        then mv ${prefix_qc}.log ${params.outdir}/logs/${prefix_qc}.log.6
    fi
    """  
}

target_qc_final_output.into { target_qc_final_output_copy1; target_qc_final_output_copy2; target_qc_final_output_copy3; target_qc_final_output_copy4; target_qc_final_output_copy5; target_qc_final_output_copy6 }

///////// PRS using plink /////////
// Update Effect Size
process plink_update_effect_size {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    file base_qc_final from base_qc_final_output

    output:
    file("${output_file}") into plink_update_effect_size_output

    script:
    base = base_qc_final.baseName
    output_file = "${base}.Transformed"
    """
    #!/usr/bin/env Rscript

    dat <- read.table(gzfile("${base_qc_final}"), header=T)
    dat\$BETA <- log(dat\$OR)
    write.table(dat, "${output_file}", quote=F, row.names=F)
    """  
}

// Clumping
process plink_clumping {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(fam), file(bed), file(bim) from target_qc_final_output_copy1
    file transformed from plink_update_effect_size_output


    output:
    tuple file("${prefix}.clumped"), file("${prefix}.valid.snp") into plink_clumping_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    plink --bfile ${prefix_qc} --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${transformed} --clump-snp-field SNP --clump-field P --out ${prefix}
    if test -f "${prefix}.log";
        then mv ${prefix}.log ${params.outdir}/logs/${prefix}.log.7
    fi
    awk 'NR!=1{print \$3}' ${prefix}.clumped >  ${prefix}.valid.snp
    """
}

// Generate PRS
process plink_generate_prs {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(fam), file(bed), file(bim) from target_qc_final_output_copy2
    tuple file(clumped), file(valid_snp) from plink_clumping_output
    file transformed from plink_update_effect_size_output

    output:
    tuple file("${prefix}.0.5.profile"), file("${prefix}.0.4.profile"), file("${prefix}.0.3.profile"), file("${prefix}.0.2.profile"), file("${prefix}.0.1.profile"), file("${prefix}.0.05.profile"), file("${prefix}.0.001.profile") into plink_generate_prs_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    awk '{print \$3,\$8}' ${transformed} > SNP.pvalue
    echo "0.001 0 0.001" > range_list 
    echo "0.05 0 0.05" >> range_list
    echo "0.1 0 0.1" >> range_list
    echo "0.2 0 0.2" >> range_list
    echo "0.3 0 0.3" >> range_list
    echo "0.4 0 0.4" >> range_list
    echo "0.5 0 0.5" >> range_list
    plink --bfile ${prefix_qc} --score ${transformed} 3 4 12 header --q-score-range range_list SNP.pvalue --extract ${valid_snp} --out ${prefix}
    if test -f "${prefix}.log";
        then mv ${prefix}.log ${params.outdir}/logs/${prefix}.log.8
    fi
    """
}

// Accounting for Population Stratification
process plink_accounting_for_population_stratification {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(fam), file(bed), file(bim) from target_qc_final_output_copy3

    output:
    tuple file("${prefix}.prune.in"), file("${prefix}.prune.out"), file("${prefix}.eigenvec"), file("${prefix}.eigenval") into plink_accounting_for_population_stratification_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    # First, we need to perform prunning
    plink --bfile ${prefix_qc} --indep-pairwise 200 50 0.25 --out ${prefix}
    if test -f "${prefix}.log";
        then mv ${prefix}.log ${params.outdir}/logs/${prefix}.log.9
    fi

    # Then we calculate the first 6 PCs
    plink --bfile ${prefix_qc} --extract ${prefix}.prune.in --pca 6 --out ${prefix}
    if test -f "${prefix}.log";
        then mv ${prefix}.log ${params.outdir}/logs/${prefix}.log.10
    fi
    """
}

plink_accounting_for_population_stratification_output.into { plink_accounting_for_population_stratification_output_copy1; plink_accounting_for_population_stratification_output_copy2;
plink_accounting_for_population_stratification_output_copy3;
plink_accounting_for_population_stratification_output_copy4 }

// Finding the "best-fit" PRS
process plink_finding_the_best_fit_prs {

    publishDir "${params.outdir}/PRS_plink_results", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy7
    tuple file(prune_in), file(prune_out), file(eigenvec), file(eigenval) from plink_accounting_for_population_stratification_output_copy1
    tuple file(profile1), file(profile2), file(profile3), file(profile4), file(profile5), file(profile6), file(profile7) from plink_generate_prs_output

    output:
    file("${prefix}_plink_result") into plink_finding_the_best_fit_prs_output

    script:
    """
    R --file=${bin}/regression_between_prs.R --args ${height} ${eigenvec} ${cov} ${prefix} > ${prefix}_plink_result
    """  
}

///////// PRS using PRSice-2 /////////
process PRSice_analysis {

    publishDir "${params.outdir}/intermediate_files", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy8
    tuple file(prune_in), file(prune_out), file(eigenvec), file(eigenval) from plink_accounting_for_population_stratification_output_copy2
    file base_qc_final from base_qc_final_output
    tuple _tmp, file(fam1), file(bed1), file(bim1) from target_qc_final_output_copy4

    output:
    file("${prefix}.covariate") into prsice_analysis_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    R --file=${bin}/create_covariate.R --args ${cov} ${eigenvec} ${prefix}.covariate

    Rscript ${params.PRSice}/PRSice.R \
    --prsice ${params.PRSice}/PRSice_linux \
    --base ${base_qc_final} \
    --target ${prefix_qc} \
    --binary-target F \
    --pheno ${height} \
    --cov ${prefix}.covariate \
    --base-maf MAF:0.01 \
    --base-info INFO:0.8 \
    --stat OR \
    --or \
    --out ${prefix}

    mkdir -p ${params.outdir}/PRSice_results
    mv ${prefix}.summary ${params.outdir}/PRSice_results/
    mv ${prefix}.prsice ${params.outdir}/PRSice_results/
    mv ${prefix}.log ${params.outdir}/PRSice_results/
    mv ${prefix}.best ${params.outdir}/PRSice_results/
    mv ${prefix}*.png ${params.outdir}/PRSice_results/
    """  
}

///////// PRS using LDpred-2 /////////
process LDpred_2_analysis {

    publishDir "${params.outdir}/LDpred_2_results", mode: 'copy'

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy9
    tuple file(prune_in), file(prune_out), file(eigenvec), file(eigenval) from plink_accounting_for_population_stratification_output_copy3
    file base_qc_final from base_qc_final_output
    tuple _tmp, file(fam1), file(bed1), file(bim1) from target_qc_final_output_copy5

    output:
    file("${prefix}_LDpred_2_result") into LDpred_2_analysis_output

    script:
    prefix_qc = "${prefix}.QC"
    """
    R --file=${bin}/LDpred_2_analysis.R --args ${height} ${cov} ${eigenvec} ${bed1} ${prefix_qc}.rds ${base_qc_final} > ${prefix}_LDpred_2_result
    """  
}

///////// PRS using lassosum /////////
process lassosum_analysis {

    input:
    tuple prefix, file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy10
    tuple file(prune_in), file(prune_out), file(eigenvec), file(eigenval) from plink_accounting_for_population_stratification_output_copy4
    file base_qc_final from base_qc_final_output
    tuple _tmp, file(fam1), file(bed1), file(bim1) from target_qc_final_output_copy6

    script:
    prefix_qc = "${prefix}.QC"
    ld_file = lassosum_prefix_to_ld_mapping["${prefix}"]
    """
    R --file=${bin}/lassosum_analysis.R --args ${height} ${cov} ${eigenvec} ${base_qc_final} ${prefix_qc} ${prefix} ${ld_file}

    mkdir -p ${params.outdir}/lassosum_results
    if test -f "Rplots.pdf";
        then mv Rplots.pdf ${params.outdir}/lassosum_results/${prefix}.Rplots.pdf
    fi
    """  
}
*/
///////// PRS-CS /////////
// NEED UPDATE
process PRS_CS_step {

    script:
    """
    mkdir -p ${params.outdir}/PRS_CS_results

    python3 ${params.prs_cs_script_dir}/PRScs.py --ref_dir=${params.prs_cs_ref_dir} --bim_prefix=${params.prs_cs_bim_prefix} --sst_file=${params.prs_cs_sst_file} --n_gwas=${params.prs_cs_n_gwas} --chrom=${params.prs_cs_chrom} --phi=1e-2 --out_dir=${params.outdir}/PRS_CS_results/
    """  
}

///////// SBayesR /////////
// NEED UPDATE
process SBayesR_step {

    script:
    """
    mkdir -p ${params.outdir}/SBayesR_results

    gctb --sbayes R --ldm ${params.sbayesr_ldm} --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ${params.sbayesr_gwas_summary} --chain-length 10000 --burn-in 2000 --out-freq 10 --out ${params.outdir}/SBayesR_results/sim_1
    """  
}