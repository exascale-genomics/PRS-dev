
// parameters
params.base_data = 'Height.gwas.txt.gz'
params.target_data = '/home/ubuntu/inputs_target'
params.outdir = '/home/ubuntu/outputs'
params.workingDir = System.getProperty("user.dir")

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
                     .map { a -> [checker(a[1]), checker(a[2]), checker(a[3]), checker(a[4]), checker(a[5])] }\
                     .multiMap { it -> copy1: copy2: copy3: copy4: copy5: copy6: copy7: copy8: copy9: copy10: copy11: it }
                     .set { target_data }

///////// base data QC /////////
// base data QC: File transfer
process base_qc_file_transfer {

    input:
    file base_data

    output:
    stdout into base_qc_file_transfer_output

    """
    md5sum $base_data
    """
}

base_qc_file_transfer_output.view { "base md5sum:\n $it" }

// base data QC: Standard GWAS QC
process base_qc_standard_gwas_qc {

    publishDir "${params.outdir}", mode: 'copy'

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

    publishDir "${params.outdir}", mode: 'copy'

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

    publishDir "${params.outdir}", mode: 'copy'

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

    publishDir "${params.outdir}", mode: 'copy'

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

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy1

    output:
    stdout into target_qc_file_transfer_output

    """
    md5sum $bed
    md5sum $bim
    md5sum $cov
    md5sum $fam
    md5sum $height
    """
}

target_qc_file_transfer_output.view { "target md5sum:\n $it" }

// target data QC: Standard GWAS QC
//// removing SNPs with low genotyping rate, low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing individuals with low genotyping rate
//// perform pruning to remove highly correlated SNPs
//// compute heterozygosity rates
process target_qc_standard_gwas_qc_1 {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy2

    output:
    tuple file("${output_base}.irem"), file("${output_base}.hh"), file("${output_base}.snplist"), \
    file("${output_base}.fam"), file("${output_base}.prune.in"), file("${output_base}.prune.out"), \
    file("${output_base}.het"), file("${output_base}.log") into target_qc_standard_gwas_qc_1_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    plink --bfile ${base} --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --make-just-fam --out ${output_base}
    plink --bfile ${base} --keep ${output_base}.fam --extract ${output_base}.snplist --indep-pairwise 200 50 0.25 --out ${output_base}
    plink --bfile ${base} --extract ${output_base}.prune.in --keep ${output_base}.fam --het --out ${output_base}
    """
}

target_qc_standard_gwas_qc_1_output.into { target_qc_standard_gwas_qc_1_output_copy1; target_qc_standard_gwas_qc_1_output_copy2; target_qc_standard_gwas_qc_1_output_copy3; target_qc_standard_gwas_qc_1_output_copy4; target_qc_standard_gwas_qc_1_output_copy5 }

//// remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean
process target_qc_standard_gwas_qc_2 {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy3
    tuple file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het), file(log1) from target_qc_standard_gwas_qc_1_output_copy1

    output:
    file("${base}.valid.sample") into target_qc_standard_gwas_qc_2_output

    script:
    base = bed.baseName
    """
    #!/usr/bin/env Rscript

    dat <- read.table("${het}", header=T) # Read in the EUR.het file, specify it has header
    m <- mean(dat\$F) # Calculate the mean  
    s <- sd(dat\$F) # Calculate the SD
    valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
    write.table(valid[,c(1,2)], "${base}.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
    """
}

// target data QC: Mismatching SNPs
process target_qc_mismatching_snps {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy4
    tuple file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het), file(log1) from target_qc_standard_gwas_qc_1_output_copy2
    file base_qc_final from base_qc_final_output

    output:
    tuple file("${base}.a1"), file("${base}.mismatch") into target_qc_mismatching_snps_output

    script:
    base = bed.baseName
    """
    R --file=${bin}/target_qc_mismatching_snps.R --args ${bim} ${snplist} ${base_qc_final} ${base}.a1 ${base}.mismatch
    """  
}

// target data QC: Sex chromosomes
process target_qc_sex_check {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy5
    tuple file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het), file(log1) from target_qc_standard_gwas_qc_1_output_copy3
    file valid_sample from target_qc_standard_gwas_qc_2_output

    output:
    tuple file("${output_base}.sexcheck"), file("${output_base}.valid"), file("${output_base}.hh") into target_qc_sex_check_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    plink --bfile ${base} --extract ${prune_in} --keep ${valid_sample} --check-sex --out ${output_base}
    R --file=${bin}/target_qc_sex_check.R --args ${valid_sample} ${output_base}.sexcheck ${output_base}.valid
    """  
}

// target data QC: Relatedness
process target_qc_relatedness {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy6
    tuple file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het), file(log1) from target_qc_standard_gwas_qc_1_output_copy4
    tuple file(sexcheck), file(valid), file("tmp_hh") from target_qc_sex_check_output

    output:
    file("${output_base}.rel.id") into target_qc_relatedness_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    plink --bfile ${base} --extract ${prune_in} --keep ${valid} --rel-cutoff 0.125 --out ${output_base}
    """  
}

// target data QC: Generate final QC'ed target data file
process target_qc_final {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy7
    tuple file(irem), file(hh), file(snplist), file(fam1), file(prune_in), file(prune_out), file(het), file(log1) from target_qc_standard_gwas_qc_1_output_copy5
    tuple file(a1), file(mismatch) from target_qc_mismatching_snps_output
    file rel_id from target_qc_relatedness_output

    output:
    tuple file("${output_base}.fam"), file("${output_base}.bed"), file("${output_base}.bim") into target_qc_final_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    plink --bfile ${base} --make-bed --keep ${rel_id} --out ${output_base} --extract ${snplist} --exclude ${mismatch} --a1-allele ${a1}
    """  
}

target_qc_final_output.into { target_qc_final_output_copy1; target_qc_final_output_copy2; target_qc_final_output_copy3 }

///////// PRS using plink /////////
// Update Effect Size
process update_effect_size {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file base_qc_final from base_qc_final_output

    output:
    file("${output_file}") into update_effect_size_output

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
process clumping {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy8
    tuple file(fam1), file(bed1), file(bim1) from target_qc_final_output_copy1
    file transformed from update_effect_size_output


    output:
    tuple file("${base}.clumped"), file("${base}.valid.snp") into clumping_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    plink --bfile ${output_base} --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ${transformed} --clump-snp-field SNP --clump-field P --out ${base}
    awk 'NR!=1{print \$3}' ${base}.clumped >  ${base}.valid.snp
    """
}

// Generate PRS
process generate_prs {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy9
    tuple file(fam1), file(bed1), file(bim1) from target_qc_final_output_copy2
    tuple file(clumped), file(valid_snp) from clumping_output
    file transformed from update_effect_size_output

    output:
    tuple file("${base}.0.5.profile"), file("${base}.0.4.profile"), file("${base}.0.3.profile"), file("${base}.0.2.profile"), file("${base}.0.1.profile"), file("${base}.0.05.profile"), file("${base}.0.001.profile") into generate_prs_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    awk '{print \$3,\$8}' ${transformed} > SNP.pvalue
    echo "0.001 0 0.001" > range_list 
    echo "0.05 0 0.05" >> range_list
    echo "0.1 0 0.1" >> range_list
    echo "0.2 0 0.2" >> range_list
    echo "0.3 0 0.3" >> range_list
    echo "0.4 0 0.4" >> range_list
    echo "0.5 0 0.5" >> range_list
    plink --bfile ${output_base} --score ${transformed} 3 4 12 header --q-score-range range_list SNP.pvalue --extract ${valid_snp} --out ${base}
    """
}

// Accounting for Population Stratification
process accounting_for_population_stratification {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy10
    tuple file(fam1), file(bed1), file(bim1) from target_qc_final_output_copy3

    output:
    tuple file("${base}.prune.in"), file("${base}.prune.out"), file("${base}.eigenvec"), file("${base}.eigenval") into accounting_for_population_stratification_output

    script:
    base = bed.baseName
    output_base = "${base}.QC"
    """
    # First, we need to perform prunning
    plink --bfile ${output_base} --indep-pairwise 200 50 0.25 --out ${base}
    # Then we calculate the first 6 PCs
    plink --bfile ${output_base} --extract ${base}.prune.in --pca 6 --out ${base}
    """
}

// Finding the "best-fit" PRS
process finding_the_best_fit_prs {

    input:
    tuple file(bed), file(bim), file(cov), file(fam), file(height) from target_data.copy11
    tuple file(prune_in), file(prune_out), file(eigenvec), file(eigenval) from accounting_for_population_stratification_output
    tuple file(profile1), file(profile2), file(profile3), file(profile4), file(profile5), file(profile6), file(profile7) from generate_prs_output

    output:
    stdout finding_the_best_fit_prs_output

    script:
    base = bed.baseName
    """
    R --file=${bin}/regression_between_prs.R --args ${height} ${eigenvec} ${cov} ${base}
    """  
}

finding_the_best_fit_prs_output.view()
