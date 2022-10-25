"""
chrX_regional_variation
Snakefile
Daniel Cotter

analyze diversity and LD across chrX and chrY from 1000 genomes data
---------------------------------------------------------------------

Requires:
    conda
        to activate the included PAB_variation.yml environment
    snakemake
        included in the Conda environment
"""

# Import statements -----------------------------------------------------------
import json
from os import path

# Global configurations -------------------------------------------------------
# declare a path to the configuration file
configfile: 'config.yml'

# parse all population codes
SUPERPOPULATIONS = sorted(json.load(open(config['POP_CODES']))['SUPERPOP'])
POPULATIONS = sorted(json.load(open(config['POP_CODES']))['POP'])

POP_LIST = ['ALL'] + SUPERPOPULATIONS + POPULATIONS
G_LIST = list(range(10,301,10))

# global variables
CHROMS = [x for x in range(1, 23)]  # list from 1 to 22

# declare rules -------------------
rule all:
    input:
        expand(path.join('data',
                         'allele_counts',
                         'chr{CHR}_counts_{CAT}.txt'),
               CHR=CHROMS, CAT=['superpops'])

rule filter_raw_data:
    """
    This file runs bcftools and points at the FTP server for 1KGP data.
    bcftools view filters for SNPS (-v snps) that are biallelic (-m2 -M2) and that
    pass the GATK quality threshold (-i '%FILTER=="PASS"'). -Ou outputs an uncompressed
    bcf file to the next command. bcftools annotate strips off extra information
    associated with each snp and each genotype (-x INFO,^FORMAT/GT), renames the chromosomes
    to drop the 'chr' prefix (--rename-chrs FILE), and sets the identifier column to the snp
    position (--set-id %POS).
    """
    input:
        vcf = path.join(config['data_path'],
                        'CCDG_13607_B01_GRM_WGS_2019-02-19_chr{chr}.recalibrated_variants.vcf.gz'),
        rename = path.join('data', 'ref_files', 'rename_chrs.txt')
    output:
        out = path.join('data', 'tmp', '1kg_nygc_chr{chr}_biallelic.snps_filt.vcf.gz')
    params:
        filter = "\'%FILTER==\"PASS\"\'"
    shell:
        "bcftools view -Ou -v snps -m2 -M2 -i {params.filter} {input.vcf} | "
        "bcftools annotate -Oz -x INFO,^FORMAT/GT --rename-chrs {input.rename} "
        "--set-id %POS > {output.out}"

rule get_pop_panel:
    """
    This rule downloads the list of 2504 individuals from phase3 of the 1000 genomes project.
    The fields in this file are sample, pop, super_pop, and sex. The file HAS a header.
    """
    params:
        URL = config['pop_panel']
    output:
        path.join('data', 'tmp', 'individual_population_codes.txt')
    shell:
        "wget -O {output} {params.URL}"

rule parse_populations:
    """
    This rule parses the list of 2504 individuals from phase3 of the 1000 genomes project.
    It separates the file into 1 (ALL) + 5 (each superpopulation) + 26 (each population) lists
    containing a 1 column file of the samples in each population or superpopulation
    """
    input:
        path.join('data', 'tmp', 'individual_population_codes.txt')
    params:
        script = path.join('src', 'parse_populations.R')
    output:
        expand(path.join('data', 'pops', '{POP}_samples.txt'),
               POP=POP_LIST)
    shell:
        "Rscript {params.script} --bid"

rule count_alleles:
    """
    This rule uses vcftools to count the number of alleles for ALL individuals and then
    for each of the 5 superpopulations and 26 populations.
    """
    input:
        vcf = path.join('data', 'tmp', '1kg_nygc_chr{chr}_biallelic.snps_filt.vcf.gz'),
        samples = path.join('data', 'pops', '{population}_samples.txt')
    params:
        outpath = lambda wildcards: expand(path.join('data', 'tmp',
                                                     '{POP}_chr{CHR}'),
                                           POP=wildcards.population,
                                           CHR=wildcards.chr)
    output:
        path.join('data', 'tmp', '{population}_chr{chr}.frq.count')
    log:
        temp(path.join('data', 'tmp', '{population}_chr{chr}.log'))
    shell:
        "vcftools --gzvcf {input.vcf} --keep {input.samples} "
        "--counts --out {params.outpath}"

rule merge_POP_allele_counts:
    """
    This rule takes the vcftools output for each POP and
    assigns the globally minor allele then organizes the allele counts
    into a single column in the format MINOR/MAJOR where minor and major are given in
    columns 4 and 5 respectively. Columns 1-3 are chromosome, position, and total alleles.
    Columns 6+ correspond to each POP's allele counts
    """
    input:
        lambda wildcards: expand(
            path.join('data', 'tmp', '{POP}_chr{CHR}.frq.count'),
            POP=POPULATIONS,
            CHR=wildcards.chr
        )
    params:
        script = path.join('src', 'merge_allele_counts.R'),
        chr = lambda wildcards: wildcards.chr
    output:
        path.join('data', 'allele_counts', 'chr{chr}_counts_pops.txt'),
        path.join('data', 'allele_counts', 'chr{chr}_equal-frequency-alleles_pops.txt')
    shell:
        "Rscript --vanilla {params.script} --pops --chr {params.chr}"

rule merge_SUPERPOP_allele_counts:
    """
    This rule takes the vcftools output for each SUPERPOP and
    assigns the globally minor allele then organizes the allele counts
    into a single column in the format MINOR/MAJOR where minor and major are given in
    columns 4 and 5 respectively. Columns 1-3 are chromosome, position, and total alleles.
    Columns 6+ correspond to each SUPERPOP's allele counts
    """
    input:
        lambda wildcards: expand(
            path.join('data', 'tmp', '{POP}_chr{CHR}.frq.count'),
            POP=SUPERPOPULATIONS,
            CHR=wildcards.chr
        )
    params:
        script = path.join('src', 'merge_allele_counts.R'),
        chr = lambda wildcards: wildcards.chr
    output:
        path.join('data', 'allele_counts', 'chr{chr}_counts_superpops.txt'),
        path.join('data', 'allele_counts', 'chr{chr}_equal-frequency-alleles_superpops.txt'),
        path.join('data', 'allele_counts', 'chr{chr}_missing-data_superpops.txt')
    shell:
        "Rscript --vanilla {params.script} --superpops --chr {params.chr}"

rule calculate_various_g_allele_patterns:
    """
    This rule takes the allele counts table for a given chromosome and calculates the 
    five-letter super population patterns for each locus. It outputs the mean probability
    of each pattern across all snps analyzed. It can operate on the whole chromosome or
    a subset of the chromosome size.
    """
    input:
        path.join('data', 'allele_counts', 'chr{chr}_counts_superpops.txt')
    params:
        script = path.join('src', 'calculate_superpop_allele_patterns.R'),
        chr = lambda wildcards: wildcards.chr,
        threshold = 0.05,
        sample = lambda wildcards: 0 if wildcards.sample_size == 'all' else wildcards.sample_size,
        singletons = lambda wildcards: '--drop-singletons' if wildcards.singletons == "no" else ''
    output:
        path.join('data', 'patterns', '{chr}_patterns_{sample_size}-snps_{singletons}Singletons.txt'),
        path.join('data', 'patterns', 
                  '{chr}_pattern-match-proportions_{sample_size}-snps_{singletons}Singletons.txt'),
        path.join('data', 'patterns',
                  '{chr}_actualPattern_{sample_size}-snps_{singletons}Singletons.txt')

    shell:
        "Rscript --vanilla {params.script} --chr {params.chr} --threshold {params.threshold} "
        "--sample {params.sample} {params.singletons}"