"""
rearefaction-project
Snakefile
Daniel Cotter

analyze geographic distribution of alleles, correcting for sample size
---------------------------------------------------------------------

Requires:
    conda
        to activate the included environment.yml environment
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
G_LIST = list(range(10,501,10))

# global variables
CHROMS = [x for x in range(1, 23)]  # list from 1 to 22

# helper functions -----------------
def plot_range(wildcards):
    """
    Function to analyze the range wildcard and define a flag
    for the plotting scripts.
    """
    if wildcards.range == "wholeChr":
        return("") 
    else:
        coords = wildcards.range.split("-")
        return(f"--range {coords[0]}:{coords[1]}")

# declare rules -------------------
rule all:
    input:
        fig1 = expand(path.join('figures', '{chr}_patterns_{sample_size}-snps_by_g.pdf'),
                      chr=22,
                      sample_size="all"),
        fig2 = expand(path.join('figures', '{chr}_g-{g}_{sample_size}-snps_{singletons}Singletons_network.{ext}'),
                      chr=22,
                      g=500,
                      singletons="no",
                      sample_size="all",
                      ext="pdf"),
        fig4 = expand(path.join('figures', '{chr}_pattern-match-probs_{sample_size}-snps.pdf'),
                      chr=22,
                      sample_size="all"),
        fig5 = expand(path.join('figures', 
                                '{chr}_g-{g}_pattern_byPosition_100kb-windows_{sample_size}-snps'
                                '_{singletons}Singletons_{range}.{ext}'),
                      chr=22,
                      g=500,
                      sample_size="all",
                      singletons="no",
                      range="wholeChr",
                      ext="pdf"),
        fig6a = expand(path.join('figures', 
                                '{chr}_g-{g}_pattern_byPosition_100kb-windows_{sample_size}-snps'
                                '_{singletons}Singletons_{range}.{ext}'),
                       chr=6,
                       g=500,
                       sample_size="all",
                       singletons="no",
                       range="20-40",
                       ext="pdf"),
        fig6b = expand(path.join('figures', 
                                '{chr}_g-{g}_pattern_byPosition_byRank_100kb-windows_{sample_size}-snps'
                                '_{singletons}Singletons_{range}.{ext}'),
                       chr=6,
                       g=500,
                       sample_size="all",
                       singletons="no",
                       range="20-40",
                       ext="pdf")


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
    a subsample of the chromosome size.
    """
    input:
        path.join('data', 'allele_counts', 'chr{chr}_counts_superpops.txt')
    params:
        script = path.join('src', 'calculate_superpop_allele_patterns.R'),
        chr = lambda wildcards: wildcards.chr,
        threshold = 0.05,
        sample = lambda wildcards: 0 if wildcards.sample_size == 'all' else wildcards.sample_size,
        singletons = lambda wildcards: '--drop-singletons' if wildcards.singletons == "no" else '',
        cores = 5
    output:
        path.join('data', 'patterns', '{chr}_patterns_{sample_size}-snps_{singletons}Singletons.txt'),
        path.join('data', 'patterns', 
                  '{chr}_pattern-match-proportions_{sample_size}-snps_{singletons}Singletons.txt'),
        path.join('data', 'patterns',
                  '{chr}_actualPattern_{sample_size}-snps_{singletons}Singletons.txt')

    shell:
        "Rscript --vanilla {params.script} --chr {params.chr} --threshold {params.threshold} "
        "--sample {params.sample} --ncores {params.cores} {params.singletons}"

rule calculate_allele_patterns_byPosition:
    """
    This rule takes the allele counts table for a given chromosome and calculates the 
    five-letter super population patterns for each locus. It outputs the probability
    of each pattern for all snps analyzed. It can operate on the whole chromosome or
    a subsample of the chromosome size.
    """
    input:
        path.join('data', 'allele_counts', 'chr{chr}_counts_superpops.txt')
    params:
        script = path.join('src', 'calculate_superpop_allele_patterns_byPosition.R'),
        chr = lambda wildcards: wildcards.chr,
        threshold = 0.05,
        sample = lambda wildcards: 0 if wildcards.sample_size == 'all' else wildcards.sample_size,
        singletons = lambda wildcards: '--drop-singletons' if wildcards.singletons == "no" else '',
        g = lambda wildcards: wildcards.g
    output:
        path.join('data', 'patterns', '{chr}_g-{g}_pattern_byPosition_{sample_size}-snps_{singletons}Singletons.txt')
    shell:
        "Rscript --vanilla {params.script} --chr {params.chr} --threshold {params.threshold} "
        "--sample {params.sample} --g_size {paramgs.g} {params.singletons}"

rule prepare_network_data:
    """
    Take the output of the pattern calculations and process it for plotting the network figure
    """
    input:
        path.join('data', 'patterns', '{chr}_patterns_{sample_size}-snps_{singletons}Singletons.txt')
    params:
        script = path.join('src', 'prepare_network_pattern_data.R'),
        g = lambda wildcards: wildcards.g
    output:
        path.join('data', 'tmp', '{chr}_g-{g}_{sample_size}-snps_{singletons}Singletons.txt')
    shell: 
        "Rscript --vanilla {params.script} --input {input} --g_filter {params.g} --output {output}"

rule plot_pattern_network:
    """
    Plot a 21 node network displaying the proportions of pattern probabilities across a chromosome
    """
    input:
        path.join('data', 'tmp', '{chr}_g-{g}_{sample_size}-snps_{singletons}Singletons.txt')
    params:
        script = path.join('src', 'network_drawing.py')
    output:
        path.join('figures', '{chr}_g-{g}_{sample_size}-snps_{singletons}Singletons_network.{ext}')
    shell:
        "python {params.script} {input} {output}"

rule plot_patterns_vs_g:
    """
    Take the output of the rule calculate_various_g_allele_patterns and plot a four panel
    figure displaying how geographic patterns change as a function of the pop sample size g.
    R script.
    """
    input:
        wSingletons = path.join('data', 'patterns',
                                '{chr}_patterns_{sample_size}-snps_wSingletons.txt'),
        noSingletons = path.join('data', 'patterns',
                                 '{chr}_patterns_{sample_size}-snps_noSingletons.txt'),
        actual_wSingletons = path.join('data', 'patterns',
                                       '{chr}_actualPattern_{sample_size}-snps_wSingletons.txt'),
        actual_noSingletons = path.join('data', 'patterns',
                                       '{chr}_actualPattern_{sample_size}-snps_noSingletons.txt')
    params:
        script = path.join('src', 'plot_patterns_vs_sample_size.R')
    output:
        plot = path.join('figures', '{chr}_patterns_{sample_size}-snps_by_g.pdf'),
        legend = path.join('figures', '{chr}_patterns_{sample_size}-snps_by_g_legend.pdf')
    shell:
        "Rscript --vanilla {params.script} --wSingletons {input.wSingletons} "
        "--noSingletons {input.noSingletons} --actual-wSingletons {input.actual_wSingletons} "
        "--actual-noSingletons {input.actual_noSingletons} --output {output.plot}"

rule plot_match_rate:
    """
    Take the match probabilities from the rule calculate_various_g_allele_patterns to plot
    how the prbobability of matching the empirical data changes as a function of the population
    sample size g
    """
    input:
        wSingletons = path.join('data', 'patterns', 
                  '{chr}_pattern-match-proportions_{sample_size}-snps_wSingletons.txt'),
        noSingletons = path.join('data', 'patterns', 
                  '{chr}_pattern-match-proportions_{sample_size}-snps_noSingletons.txt')
    params:
        script = path.join('src', 'plot_match_probs.R')
    output:
        path.join('figures', '{chr}_pattern-match-probs_{sample_size}-snps.pdf')
    shell:
        "Rscript --vanilla {params.script} --wSingletons {input.wSingletons} "
        "--noSingletons {input.noSingletons} --output {output}"

rule plot_probs_byPosition:
    """
    Takes a file with probabilities across a chromosome and recodes the patterns as summaries.
    Then plots the mean probabilities of these summaries in 100kb windows across the chromosome
    or across a designated range.
    """
    input:
        path.join('data', 'patterns', '{chr}_g-{g}_pattern_byPosition_{sample_size}-snps_{singletons}Singletons.txt')
    params:
        script = path.join('src', 'plot_genomic_positions.R'),
        plot_range = plot_range
    output:
        path.join('figures', '{chr}_g-{g}_pattern_byPosition_100kb-windows_{sample_size}-snps_{singletons}Singletons_{range}.{ext}')
    shell:
        "Rscript --vanilla {params.script} --input {input} "
        "--output {output} {params.plot_range}"

rule plot_ranks_byPosition:
    """
    Takes a file with probabilities across a chromosome and recodes the patterns as summaries.
    Then plots the ranks of these summaries in each 100kb window across the chromosome
    or across a designated range.
    """
    input:
        path.join('data', 'patterns', '{chr}_g-{g}_pattern_byPosition_{sample_size}-snps_{singletons}Singletons.txt')
    params:
        script = path.join('src', 'plot_genomic_positions_ranks.R'),
        plot_range = plot_range,
        rank_cutoff = 2
    output:
        path.join('figures',
                  '{chr}_g-{g}_pattern_byPosition_byRank_100kb-windows_{sample_size}-snps_{singletons}Singletons_{range}.{ext}')
    shell:
        "Rscript --vanilla {params.script} --input {input} "
        "--output {output} --rank_cutoff {params.rank_cutoff} "
        "{params.plot_range}"