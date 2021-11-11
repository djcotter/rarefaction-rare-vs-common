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
from os import path

# Global configurations -------------------------------------------------------
# declare a path to the configuration file
configfile: 'config.yml'

# global variables
CHROMS = [x for x in range(1, 23)]  # list from 1 to 22

# declare rules -------------------
rule all:
    input:
        expand(path.join('data',
                         '1kg_nygc_chr{chr}_biallelic.snps_filt.vcf.gz'),
               chr=CHROMS),
        path.join('data', 'individual_population_codes.txt')

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
        rename = path.join('data', 'ref_files', 'rename_chrs.txt')
    output:
        out = path.join('data', '1kg_nygc_chr{chr}_biallelic.snps_filt.vcf.gz')
    params:
        data_path = lambda wildcards: path.join(config['data_path'],
                                                'CCDG_13607_B01_GRM_WGS_2019-02-19_chr'
                                                '{}.recalibrated_variants.vcf.gz'.format(wildcards.chr)),
        filter = "\'%FILTER==\"PASS\"\'"
    shell:
        "bcftools view -Ou -v snps -m2 -M2 -i {params.filter} {params.data_path} | "
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
        path.join('data', 'individual_population_codes.txt')
    shell:
        "wget -O {output} {params.URL}"
