# Rarefaction Project for Common and Rare Alleles
## 11 May 2021

# Analysis
This `snakemake` pipeline will filter each vcf file from the new New York Genome Center for biallelic SNPS that pass the GATK filter. These filtered files will be stored in the `data/tmp` directory.

The pipeline will then convert the filtered VCF files into raw allele counts and organize these allele counts by site and by population. These summary tables will be stored in `data/allele_counts`

## Steps
### Step 1: Create Conda Environment
Download [conda](https://conda.io/docs/user-guide/install/index.html) (miniconda3 is usually best) and create an environment using the provided `environment.yml` file.
```shell
conda env create -f environment.yml
```

### Step 2: Point the pipeline at the data directory
Specify the `data_path` flag in the `config.yml` file so the pipeline knows where to look for the data. The data can be downloaded in bulk from the 1kg project via a tool such as Globus. It is available at [the 1kg FTP site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/).

### Step 2: Run the Snakemake workflow
Note: this should probably be done in parallel to speed up the generation of all allele counts tables.
```shell
conda activate snakemake_default
snakemake -j 1
```
It can be parallelized by changing `j>1` and/or via using cluster submission commands.
