# NYGC 1000 Genomes High Coverage Samples
## 11 May 2021

# Analysis
This `snakemake` pipeline will filter each vcf file from the new New York Genome Center for biallelic SNPS that pass the GATK filter. These filtered files will be stored in the `data/tmp` directory.

The pipeline will then convert the filtered VCF files into raw allele counts and organize these allele counts by site and by population. These summary tables will be stored in `data/allele_counts`

## Steps
### Step 1: Create Conda Environment
Download conda](https://conda.io/docs/user-guide/install/index.html) (miniconda3) is usually best and create an environment using the provided `environment.yml` file.
```shell
conda env create -f environment.yml
```

### Step 2: Point the pipeline at the data directory
Specify the `data_path` flag in the `config.yml` file so the pipeline knows where to look for the data. The data can be downloaded in bulk from 1kg project via a tool such as Globus.

### Step 2: Run the Snakemake workflow
Note: this should probably be done on the cluster.
```shell
conda activate snakemake_default
snakemake -j 1
```
It can be parallelized by changing `j>1` and/or via using cluster submission commands.
