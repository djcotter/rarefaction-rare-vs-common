# NYGC 1000 Genomes High Coverage Samples
## 11 May 2021

# Analysis
This `snakemake` pipeline will download and filter each vcf file from the new New York Genome Center for biallelic SNPS that pass the GATK filter. These filtered files will be stored in the `data/` directory.

## Steps
### Step 1: Create Conda Environment
Download conda](https://conda.io/docs/user-guide/install/index.html) (miniconda3) is usually best and create an environment using the provided `environment.yml` file.
```shell
conda env create -f environment.yml
```

### Step 2: Run the Snakemake workflow
Note: this should probably be done on the cluster. It can be parallelized using cluster submission commands.
```shell
conda activate snakemake_default
snakemake -j 1
```
