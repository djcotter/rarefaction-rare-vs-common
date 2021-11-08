# High-coverage 1000 Genomes Data

This directory contains data from the higher-coverage sequencing of the
individuals in the 1000 Genomes Project. The sequencing was conducted by the
New York Genome Center (NYGC) and variants have been genotyped using GATK v4.0.

NOTE : this dataset consists of data on genome build hg38 (not hg19). If you
are interested in repositioning variants to hg19 coordinates use LiftOver.

Original data can be found at this URL:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/

In this directory we have filtered to biallelic single-nucleotide variants
(SNVs) on the autosomes using the `bcftools` command:
```
bcftools view -v snps -m2 -M2 -i '%FILTER == "PASS"'
```
The files either have:

* .filt.vcf.gz : a version of the VCF file that requires the `FILTER==PASS`,
limiting the inclusion of poorly called genotypes

For any questions regarding data generation or error modes, contact Arjun
Biddanda via Slack
or at aabiddanda@gmail.com
