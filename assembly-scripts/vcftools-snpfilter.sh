#!/bin/bash
#SBATCH -p owners,spalumbi
#SBATCH -c 4

#usage: sbatch vcftools-snpfilter.sh combined.vcf
# max-missing 1 = no missing

vcftools --vcf $1 --remove-indels --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --minGQ 30 --minDP 10 --max-missing 1 --max-maf 0.99 --maf 0.01 --out $(basename $1 .vcf)_minMAF01_GQ30_DP10_noNA
