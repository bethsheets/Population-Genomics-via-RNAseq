#!/bin/bash
#SBATCH -p owners
#SBATCH --mem 8000

#usage: sbatch vcf-filter-nomissing-maf05-allgq30.sh #samples outfilename input.vcf

#made by Noah & Beth April 2016

#we filter for snps with no NAs, a minimum MAF=0.05, and across all individuals SNP probability <0.001
#had to add an extra filter to take out sites with no genotype calls after genotype filtering, which led to nans in the AF field
vcffilter -g "GQ > 30" $3 | vcffixup - | vcffilter -f "AC > 0" | vcffilter -f "TYPE = snp & NS = $1 & AF > 0.05 & AF < 0.95" >  ${2}_GQ30_MAF05_NONA.vcf
