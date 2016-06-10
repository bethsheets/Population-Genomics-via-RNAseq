#!/bin/bash
#SBATCH -p owners
#SBATCH --mem 8000

#usage: sbatch vcf-filter-nomissing-maf05-qual30.sh #samples outfilename *.vcf

#made by Noah & Beth April 2016

#we filter for snps with no NAs, a minimum MAF=0.05, and a locus-wide SNP probability <0.001
vcffilter -f "TYPE = snp & QUAL > 30 & NS = $1 & AF > 0.05 & AF < 0.95" $3 > $2_QUAL30_MAF05_NONA.vcf

#can run this in parallel with genotype-based quality cutoff to see results of different filtering techniques
