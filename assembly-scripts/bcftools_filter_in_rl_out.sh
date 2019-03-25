#!/bin/bash
bcftools view -g ^miss -q 0.05 -Q 0.95 -m2 -M2 $1 | bcftools filter -g $2 -i 'TYPE="snp"' > $3
