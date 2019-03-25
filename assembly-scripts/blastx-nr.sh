#!/bin/bash
#SBATCH -p owners 
#SBATCH --time=48:00:00
#SBATCH --mem=32000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
########################
# -outfmt 5 XML
# to call: sbatch blastx-nr.sh input.fa databasepath

echo $1
#blast against blastnr
/home/groups/spalumbi/programs/ncbi-blast-2.8.1+/bin/blastx -db $2 -query $1 -out $1.blast.out -evalue 0.001 -max_hsps 1 -num_threads 12 -outfmt 5

