#!/bin/bash
#SBATCH --partition spalumbi,owners
#SBATCH --mem 16000

bwa index $1 -p $2
