#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH -t 24:00:00
#SBATCH -p owners
#usage: the directory that you want to upload to google drive
module spider gdrive
ml load system
ml load gdrive
gdrive upload --recursive $1
