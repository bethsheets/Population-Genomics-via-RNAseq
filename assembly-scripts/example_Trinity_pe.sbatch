#!/bin/bash 
#SBATCH -p spalumbi,hns
#################
#set a job name  
#SBATCH --job-name=Astroides_GM
#################  
#a file for job output, you can check job progress
#SBATCH --output=Astroides_GM
#################
# a file for errors from the job
#SBATCH --error=Astroides_GM.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=01:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#SBATCH --ntasks=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-tasks=16" instead
#################
#SBATCH -p spalumbi,hns
#SBATCH --dependency=singleton
#SBATCH -c 8

module load biology
module load trinity/2.8.4/

Trinity --seqType fq --max_memory 250G --CPU 8 --min_contig_length 300 --trimmomatic --left <PATH TO FILES_1.fq>, --right <PATH TO FILES_2.fq>,
