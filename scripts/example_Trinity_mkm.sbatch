#!/bin/bash 
#################
#set a job name  
#SBATCH --job-name=LOPaired
#################  
#a file for job output, you can check job progress
#SBATCH --output=LOPaired.out
#################
# a file for errors from the job
#SBATCH --error=LOPaired.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=48:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#SBATCH --ntasks=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=160000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-tasks=16" instead
#################
#SBATCH -p hns
#SBATCH --dependency=singleton
#SBATCH -c 8

/home/mkm44/Programs/trinityrnaseq-2.1.1/Trinity --seqType fq --max_memory 160G --CPU 8 --left /scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO01_03_H_69_P1_TrimClip_FLASH.notCombined_1.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO02_03_H_70_P1_TrimClip_FLASH.notCombined_1.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO03_03_R_75_P1_TrimClip_FLASH.notCombined_1.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO01_03_H_69_P1_TrimClip_FLASH.extendedFrags.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO02_03_H_70_P1_TrimClip_FLASH.extendedFrags.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO03_03_R_75_P1_TrimClip_FLASH.extendedFrags.fastq,/scratch/users/mkm44/Multispecies/30bp/30bp_LO01_03_H_69_P1_TrimClip_forward_unpaired.fq,/scratch/users/mkm44/Multispecies/30bp/30bp_LO01_03_H_69_P1_TrimClip_reverse_unpaired.fq,/scratch/users/mkm44/Multispecies/30bp/30bp_LO02_03_H_70_P1_TrimClip_forward_unpaired.fq,/scratch/users/mkm44/Multispecies/30bp/30bp_LO02_03_H_70_P1_TrimClip_reverse_unpaired.fq,/scratch/users/mkm44/Multispecies/30bp/30bp_LO03_03_R_75_P1_TrimClip_forward_unpaired.fq,/scratch/users/mkm44/Multispecies/30bp/30bp_LO03_03_R_75_P1_TrimClip_reverse_unpaired.fq --right /scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO01_03_H_69_P1_TrimClip_FLASH.notCombined_2.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO02_03_H_70_P1_TrimClip_FLASH.notCombined_2.fastq,/scratch/users/mkm44/Multispecies/30bp/merged/30bp_LO03_03_R_75_P1_TrimClip_FLASH.notCombined_2.fastq
