#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=fastqc
#################  
#a file for job output, you can check job progress
#SBATCH --output=barnacle_fastqc.out
#################
# a file for errors from the job
#SBATCH --error=barnacle_fastqc.err
#################
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#SBATCH --ntasks=1
#################
#memory per node; default is 4000 MB per CPU

/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X9_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X9_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X9_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X9_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X5_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X5_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X5_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X5_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X8_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X8_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X8_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X8_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X16_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X16_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X16_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X16_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X12_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X12_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X12_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X12_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X4_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X4_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X4_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X4_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X7_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X7_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X7_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X7_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X15_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X15_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X15_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X15_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X1_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X1_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X1_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X1_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X3_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X3_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X3_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X3_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X14_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X14_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X14_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X14_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X6_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X6_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X6_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X6_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X11_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X11_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X11_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X11_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X2_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X2_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X2_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X2_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X13_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X13_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X13_160204_D00550_0321_BC8RYTANXX_7_2_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X13_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X10_160204_D00550_0321_BC8RYTANXX_7_1_paired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X10_160204_D00550_0321_BC8RYTANXX_7_1_unpaired.fastq
/share/PI/spalumbi/programs/FastQC/fastqc -o qc/ --extract -f fastq 12132X10_160204_D00550_0321_BC8RYTANXX_7_2_unpaired.fastq
