#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH -t 24:00:00
#SBATCH -p spalumbi,owners
#usage: sbatch rclone.sh sherlockpath palumbiteam:directory/
ml system rclone
rclone --checkers=16 --drive-chunk-size=16384k copy $1 $2
