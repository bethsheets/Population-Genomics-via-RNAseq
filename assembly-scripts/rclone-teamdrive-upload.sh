#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH -t 24:00:00
#SBATCH -p owners

ml system rclone
rclone

