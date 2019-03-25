#!/bin/bash
#SBATCH -p owners,spalumbi
#SBATCH -c 4

#created by Beth March 2019
#before using this script, you need add this path to the configfile in your home .bashrc file: 
# export BUSCO_CONFIG_FILE="/home/groups/spalumbi/programs/busco-master/config/config.ini"

#usage
#busco.sh inputfile outputfile genome

#genome= genome, transcriptome, proteins
# I didn't install the Augustus package needed for genome assessment yet. You will need to do that to use this script on a genome.


ml biology py-busco
run_BUSCO.py -i $1 -o $2 -m $3 -l /home/groups/spalumbi/programs/busco-master/metazoa_odb9/

