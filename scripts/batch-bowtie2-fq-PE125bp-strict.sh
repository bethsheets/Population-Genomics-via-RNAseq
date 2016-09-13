#!/bin/bash
#USAGE: bash batch-bowtie2-fq-PE125bp-strict.sh bt2index chunksize *.fastq
#if you don't have a bowtie2 index, build it with "bowtie2-build <reference>.fa basename"
#for 125bp PE reads, allows 4 mismatches per 125 bp. The x value in the --score-min is based on 125 bp length, so the current x value would not be optimum with FLASH merged PE reads. 
CHUNK=$2
COUNTER=0
FQ="${@:3}"
for i in $FQ; do
    if [ $COUNTER -eq 0 ]; then
    echo -e "#!/bin/bash\n#SBATCH -p owners\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=3\n#SBATCH -t 12:00:00\n#SBATCH --mem 24000" > TEMPBATCH.sbatch; fi
    BASE=$( basename $i _1_paired.fastq )
    NAME=${BASE%%_1_paired.fastq}
    echo "srun bowtie2 -p 3 -X 1500 --rg-id $NAME --rg SM:$NAME -D 20 -R 3 -N 1 --score-min L,0,-0.384 -L 20 -i S,1,0.50 -x $1 -1 ${NAME}_1_paired.fastq -2 ${NAME}_2_paired.fastq $
    echo "samtools view -bSq 10 ${NAME}.sam > ${NAME}_BTVS-UNSORTED.bam " >> TEMPBATCH.sbatch
    echo "srun samtools sort ${NAME}_BTVS-UNSORTED.bam ${NAME}_UNDEDUP" >> TEMPBATCH.sbatch
    echo "srun java -Xmx4g -jar /share/PI/spalumbi/programs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=${NAME}_UNDEDUP.bam OUTPUT=${NAME}.bam METRICS_FILE=${NAME}-metr$
    echo "srun samtools index ${NAME}.bam" >> TEMPBATCH.sbatch
    echo "rm ${NAME}.sam" >> TEMPBATCH.sbatch
    echo "rm ${NAME}_BTVS-UNSORTED.bam" >> TEMPBATCH.sbatch
    echo "rm ${NAME}_UNDEDUP.bam" >> TEMPBATCH.sbatch
    let COUNTER=COUNTER+1
    if [ $COUNTER -eq $CHUNK ]; then
    sbatch TEMPBATCH.sbatch
    COUNTER=0; fi
done
if [ $COUNTER -ne 0 ]; then
sbatch TEMPBATCH.sbatch; fi



