
#!/bin/bash
#USAGE: give ref, chunksize, and a glob of reads
#for 50bp SE reads, this script allows 2 mismatches or 1 mismatch/1 indel for every 50 bp. see --score-min section
#EXAMPLE (from within directory with fastqs): bash batch-bowtie2-fq-single.sh $PI_HOME/ref/ahy_nosym 1 *.txt.gz
CHUNK=$2
COUNTER=0
FQ="${@:3}"
echo $FQ
for i in $FQ; do
    if [ $COUNTER -eq 0 ]; then
    echo -e "#!/bin/bash\n#SBATCH -p owners\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=3\n#SBATCH -t 6:00:00\n#SBATCH --mem 12000" > TEMPBATCH.sbatch; fi
    BASE=$(basename $(basename $(basename $( basename $i .gz) .txt) .fq) .fastq)
    echo submitting $BASE
    echo "srun bowtie2 -p 3 --rg-id $BASE --rg SM:$BASE -D 20 -R 3 -N 1 --score-min L,0,-0.32 -L 20 -i S,1,0.50 -x $1 -U $i > temp$BASE.sam" >> TEMPBATCH.sbatch
    echo "samtools view -bSq 20 temp$BASE.sam > ${BASE}-UNSORTED.bam " >> TEMPBATCH.sbatch
    echo "rm temp$BASE.sam" >> TEMPBATCH.sbatch
    echo "srun samtools sort ${BASE}-UNSORTED.bam ${BASE}" >> TEMPBATCH.sbatch
    echo "srun samtools index ${BASE}.bam" >> TEMPBATCH.sbatch
    echo "rm ${BASE}-UNSORTED.bam" >> TEMPBATCH.sbatch
    let COUNTER=COUNTER+1
    if [ $COUNTER -eq $CHUNK ]; then
    sbatch TEMPBATCH.sbatch
    COUNTER=0; fi
done
if [ $COUNTER -ne 0 ]; then
sbatch TEMPBATCH.sbatch; fi
