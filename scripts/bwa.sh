#!/bin/bash
#SBATCH -n 4
#SBATCH -t 0-00:30:00

module load bwa/0.7.17
module load samtools/1.9

# bwa.bash
# Align files

if [ $# -lt 6 ]
then
    echo "Usage: bash bwa.sh ref_genome 1.fastq 2.fastq out.sam out.bam out.sorted.bam"
    exit
fi

#$1 reference genome
#$2 _1.fastq
#$3 _2.fastq
#$4 output.sam
#$5 output.bam
#$6 output.sorted.bam

bwa mem $1 $2 $3 > $4

# convert aligned file to bam
samtools view -S -b $4 > $5

# sort and index
samtools sort -o $6 $5
samtools index $6
