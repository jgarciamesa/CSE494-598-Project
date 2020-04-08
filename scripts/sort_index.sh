#!/bin/bash
#SBATCH -n 2
#SBATCH -t 0-04:00:00

module load samtools/1.9

bamfile=results/bam/NA18522_SRR107025_SRP004074.bam
stem=$(basename $bamfile .bam)

samtools sort -o "${stem}.sorted.bam" ${bamfile}

samtools index "${stem}.sorted.bam"

