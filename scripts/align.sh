#!/bin/bash
#SBATCH -n 2
#SBATCH -t 0-01:00:00

module load bwa/0.7.17

fastq1=data/exome/NA18522_SRR107025_SRP004074_1.hla_filt.fastq
fastq2=data/exome/NA18522_SRR107025_SRP004074_1.hla_filt.fastq
stem=$(basename $fastq1 _1.hla_filt.fastq)

bwa mem data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa $fastq1 $fastq2 > results/bam/${stem}.bam
