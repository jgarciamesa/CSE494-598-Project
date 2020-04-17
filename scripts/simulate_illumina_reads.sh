#!/bin/bash
#SBATCH -n 2
#SBATCH -t 0-01:00:00


# simulate_illumina_reads.sh
# Simulate Illumina sequence reads using ART

mkdir -p ../results/simulate_reads

#$1 fasta file
#$2 prefix for output reads

tools/art_illumina -ss HS25 -i $1 -o $2 -l 150 -f 50 -p -m 500 -s 10 -sp -sam -na
