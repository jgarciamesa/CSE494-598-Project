#!/bin/bash
#SBATCH -n 2
#SBATCH -t 0-01:00:00

# simulate_illumina_reads.sh
# Simulate Illumina sequence reads using ART

mkdir -p ./results/simulate_reads

#$1 fasta file
#$2 prefix for output reads

art_output="$(art_illumina -ss HS25 -i $1 -o $2 -l 150 -f 50 -p -m 500 -s 10 -sp -sam -na)"

# grab fastq and sam output file names from output of art command
fq1name="$(echo "$art_output" | grep "1st reads" | sed 's/.*reads:[[:space:]]\(.*\.fq\).*/\1/')"
fq2name="$(echo "$art_output" | grep "2nd reads" | sed 's/.*reads:[[:space:]]\(.*\.fq\).*/\1/')"
samname="$(echo "$art_output" | grep "\.sam" | sed 's/[[:space:]]\(.*\.sam\)/\1/')"

mv "./$fq1name" "./results/simulate_reads/$fq1name"
mv "./$fq2name" "./results/simulate_reads/$fq2name"
mv "./$samname" "./results/simulate_reads/$samname"

echo "Illumina read results can be found at ./results/simulate_reads/"
