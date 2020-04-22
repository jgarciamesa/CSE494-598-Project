#!/bin/bash
#SBATCH -n 10
#SBATCH -t 0-05:00:00

module load bwa/0.7.17
module load samtools/1.9
module load bamUtil/1.0.14

# run.sh
# Usage: bash run.sh n_individuals

if [ $# -lt 1 ]
then
	echo 'Argument with number of individuals is needed. Exiting!'
	exit
fi

n=$1

# create folders
mkdir -p data/fasta
mkdir -p results/{bam/fasta/sam/simulated_reads}
mkdir tools

# install ART and Kuorami inside the folder tools
#https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
#https://github.com/Kingsford-Group/kourami

bash simulate_individuals.sh $n

bash simulated_reads.sh results/fasta/mixed_${n}.fasta results/simulated_reads/mixed_{n}

bash bwa.sh
