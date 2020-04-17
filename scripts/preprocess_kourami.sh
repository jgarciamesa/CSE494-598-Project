#!/bin/bash
#SBATCH -n 2
#SBATCH -t 0-00:30:00

module load bamutil/1.0.14
module load samtools/1.9
module load bwa/0.7.17

#$1 output prefix
#$2 input bam file

./tools/kourami-0.9.6/scripts/alignAndExtract_hs38DH_NoAlt.sh $1 $2
