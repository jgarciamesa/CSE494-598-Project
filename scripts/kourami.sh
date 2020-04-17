#!/bin/bash
#SBATCH -n 10
#SBATCH -t 0-01:00:00

module load bamutil/1.0.14

#$1 output prefix
#$2 input BAM preprocessed for kourami

java -jar tools/kourami-0.9.6/build/Kourami.jar -o $1 $2 -d tools/kourami-0.9.6/db/
