#!/bin/bash
#wrapper for hla_tally.py

module load python/3.7.1
module load r/latest

python ./scripts/hla_tally.py
Rscript --vanilla ./scripts/mtdna_tally.R
