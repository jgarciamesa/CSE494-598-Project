#!/bin/bash
#SBATCH -n 2
#SBATCH -t 0-01:00:00

fastq1=data/fastq/SRR062634_1.filt.fastq
fastq2=data/fastq/SRR062634_2.filt.fastq
stem=$(basename $fastq1 _1.filt.fastq)

# if fast1 is not decompressed, do so
if [ ! -f $fastq1 ]
then
    gunzip ${fastq1}.gz
fi

# if fastq2 is not decompressed, do so
if [ ! -f $fastq2 ]
then
    gunzip ${fastq2}.gz
fi

# run hlascan
./tools/hlascan/hla_scan_r_v2.1.4 -l $fastq1 -r $fastq2 -d tools/hlascan/db/HLA-ALL.IMGT -v 38 -g HLA-A -t 4 > hlascan_${stem}_HLAA.log
