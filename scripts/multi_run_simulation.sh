#!/bin/bash
#SBATCH -n 1
#SBATCH -t 0-10:00:00

# multi_run_simulation.sh

if [ $# -lt 2 ]
then
    echo 'Usage: multi_run_simulations.sh number_of_individuals number_of_runs'
    exit
fi

n=$1
num_runs=$2

for run_num in $(seq $num_runs) 
do
    sbatch ./scripts/run_simulation.sh $n
done
