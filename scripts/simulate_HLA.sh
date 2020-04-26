#!/bin/bash



# Simulate n individuals and output a fasta file containing the HLA alleles.
#  Fasta contains two alleles per gene per individual, using genes: A, B, C,
#  DQA1, DQB1, DRB1.

if [ $# -lt 1 ]
then
	echo 'Number of individuals needed. Exiting!'
	exit
fi
n=$1   # number of individuals to simulate

gene_a=https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/A_gen.fasta
gene_b=https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/B_gen.fasta
gene_c=https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/C_gen.fasta
gene_dqa1=https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/DQA1_gen.fasta
gene_dqb1=https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/DQB1_gen.fasta
gene_drb1=https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/DRB1_gen.fasta

genes=($gene_a $gene_b $gene_c $gene_dqa1 $gene_dqb1 $gene_drb1)

# Download allele sequences
mkdir -p ./data/fasta
for gene in ${genes[*]}
do
	stem=$(basename $gene)
	if [ ! -e data/fasta/$stem ]
	then
		echo 'Downloading' $stem
		wget --quiet -O ./data/fasta/$stem $gene
	fi

	# create list of sequence ID for each gene
	grep -e '>' data/fasta/${stem} | grep -o '^\S*' > data/fasta/id_${stem}

done


# initialize mixed sample
# echo -n > results/fasta/mixed_${n}.fasta
# simulate diploid alleles for each individual and merge
mkdir -p ./results/fasta
for i in $(seq 1 $n)
do
	for gene in ${genes[*]}
	do
		for j in $(seq 1 2)
		do
			stem=$(basename $gene)  									# get filename
			length=$(wc -l data/fasta/id_${stem} | grep -E -o '^[^[:space:]]*')
			rand=$((1 + RANDOM % $length))  							# random ID position
			randp1=$(( $rand + 1 ))    									# next position
			id=$(head -n $rand data/fasta/id_${stem} | tail -n 1)  		# ID
			bound=$(head -n $randp1 data/fasta/id_${stem} | tail -n 1)  # ID of next
			pos=$(grep -n $id data/fasta/${stem} | grep -o -E ^[[:digit:]]+)  # position
			posp1=$(grep -n $bound data/fasta/${stem} | grep -o -E ^[[:digit:]]+) # position of next
			pos2=$(($posp1 - 1))
			l=$(($pos2 - $pos + 1))	  									# length of sequence
			# concatenate sequence to mixed fasta
			head -n $pos2 data/fasta/${stem} | tail -n ${l} >> results/fasta/mixed_${n}.fasta
		done
	done
done

# save information for "known individual"
# echo -n > results/fasta/mixed_${n}.txt
echo 'HLA alleles' >> results/fasta/mixed_${n}.txt

grep 'A\*' results/fasta/mixed_${n}.fasta | head -n 2 >> results/fasta/mixed_${n}.txt
grep 'B\*' results/fasta/mixed_${n}.fasta | head -n 2 >> results/fasta/mixed_${n}.txt
grep 'C\*' results/fasta/mixed_${n}.fasta | head -n 2 >> results/fasta/mixed_${n}.txt
grep 'DQA1\*' results/fasta/mixed_${n}.fasta | head -n 2 >> results/fasta/mixed_${n}.txt
grep 'DQB1\*' results/fasta/mixed_${n}.fasta | head -n 2 >> results/fasta/mixed_${n}.txt
grep 'DRB1\*' results/fasta/mixed_${n}.fasta | head -n 2 >> results/fasta/mixed_${n}.txt
