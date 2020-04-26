#!/bin/bash
#SBATCH -n 10
#SBATCH -t 0-05:00:00

module load bwa/0.7.17
module load samtools/1.9
module load bcftools/1.9
module load bamutil/1.0.14
module load r/latest
module load maven/3.3.9

# run.sh
# Usage: bash run.sh n_individuals

if [ $# -lt 1 ]
then
	echo 'Argument with number of individuals is needed. Exiting!'
	exit
fi

n=$1

################################################################################
# Setup                                                                        #
################################################################################

# install ART, Kourami, and haplogrep-cmd inside the folder tools
#https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
#https://github.com/Kingsford-Group/kourami - use hs38NoAltDH
#https://github.com/seppinho/haplogrep-cmd

# Setup ART if not downloaded
if [ ! -f tools/art_illumina ]
then
	echo 'Downloading ART'
	cd tools/
	wget -nv https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
	tar zxvf artbinmountrainier2016.06.05linux64.tgz
	rm artbinmountrainier2016.06.05linux64.tgz
	cp art_bin_MountRainier/art_illumina .
	cd ..
fi

# Download haplogrep jar executable if it doesn't exist inside tools
if [ ! -f tools/haplogrep-2.1.25.jar ]
then
	echo 'Downloading haplogrep'
	cd tools/
	wget https://github.com/seppinho/haplogrep-cmd/releases/download/2.1.25/haplogrep-2.1.25.jar
	cd ..
fi

# Download Kourami if not found
if [ ! -f tools/kourami-0.9.6/build/Kourami.jar ]
then
    echo 'Downloading Kourami'
    cd tools/
    rm -rf ./kourami* # clean slate
    wget https://github.com/Kingsford-Group/kourami/releases/download/v0.9.6/kourami-0.9.6_bin.zip
    unzip kourami-0.9.6_bin.zip
    cd kourami-0.9.6
    mvn install
    ./scripts/download_panel.sh
    ./scripts/download_grch38.sh hs38NoAltDH
    echo 'Indexing; this could take a while...'
    bwa index ./resources/hs38NoAltDH.fa
    cd ..
fi

# if mtDNA_dataset is not been created, do so
if [ ! -f data/mtDNA/mtDNA_dataset.tsv ]
then
	echo 'Creating mtDNA database'
	Rscript --vanilla --quiet scripts/create_mtDNA_database.R
fi

################################################################################
# Simulate mtDNA haplogroups and HLA alleles                                   #
################################################################################

# mtDNA haplogroups
echo 'Simulating mtDNA haplogroups'
Rscript --vanilla --quiet scripts/simulate_mtDNA.R $n

# HLA alleles
echo 'Simulating HLA alleles'
bash scripts/simulate_HLA.sh $n

################################################################################
# Simulate Illumina reads
################################################################################

echo 'Simulating illumina reads with ART'
art_output="$(tools/art_illumina -ss HS25 -i results/fasta/mixed_${n}.fasta -o mixed_${n} -l 150 -f 50 -p -m 500 -s 10 -sp -sam -na)"

fq1name="$(echo "$art_output" | grep "1st reads" | sed 's/.*reads:[[:space:]]\(.*\.fq\).*/\1/')"
fq2name="$(echo "$art_output" | grep "2nd reads" | sed 's/.*reads:[[:space:]]\(.*\.fq\).*/\1/')"
samname="$(echo "$art_output" | grep -o "[[:alnum:]_]*\.sam")"

mv "./$fq1name" "./results/simulate_reads/$fq1name"
mv "./$fq2name" "./results/simulate_reads/$fq2name"
mv -t "./results/simulate_reads/" ./$fq1name ./$fq2name ./$samname

echo "Completed; Illumina read results can be found at ./results/simulate_reads/"

################################################################################
# Align fastq to human reference genome , sort, and index                      #
################################################################################

# get path to human reference genome
#  (should have been downloaded and indexed when seting up Kourami)
ref_genome=$(find -name hs38NoAltDH.fa)

# align and sort
echo 'Aligning and sorting fastq reads to human reference genome for HLA typing'
bwa mem $ref_genome results/simulate_reads/${fq1name} results/simulate_reads/${fq2name} |\
samtools view -S -b | samtools sort -o results/bam/mixed_${n}.bam

# index
echo 'Indexing bam file'
samtools index results/bam/mixed_${n}.bam

################################################################################
# HLA typing with Kourami                                                      #
################################################################################

# preprocess bam file
echo 'Preprocessing bam file for kourami'
./tools/kourami-0.9.6/scripts/alignAndExtract_hs38DH_NoAlt.sh results/bam/mixed_${n} results/bam/mixed_${n}.bam

# run Kourami
echo 'HLA typing with Kourami'
java -jar tools/kourami-0.9.6/build/Kourami.jar -o results/kourami_${n}.txt results/bam/mixed_${n}_on_KouramiPanel.bam -d tools/kourami-0.9.6/db/

################################################################################
# Align fastq to reference sequence (rCRS) of Human mtDNA & variant calling    #
################################################################################

if [ ! -f data/mtDNA/mtDNA/mtDNA_reference.fasta.bwt ]
then
	bwa index data/mtDNA/mtDNA_reference.fasta
fi

# align and sort
echo 'Aligning and sorting fastq reads to human mtDNA reference sequence for mtDNA typing'
bwa mem data/mtDNA/mtDNA_reference.fasta results/simulate_reads/${fq1name} results/simulate_reads/${fq2name} |\
samtools view -S -b | samtools sort -o results/bam/mixed_${n}_mtDNA.bam

# index
echo 'Indexing bam file'
samtools index results/bam/mixed_${n}_mtDNA.bam

# variant calling
echo 'Variant calling'
bcftools mpileup -O v -f data/mtDNA/mtDNA_reference.fasta results/bam/mixed_${n}_mtDNA.bam | \
bcftools call --ploidy 1 -m -v -o results/vcf/mixed_${n}.vcf

################################################################################
# mtDNA typing with haplogrep-cmd
################################################################################

echo 'mtDNA typing with haplogrep'
java -jar tools/haplogrep-2.1.25.jar classify --in results/vcf/mixed_${n}.vcf --format vcf --out results/haplogrep_${n}.txt --extend-report


###############################################################################
# Make simulation folder for new data
###############################################################################

cd simulations
declare -i numFiles 
numFiles=$(ls -1 | wc -l)
((numFiles = numFiles + 1))
mkdir "simulation_$numFiles"
cd ../

###############################################################################
# Move results data into new folder
###############################################################################

mv "./results" "simulations/simulation_$numFiles"

###############################################################################
# Remake moved folders
###############################################################################

mkdir results
mkdir results/vcf
mkdir results/simulate_reads
mkdir results/bam
mkdir results/fasta
