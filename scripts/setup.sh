#!/bin/bash                                                                     
#SBATCH -n 2                                                                    
#SBATCH -t 0-06:00:00                                                           
                                                                                
module load trimmomatic/0.33                                                    
module load bwa/0.7.17                                                          
module load samtools/1.9       

################################################################################
# Set up project folder                                                        #
################################################################################

mkdir -p data/{fastq,ref_genome}
mkdir -p results/{sam,bam,simulated_reads,fasta}
mkdir tools

exome_hla=https://s3.amazonaws.com/omixon-publication/hapmap_hla/HapMap_1KG_HLA_suppl_filtered_reads.tgz
hlascan=https://github.com/SyntekabioTools/HLAscan/releases/download/v2.1.4/hla_scan_r_v2.1.4

# download whole-exome sequencing reads (fastq)
cd data/fastq
echo "Downloading whole-exome sequencing reads" 
wget $exome_hla
echo "Decompressing whole-exome sequencing reads"
tar xzfv $(basename $exome_hla)
cd ../..

cd tools
echo "Downloading HLAscan"
wget $hlascan
cd ../..


# if the human reference genome is not downloaded do so                         
if [ ! -f "data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" ]        
then                                                                            
    echo "Downloading refernce genome" > progress.log                           
    cd data/ref_genome/                                                         
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
    cd ../..                                                                    
fi                                                                              
                                                                                
ref=./data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa                
                                                                                
# if the human reference genome is indexed do so                                
if [ ! -f "data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt" ]    
then                                                                            
    echo "Indexing reference genome" >> progress.log                            
    bwa index $ref                                                              
fi  

