#!/bin/bash
## 1 - species
## 2 - order of the mm (default 4)
## 3 analysis_dir 
analysis_dir=$3
PATH_TO_FASTA="${analysis_dir}/00.2_QC/filtered_sequences/${1}/toSelf_filtered_0.08mm_final_filtered.fa"

PATH_TO_SAVE_DIR="${analysis_dir}/99.4_kmer_count_filtered/${1}_${2}"

mkdir $PATH_TO_SAVE_DIR
#module load meme/4.10.2 

date
fasta-get-markov -m ${2} -norc $PATH_TO_FASTA ${PATH_TO_SAVE_DIR}/count.txt
date


##splitting the file into several
cd $PATH_TO_SAVE_DIR
csplit -f kmer count.txt -z '/# order/' '{*}'
