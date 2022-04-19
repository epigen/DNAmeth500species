#!/bin/bash
## 1 - species
## 2 - order of the mm (default 4)
## 3 analysis_dir 

##this script is tun on any set of sequences provided and generates a kmer count matrix


resource_dir=$3

PATH_TO_FASTA="${resource_dir}/${1}/${1}_concatinated.fa"
PATH_TO_SAVE_DIR="${resource_dir}/${1}/"

mkdir $PATH_TO_SAVE_DIR
#module load meme/4.10.2 

date
fasta-get-markov -m ${2} -norc $PATH_TO_FASTA ${PATH_TO_SAVE_DIR}/count.txt
date


##splitting the file into several
cd $PATH_TO_SAVE_DIR
csplit -f kmer count.txt -z '/# order/' '{*}'
