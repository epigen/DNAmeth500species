#!/bin/bash
## 1 - species
## 2 - order of the mm (default 4)
analysis_dir=$3

PATH_TO_FASTA="${analysis_dir}/05_predict_meth/05.1_within_species/screen/${1}/sequences/"

PATH_TO_SAVE_DIR="${analysis_dir}/06_inv/06.3_kmer_freq/${1}/"

#mkdir $PATH_TO_SAVE_DIR
#module load meme/4.10.2 

date
#fasta-get-markov -m ${2} -norc #$PATH_TO_FASTA/${1}_high.fa #${PATH_TO_SAVE_DIR}/high.txt
#fasta-get-markov -m ${2} -norc #$PATH_TO_FASTA/${1}_low.fa #${PATH_TO_SAVE_DIR}/low.txt
date


##splitting the file into several
cd $PATH_TO_SAVE_DIR
csplit -f kmer_high high.txt -z '/# order/' '{*}' 
csplit -f kmer_low low.txt -z '/# order/' '{*}' 
