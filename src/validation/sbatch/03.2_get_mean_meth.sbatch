#!/bin/bash
#SBATCH --job-name WGBS
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/WGBS_mean_%A_%a.txt
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1

source ~/.bashrc

export SCRATCH_PATH=/binfl/lv71484/droman/DNAmeth500species/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

#ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/01_crossMapping/genomes_ucsc_todownload.tsv"

date

/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript $CODEBASE/DNAmeth500species/src/validation/03.2_process_WGBS_data.R  
date
