#!/bin/bash
#SBATCH --job-name kmercount
#SBATCH --output ../../log/01.5.2_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH --array=183,472
#SBATCH -N 1

ANALYSIS_DIR=/binfl/lv71484/droman/DNAmeth500species/results_analysis
export SCRATCH_PATH=/binfl/lv71484/droman/DNAmeth500species/log/mytmp
export TMPDIR=/binfl/lv71484/droman/tmp

echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`

echo "$SP"
k=10

date
bash /home/lv71484/droman/DNAmeth500species/src/01.5.2_kmer_count.sh $SP $k $ANALYSIS_DIR
date
