#!/bin/bash
#SBATCH --job-name kmercount
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/01.5.2_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH --array=1-76
#SBATCH -N 1

ANALYSIS_DIR=/binfl/lv71484/droman/DNAmeth500species/resources/reference_genomes

export SCRATCH_PATH=/binfl/lv71484/droman/DNAmeth500species/log/mytmp
export TMPDIR=/binfl/lv71484/droman/tmp

echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' /binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/02_insilico_digest/genomes_to_run.tsv`

echo "$SP"
k=1

date
bash /home/lv71484/droman/DNAmeth500species/src/validation/01.6.1_kmer_count.sh $SP $k $ANALYSIS_DIR
date
