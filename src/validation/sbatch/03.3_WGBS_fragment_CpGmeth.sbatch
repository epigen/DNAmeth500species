#!/bin/bash
#SBATCH --job-name WGBS
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/WGBS_%A_%a.txt
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=3

source ~/.bashrc
conda activate Zoo_R

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/03_WGBS/WGBS_prediction_selection.tsv"
genome_id=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $3}' ${ANNOT_TABLE}`

date
/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ${CODEBASE}/DNAmeth500species/src/validation/03.3_WGBS_fragment_CpGmeth.R ${genome_id}
date
