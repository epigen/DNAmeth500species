#!/bin/bash
#SBATCH --job-name kebabs_WGBS
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/05.1_%A_%a.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=4,7

source ~/.bashrc
conda activate compEpi

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/03_WGBS/WGBS_prediction_selection.tsv"

sp=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' ${ANNOT_TABLE}`
genomeid=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $2}' ${ANNOT_TABLE}`
cov_min=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`
matched=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $3}' ${ANNOT_TABLE}`
echo "$sp"

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/validation/03.4_within_WGBS.R ${sp} ${genomeid} ${cov_min} ${matched}
date
