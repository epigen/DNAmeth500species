#!/bin/bash
#SBATCH --job-name kebab
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/05.4_%A_%a.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=183,472

source ~/.bashrc
conda activate Zoo_SVM_4.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' $CODEBASE/DNAmeth500species/meta/species_list.txt`

echo "$SP"

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/05.4_predict_crossSp.R ${SP}
date
