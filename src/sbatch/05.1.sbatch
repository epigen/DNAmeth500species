#!/bin/bash
#SBATCH --job-name kebabs
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/05.1_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=1-580%50

source ~/.bashrc
conda activate Zoo_SVM_4.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' $CODEBASE/DNAmeth500species/meta/species_list.txt`

echo "$SP"

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/05.1_predict_meth_kebabs.R $SP
date
