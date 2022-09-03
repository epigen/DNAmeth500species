#!/bin/bash
#SBATCH --job-name Zoo_inv
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/06.2.1_%j.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1

source ~/.bashrc
conda activate Zoo_SVM_4.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=${1}

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/06.2.1_verify_inv_tissue.R $SP
date
