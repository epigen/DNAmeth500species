#!/bin/bash
#SBATCH --job-name Zoo_inv
#SBATCH --output ../../log/06.2.2_%j.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1

source ~/.bashrc
conda activate Zoo_SVM_4.1

export SCRATCH_PATH=${LOGDIR}
export TMPDIR=${LOGDIR}/tmp


SP=${1}

date
<your env path here>/Rscript ${CODEBASE}/DNAmeth500species/src/06.2.2_verify_inv_sample.R $SP
date
