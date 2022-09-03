#!/bin/bash
#SBATCH --job-name Zoo_inv
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/06.2.3_%j.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1

export SCRATCH_PATH=/binfl/lv71484/droman/DNAmeth500species/log/mytmp
export TMPDIR=/binfl/lv71484/droman/tmp

SP=${1}

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/06.2.3_verify_inv_subsampling.R $SP
date
