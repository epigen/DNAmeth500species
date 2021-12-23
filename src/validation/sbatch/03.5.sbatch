#!/bin/bash
#SBATCH --job-name kebabs_WGBS
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/05.1_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=1

source ~/.bashrc
conda activate compEpi

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

sp="DOG"

echo "$sp"

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/validation/03.5_RRBS_test_in_WGBS.R ${sp}
date
