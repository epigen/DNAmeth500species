#!/bin/bash
#SBATCH --job-name kebabs_WGBS
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/03.5_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=6-580

source ~/.bashrc
conda activate compEpi

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

sp=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' $CODEBASE/DNAmeth500species/meta/species_list.txt`

echo "$sp"

date
/home/lv71484/droman/miniconda3/envs/Zoo_SVM_4.1/bin/Rscript ${CODEBASE}/DNAmeth500species/src/validation/03.5_RRBS_test_in_WGBS.R ${sp}
date
