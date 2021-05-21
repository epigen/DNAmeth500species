#!/bin/bash
#SBATCH --job-name kebabs_WGBS
#SBATCH --output ../../log/WGBS_kebabs_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --qos shortq
#SBATCH --mem 128G
#SBATCH --cpus-per-task=8
#SBATCH --array=6,7

source ~/.bashrc
conda activate compEpi

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE=${1}

SPECIES=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' ${ANNOT_TABLE}`
echo ${SPECIES}
ASSEMBLY=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $3}' ${ANNOT_TABLE}`
ABBR=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`
COV=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $7}' ${ANNOT_TABLE}`
echo "$SPECIES"
date
Rscript 05.3_WGBS.R ${SPECIES} ${ASSEMBLY} ${ABBR} ${COV}
date
