#!/bin/bash
#SBATCH --job-name SRA
#SBATCH --output ../../log/SRA_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --qos shortq
#SBATCH --array=3

source ~/.bashrc
conda activate crossSp

ANNOT_TABLE=${1}

SRR_id=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`

date
fasterq-dump --outdir $2 ${SRR_id}
date
