#!/bin/bash
#SBATCH --job-name SRA
#SBATCH --output /nobackup/lab_bock/projects/DNAmeth500species/log/SRA_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --qos shortq
#SBATCH --array=3

source ~/.bash_profile
conda activate compEpi

#ANNOT_TABLE=${1}

#SRR_id=`awk -F";" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`

mkdir ${1}

SRR_id="SRR2070790"

date
fastq-dump --outdir ${1} ${SRR_id}
date
