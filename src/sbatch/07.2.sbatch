#!/bin/bash
#SBATCH --job-name TFBS_det
#SBATCH --output ${LOGDIR}/compEpi/03.2_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --mem 8000
#SBATCH --cpus-per-task=4
#SBATCH --time 2:00:00
#SBATCH --array=4-287%20

module load R/3.5.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"
SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ${CODEBASE}/DNAmeth500species/meta/species_list.txt`
echo "$SP"

date
${CODEBASE}/DNAmeth500species/src/03.2_TF_freq.R ${SP}
date
