#!/bin/bash
#SBATCH --job-name TFBS_det
#SBATCH --output $LOGDIR/compEpi/03.1_%A_%a.txt
#SBATCH --partition mediumq
#SBATCH --mem 240000
#SBATCH --cpus-per-task=16
#SBATCH --time 24:00:00
#SBATCH --array=404-583%20

module load R/3.5.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"
SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ${CODEBASE}/DNAmeth500species/meta/species_list.txt`
echo "$SP"

date
${CODEBASE}/DNAmeth500species/src/07.1_TFBS_detection.R ${SP}
date
