#!/bin/bash
#SBATCH --job-name diffMeth
#SBATCH --output $LOGDIR/compEpi/03.98_%j.txt
#SBATCH --partition mediumq
#SBATCH --mem 8000
#SBATCH --cpus-per-task=16
#SBATCH --time 8:00:00
#SBATCH --array=1-174

export SCRATCH_PATH=$LOGDIR/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"
SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ${CODEBASE}/compEpi/meta/species_list_TF_analysis.txt`
echo "$SP"

date
${CODEBASE}/DNAmeth500species/src/07.4.1_run_diffMeth.sh $SP
date
