#!/bin/bash
#SBATCH --job-name AME
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --output ../../log/07.4.2_%A_%a.log
#SBATCH --array=1-2

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp


echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ${CODEBASE}/DNAmeth500species/meta/species_list_TF_analysis.txt`
echo "$SP"

date
bash ${CODEBASE}/DNAmeth500species/src/07.4.2_ame_run.sh $SP
date
