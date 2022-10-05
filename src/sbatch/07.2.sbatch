#!/bin/bash
#SBATCH --job-name JASPAR
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --output ../../log/07.2_%A_%a.log
#SBATCH --array=1-2

source ~/.bashrc
conda activate Zoo_R

echo "TEST"


export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp


echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ${CODEBASE}/DNAmeth500species/meta/species_list.txt`
echo "$SP"

date
<your env path here>/Rscript  ${CODEBASE}/DNAmeth500species/src/07.2_TF_freq.R ${SP}
date
