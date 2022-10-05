#!/bin/bash
#SBATCH --job-name kebabs
#SBATCH --output ../../log/05.1.1_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=7

source ~/.bash_profile

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRAT=/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' $CODEBASE/DNAmeth500species/meta/species_list.txt`

echo "$SP"


date
<your env path here>/Rscript ${CODEBASE}/DNAmeth500species/src/05.1.1_predict_meth_kebabs_alt.R $SP
date
