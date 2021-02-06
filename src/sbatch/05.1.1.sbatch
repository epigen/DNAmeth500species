#!/bin/bash
#SBATCH --job-name kebabs
#SBATCH --output ${LOGDIR}/compEpi/02.13_%A_%a.txt
#SBATCH --partition longq
#SBATCH --mem 160000
#SBATCH --cpus-per-task=16
#SBATCH --time 24:00:00
#SBATCH --array=1-583%15

module load R/3.5.1

export SCRATCH_PATH=/scratch/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' $CODEBASE/compEpi/meta/species_list.txt`

echo "$SP"

date
Rscript ${CODEBASE}/compEpi/src/05.1.1_predict_meth_kebabs_alt.R $SP
date
