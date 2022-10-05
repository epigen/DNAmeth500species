#!/bin/bash

#SBATCH --job-name meme
#SBATCH --output ../../log/06.3.2_%A_%a.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=1-192%100

ANALYSIS_DIR=../../results_analysis/
export SCRATCH_PATH=${LOGDIR}
export TMPDIR=${LOGDIR}/tmp

echo "$SLURM_ARRAY_TASK_ID"
SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' ${ANALYSIS_DIR}/06_inv/classes_with_inv.txt`
echo "$SP"

for k in 3 6 9 12
do
date
bash ${CODEBASE}/DNAmeth500species/src/06.3.2_kmer_count.sh $SP $k $ANALYSIS_DIR
date
done
