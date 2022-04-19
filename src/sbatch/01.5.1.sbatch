#!/bin/bash
#SBATCH --job-name filterSequnces
#SBATCH --output ../../log/01.5.1_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=183,472

#setup
source ~/.bashrc
conda activate Zoo_R

export SCRATCH_PATH=$LOGDIR/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`

echo "$SP"

date
Rscript ${CODEBASE}/DNAmeth500species/src/01.5.1_prefilter_for_kmer_count.R $SP
date
