#!/bin/bash
#SBATCH --job-name dedRefCount
#SBATCH --output ../../log/01.4_dedRefCount.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1

source ~/.bashrc
conda activate Zoo_R

export SCRATCH_PATH=$LOGDIR/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp


date
Rscript ${CODEBASE}/DNAmeth500species/src/01.4_dedRef_CpG_count.R
date
