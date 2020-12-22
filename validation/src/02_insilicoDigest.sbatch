#!/bin/bash
#SBATCH --job-name insilicoDigest
#SBATCH --output /scratch/users/dromanovskaia/compEpi/validation_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --mem 32000
#SBATCH --cpus-per-task=4
#SBATCH --time 08:00:00
#SBATCH --array=1-47%10

conda activate crossSp

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

date
Rscript ${CODEBASE}/DNAmeth500species/validation/src/02_insilicoDigest.R $SLURM_ARRAY_TASK_ID
date
