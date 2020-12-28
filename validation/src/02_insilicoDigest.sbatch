#!/bin/bash
#SBATCH --job-name insilicoDigest
#SBATCH --output /scratch/users/dromanovskaia/compEpi/validation_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --mem 16000
#SBATCH --cpus-per-task=4
#SBATCH --time 02:00:00
#SBATCH --array=20

conda activate crossSp

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

date
Rscript ${CODEBASE}/DNAmeth500species/validation/src/02_insilicoDigest.R $SLURM_ARRAY_TASK_ID
date
