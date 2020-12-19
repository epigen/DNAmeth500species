#!/bin/bash

#SBATCH --job-name insilicoDigest
#SBATCH --output out.txt
#SBATCH --partition shortq
#SBATCH --mem 8000
#SBATCH --cpus-per-task=4
#SBATCH --time 08:00:00


conda activate crossSp

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

date
Rscript ${CODEBASE}/DNAmeth500species/validation/insilicoDigest.R Mouse
date
