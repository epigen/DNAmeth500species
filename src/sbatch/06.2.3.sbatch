#!/bin/bash
#SBATCH --job-name kebabs
#SBATCH --output ${LOGDIR}/compEpi/02.43_%j.txt
#SBATCH --partition longq
#SBATCH --mem 160000
#SBATCH --cpus-per-task=16
#SBATCH --time 24:00:00

module load R/3.5.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

SP=${1}

date
Rscript ${CODEBASE}/compEpi/src/06.1.3_verify_inv_subsampling.R $SP
date
