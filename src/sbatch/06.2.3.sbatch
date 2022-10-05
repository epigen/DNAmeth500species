#!/bin/bash
#SBATCH --job-name Zoo_inv
#SBATCH --output ../../log/06.2.3_%j.log
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1

export SCRATCH_PATH=${LOGDIR}
export TMPDIR=${LOGDIR}

SP=${1}

date
<your env path here>/Rscript ${CODEBASE}/DNAmeth500species/src/06.2.3_verify_inv_subsampling.R $SP
date
