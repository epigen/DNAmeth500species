#!/bin/bash
#SBATCH --job-name CpG
#SBATCH --output ${LOGDIR}/compEpi/compEpi_01.8_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --mem 80000
#SBATCH --cpus-per-task=16
#SBATCH --time 8:00:00
#SBATCH --array=1-583
module load R/3.5.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`

date
Rscript ${CODEBASE}/DNAmeth500species/src/01.8_CpG_islands_per_ref.R ${SP}
date
