mkdir  $LOGDIR/compEpi/

#!/bin/bash

#SBATCH --job-name compEpiQC
#SBATCH --output $LOGDIR/compEpi/00.21_%j\A_%a.txt
#SBATCH --partition shortq
#SBATCH --mem 80000
#SBATCH --cpus-per-task=16
#SBATCH --time 8:00:00
#SBATCH --array=4-582

module load R/3.5.1

export SCRATCH_PATH=/scratch/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`

echo "$SP"

date
Rscript ${CODEBASE}/compEpi/src/00.21_filter_for_QC.R $SP
date
