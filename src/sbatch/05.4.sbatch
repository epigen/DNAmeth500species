#!/bin/bash
#SBATCH --job-name compEpi
#SBATCH --output /scratch/users/dromanovskaia/compEpi/02.2_%j.txt
#SBATCH --partition longq
#SBATCH --mem 240000
#SBATCH --cpus-per-task=16
#SBATCH --time 24:00:00
#SBATCH --array=1-6

module load R/3.5.1

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`
echo "species: ${SP}"

date
$CODEBASE/compEpi/src/05.2_predict_crossSp.R ${SP}
date
