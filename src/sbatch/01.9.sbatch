#!/bin/bash
#SBATCH --job-name coverage
#SBATCH --output ../../log/coverage_%A_%a.txt
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=183

export SCRATCH_PATH=<your path here>
export TMPDIR=$SCRATCH_PATH/Rtmp

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`

echo "$SP"

source ~/.bashrc
conda activate Zoo_R

date
/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ${CODEBASE}/DNAmeth500species/src/01.9_coverage_analysis.R ${SP} 
date

