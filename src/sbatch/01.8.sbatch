#!/bin/bash
#SBATCH --job-name CpGislands
#SBATCH --output ../../log/01.8_%A_%a.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --array=1-581

#setup
source ~/.bashrc
conda activate Zoo_R

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/species_list.txt`

date
Rscript ${CODEBASE}/DNAmeth500species/src/01.8_CpG_islands_per_ref.R ${SP}
date
