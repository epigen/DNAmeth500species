#!/bin/bash
#SBATCH --job-name concat
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/crossMapping_%A_%a.txt
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=42,50,967

export SCRATCH_PATH=/binfl/lv71484/droman/DNAmeth500species/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/01_crossMapping/01.4_crossMapping_stats/cross_mapping_stats.tsv"

species=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $9}' ${ANNOT_TABLE}`
genome_id=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $10}' ${ANNOT_TABLE}`
echo "$species"

source ~/.bashrc
conda activate Zoo_R

date
#/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript 01.2_concat_genome.R ${db}
/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ${CODEBASE}/DNAmeth500species/src/validation/01.5_cross-mapping_analysis.R ${species} ${genome_id}
date

