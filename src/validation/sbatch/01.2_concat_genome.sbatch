#!/bin/bash
#SBATCH --job-name concat
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/concat_%A_%a.txt
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=1

source ~/.bashrc

export SCRATCH_PATH=/binfl/lv71484/droman/DNAmeth500species/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

#ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/01_crossMapping/genomes_ucsc_todownload.tsv"
ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/02_insilico_digest/genomes_ucsc_todownload_for_insilico.tsv"
db=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' ${ANNOT_TABLE}`
#genome_id=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`

echo ${db}

date

/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ../01.2_concat_genome.R AmexG_v3.0.0
#/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ../01.2_concat_genome.R ${db}
date
