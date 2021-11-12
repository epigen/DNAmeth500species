#!/bin/bash
#SBATCH --job-name insilicoDigest
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/insilicoRRBS_%A_%a.txt
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=9,66,70

source ~/.bashrc

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE="/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/02_insilico_digest/genomes_to_run.tsv"
genome_id=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' ${ANNOT_TABLE}`


date
/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ${CODEBASE}/DNAmeth500species/src/validation/02.2_insilicoDigest.R ${genome_id}
date
