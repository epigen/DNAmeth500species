#!/bin/bash
#SBATCH --job-name concat
#SBATCH --output ../../log/concat_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --qos shortq
#SBATCH --mem 16000
#SBATCH --time 08:00:00
#SBATCH --array=1-5

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE="../meta/genomes_to_concatinate.txt"

db=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $2}' ${ANNOT_TABLE}`
genome_id=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`

echo ${db}

date
/nobackup/lab_bock/users/dromanovskaia/miniconda3/envs/crossSp/bin/Rscript 04_concat_genome.R ${db} ${genome_id}
date
