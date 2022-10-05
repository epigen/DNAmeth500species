#!/bin/bash
#SBATCH --job-name concat
#SBATCH --output ../../log/concat_%A_%a.txt
#SBATCH --partition vsc3plus_0256
#SBATCH --qos vsc3plus_0256
#SBATCH -N 1
#SBATCH --array=1

source ~/.bashrc

export SCRATCH_PATH=$P{LOGDIR}
export TMPDIR=$SCRATCH_PATH/Rtmp

ANNOT_TABLE="../../results_analysis/validation/02_insilico_digest/genomes_ucsc_todownload_for_insilico.tsv"
db=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $1}' ${ANNOT_TABLE}`
#genome_id=`awk -F"\t" 'NR=='$SLURM_ARRAY_TASK_ID'{print $5}' ${ANNOT_TABLE}`

echo ${db}

date

<your env path here>/Rscript ../01.2_concat_genome.R AmexG_v3.0.0

date
