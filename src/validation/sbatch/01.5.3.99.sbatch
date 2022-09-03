#!/bin/bash
#SBATCH --job-name filterSeq
#SBATCH --output ../../log/01.5.3.99_%j.log
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1


#setup
source ~/.bashrc
conda activate Zoo_R

LOGDIR="/binfl/lv71484/droman/DNAmeth500species/log/mytmp"

export SCRATCH_PATH=$LOGDIR/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"

bash 01.5.3.99_refseq_bash.sh