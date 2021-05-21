#!/bin/bash
#SBATCH --job-name WGBS
#SBATCH --output ../../log/WGBS_fragments_%A_%a.txt
#SBATCH --partition shortq
#SBATCH --qos shortq
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=7

source ~/.bashrc
conda activate compEpi

export SCRATCH_PATH=/nobackup/lab_bock/users/dromanovskaia/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

date
Rscript 05.2_WGBS_fragment_CpGmeth.R $SLURM_ARRAY_TASK_ID
date
