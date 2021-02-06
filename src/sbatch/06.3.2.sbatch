#!/bin/bash
#SBATCH --job-name kmercount
#SBATCH --output $LOGDIR/compEpi/99.4_%j.txt
#SBATCH --partition shortq
#SBATCH --mem 80000
#SBATCH --cpus-per-task=16
#SBATCH --time 8:00:00
#SBATCH --array=1-583

ANALYSIS_DIR=<data_dir>/analysis_dir/
export SCRATCH_PATH=$LOGDIR/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp
echo "$SLURM_ARRAY_TASK_ID"
SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' ../../meta/classes_with_inv.txt`
echo "$SP"
for k in 3 6 9 12
do
date
${CODEBASE}/compEpi/src/06.3.2_kmer_count.sh $SP $k $ANALYSIS_DIR
date
done
