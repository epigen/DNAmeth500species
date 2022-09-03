#!/bin/bash
#SBATCH --job-name JASPAR
#SBATCH --partition vsc3plus_0064
#SBATCH --qos vsc3plus_0064
#SBATCH -N 1
#SBATCH --output /binfl/lv71484/droman/DNAmeth500species/log/07.1_%A_%a.log
#SBATCH --array=1-140

source ~/.bashrc
conda activate Zoo_R


LOGDIR="/binfl/lv71484/droman/DNAmeth500species/log/mytmp"

export SCRATCH_PATH=${LOGDIR}/mytmp
export TMPDIR=$SCRATCH_PATH/Rtmp

echo "$SLURM_ARRAY_TASK_ID"

SP=`awk 'NR=='$SLURM_ARRAY_TASK_ID'{print $0}' /binfl/lv71484/droman/DNAmeth500species/results_analysis/07_motifAnalysis/missing_species.txt`

#${CODEBASE}/DNAmeth500species/meta/species_list.txt`
echo "$SP"


pwd
date
#find motifs

/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ${CODEBASE}/DNAmeth500species/src/07.1_TFBS_detection.R ${SP}
##summarize
/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript  ${CODEBASE}/DNAmeth500species/src/07.2_TF_freq.R ${SP}

cd ${CODEBASE}/DNAmeth500species/results_analysis/07_motifAnalysis/07.1_TFBS_detection_2020
tar -zcvf ${SP}.tar.gz ${SP}
rm -r ${SP}
date
