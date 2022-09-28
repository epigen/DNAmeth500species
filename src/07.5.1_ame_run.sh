#!/bin/bash

## $1 - the species to run the analysis on
DATA_DIR="/binfl/lv71484/droman/DNAmeth500species/"
MOTIF_FILE=${DATA_DIR}"resources/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms.meme"
FASTA_DIR=$DATA_DIR"/results_analysis/07_motifAnalysis/07.4_diffMeth_additional_run/${1}/motifAnalysis/"
OUTPUT_DIR=$DATA_DIR"/results_analysis/07_motifAnalysis/07.5_AME/screen/${1}"
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR"/Heart_vs_Liver_top" 
mkdir "$OUTPUT_DIR/Liver_vs_Heart_top" 

echo "running Heart vs Liver ${1}"

ame --o $OUTPUT_DIR/Heart_vs_Liver_top --control $FASTA_DIR/ded_top500_cov2_Liver.fa $FASTA_DIR/ded_top500_cov2_Heart.fa $MOTIF_FILE

echo "running Liver vs Heart ${1}"

ame --o $OUTPUT_DIR/Liver_vs_Heart_top --control $FASTA_DIR/ded_top500_cov2_Heart.fa $FASTA_DIR/ded_top500_cov2_Liver.fa $MOTIF_FILE