#!/bin/bash

## $1 - the species to run the analysis on
DATA_DIR=<your path here>
MOTIF_FILE=<your path here>"/JASPAR/JASPAR_CORE_2020_vertebrates.meme"
FASTA_DIR=$DATA_DIR"/results_analysis/03_motifAnalysis/diffMeth_additional_run/${1}/motifAnalysis/"
OUTPUT_DIR=$DATA_DIR"/results_analysis/03_motifAnalysis/03.7_AME/screen/${1}"
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR"/Heart_vs_Liver_top" 
mkdir "$OUTPUT_DIR/Liver_vs_Heart_top" 

echo "running Heart vs Liver ${1}"

ame --o $OUTPUT_DIR/Heart_vs_Liver_top --control $FASTA_DIR/ded_top500_cov2_Liver.fa $FASTA_DIR/ded_top500_cov2_Heart.fa $MOTIF_FILE

echo "running Liver vs Heart ${1}"

ame --o $OUTPUT_DIR/Liver_vs_Heart_top --control $FASTA_DIR/ded_top500_cov2_Heart.fa $FASTA_DIR/ded_top500_cov2_Liver.fa $MOTIF_FILE
