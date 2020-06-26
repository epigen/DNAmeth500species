##specifying the species:
DATA_DIR=<your path here>
wd=${DATA_DIR}"/results_pipeline/${1}"
RRBSdir="toSelf_filtered_0.08mm_final_concat"
species=${1}
dedRef_id="toSelf_filtered_0.08mm_final_concat"
sampleAnnotation=$DATA_DIR"/results_pipeline/${species}/meta/sampleAnnotation.tsv"
comp_col="new_comp"
groups_col="Tissue"
nTopDiffMeth=500
scripts="${CODEBASE}/compEpi/src"
motif="cpg"
uc_tag="_uc"
output_dir=$DATA_DIR"/results_analysis/03_motifAnalysis/diffMeth_additional_run/${1}"
echo $output_dir

${CODEBASE}/compEpi/src/diffMeth.R ${wd} ${RRBSdir} $species $dedRef_id $sampleAnnotation $comp_col $groups_col $nTopDiffMeth $scripts $motif $uc_tag $output_dir

