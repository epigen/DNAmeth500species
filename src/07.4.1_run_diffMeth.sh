##specifying the species:
DATA_DIR="/binfl/lv71484/droman/DNAmeth500species/"
wd=${DATA_DIR}"/results_pipeline/${1}"
RRBSdir="toSelf_filtered_0.08mm_final_concat"
species=${1}
dedRef_id="toSelf_filtered_0.08mm_final_concat"
sampleAnnotation=$DATA_DIR"/results_pipeline/${species}/meta/sampleAnnotation.tsv"
comp_col="new_comp"
groups_col="Tissue"
nTopDiffMeth=500
scripts="${CODEBASE}/DNAmeth500species/src"
motif="cpg"
uc_tag="_uc"
output_dir=$DATA_DIR"/results_analysis/07_motifAnalysis/07.4_diffMeth_additional_run/${1}"
echo $output_dir

/home/lv71484/droman/.conda/envs/Zoo_R/bin/Rscript ${CODEBASE}/DNAmeth500species/src/diffMeth.R ${wd} ${RRBSdir} $species $dedRef_id $sampleAnnotation $comp_col $groups_col $nTopDiffMeth $scripts $motif $uc_tag $output_dir

