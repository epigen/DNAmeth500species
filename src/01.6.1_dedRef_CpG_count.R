#!/usr/bin/env Rscript
##classical initiation
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

#redirecting to the stats_detail (where it is going to be analyzed)
wd=file.path(analysis_dir,"01_basicStats/01.6_stats_detail")
setwd(wd)

#function
get_cpg_count <- function(SP){
  print(SP)
  ##uploading fasta file
  ded_ref=readDNAStringSet(file.path(processed_dir, SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  ##initial count of dedRef fragments
  N <- length(ded_ref)
  #counting the CG nuclotides 
  CpG_count <-  sum(dinucleotideFrequency(ded_ref)[,"CG"])
  #filtering based on the restriction enzyme sequence
  ded_ref_cgg=as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names=names(ded_ref[ded_ref_cgg])
  ded_ref = ded_ref[ded_ref_cgg]
  #counting the CG nuclotides after filtering
  CpG_count_filtered <-  sum(dinucleotideFrequency(ded_ref)[,"CG"])
  ##count of the filtered dedRef fragments
  N_new <- length(ded_ref)
  return(c("dedRef_count_full" = N, "dedRef_count_filtered" = N_new, 
           "CpG_count" = CpG_count, "CpG_count_filtered" = CpG_count_filtered ))
}
#lazy implementation - not the fastest, but no need to create additional files and merging them after
ans <- sapply(sp_df$species, get_cpg_count)
#saving 
write.table(t(ans), "dedRef_CpG_count.csv", quote = F, sep = ";")
