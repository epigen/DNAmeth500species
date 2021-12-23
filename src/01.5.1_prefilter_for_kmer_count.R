#source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/02.0_predict_initialization.R"))
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

##which species to run on 
args = commandArgs(trailingOnly = TRUE)
species = args[1]

##create output
subdir = file.path(analysis_dir, "01_basicStats/01.99_filtered_sequences", species) ###not saving it to a specific folder, as might be needed for other analysis
dir.create(subdir,recursive = TRUE)

wd = file.path(processed_dir,species)
setwd(wd)

uc_map = fread(list.files(path="toSelf_filtered_0.08mm_final_concat",
                        pattern=".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$", recursive=TRUE, full.names=TRUE))

ded_ref = readDNAStringSet("reduced/consensus//toSelf_filtered_0.08mm_final.fa")
ded_ref_cgg = as.character(substring(ded_ref,1,3))=="CGG"
cgg_names = names(ded_ref[ded_ref_cgg])
cgg_names_match_uc <- cgg_names[cgg_names %in% uc_map[V13==1]$V4]

writeXStringSet(ded_ref[cgg_names_match_uc], file.path(subdir, "toSelf_filtered_0.08mm_final_filtered.fa"))

df<-data.frame(length(ded_ref), length(cgg_names), length(cgg_names_match_uc))
colnames(df)<-c("initial", "only_cgg", "cgg_and_match_uc")

write.csv(df, file.path(subdir, "filter_stats.csv"), quote = F, row.names = F)
