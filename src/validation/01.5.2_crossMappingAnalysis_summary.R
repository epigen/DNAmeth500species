## ---------------------------
##
## cross-mapping analysis
##
## Summary of the per genome-species pair analysis
##
## Authors: Daria Romanovskaia
##
## Date Created: 2021-10-28
##
##
## ---------------------------

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd = file.path(analysis_dir, "validation", "01_crossMapping","01.5_analysis")
setwd(wd)

dir.create("summary")

cm_stats <- fread("../01.4_crossMapping_stats/cross_mapping_stats.tsv")
head(cm_stats)

cross_mapping_stats_best <- cm_stats %>% arrange(desc(perc_mapped_passed)) %>%
group_by(species) %>%
filter(perc_mapped_passed == first(perc_mapped_passed))
setDT(cross_mapping_stats_best)

head(cross_mapping_stats_best)

cross_mapping_stats_best[species == "AX"]

cross_mapping_stats_best[, folder:=paste0(species, "_", mapped_genome), by = row.names(cross_mapping_stats_best)]

summary_stats = data.table()
for(i in 1:NROW(cross_mapping_stats_best)){
    
    df <- try(fread(file.path(cross_mapping_stats_best$species[[i]], cross_mapping_stats_best$folder[[i]], 
                          paste0(cross_mapping_stats_best$species[[i]], "_methylation_profile_stats.tsv"))))
    if(is.list(df)){
    df$mapped_genome <- cross_mapping_stats_best$mapped_genome[[i]]
    df$color_class <- cross_mapping_stats_best$color_class[[i]]
        
    summary_stats <- rbind(summary_stats, df)
        
    }else{
        print(cross_mapping_stats_best$folder[[i]])
    }
    
}

head(summary_stats)

summary_stats <- summary_stats[label %in% stats_annot$Sample_Name,]

summary_stats[, species:=strsplit(label, "_")[[1]][1], by = row.names(summary_stats)]

head(cross_mapping_stats_best)

summary_stats <- left_join(summary_stats, cross_mapping_stats_best[, c("mapped_genome", "species", "perc_mapped_passed")])

my_wt(summary_stats, "summary/Methylation_profile_stats_summary.tsv")

getwd()

"/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/01_crossMapping/01.5_analysis/CM/CM_calJac4/"

full_stats = data.table()
for(i in 1:NROW(cross_mapping_stats_best)){
    
    df <- try(fread(file.path(cross_mapping_stats_best$species[[i]], cross_mapping_stats_best$folder[[i]], 
                          paste0(cross_mapping_stats_best$species[[i]], "_methylation_profile.tsv"))))
    if(is.list(df)){
    df$mapped_genome <- cross_mapping_stats_best$mapped_genome[[i]]
    df$color_class <- cross_mapping_stats_best$color_class[[i]]
        
    full_stats <- rbind(full_stats, df)
        
    }else{
        print(cross_mapping_stats_best$folder[[i]])
    }
    
}

object.size(full_stats)

NROW(full_stats)

my_wt(full_stats, "/binfl/lv71484/droman/DNAmeth500species/results_analysis/validation/01_crossMapping/01.5_analysis/summary/Methylation_profiles_merged.tsv")


