source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"validation", "01_crossMapping", "summary")
dir.create(wd,recursive=TRUE)
setwd(wd)

##collecting all the files
file_list <- system(paste0("ls ", processed_dir, "/*/crossMapping/*_concatinated/toSelf_filtered_0.08mm.bsmap.0.2.mism.1.r.0.n.bam.flagstat"), intern = TRUE)

cross_mapping_stats <- data.table()

#parsing them one by one
for(filepath in file_list){
test_file <-fread(filepath, fill = TRUE)

perc_mapped = as.character(test_file[5,5])
perc_mapped <- gsub(":", "", gsub("%", "", gsub("[()]","",perc_mapped)))

temp_df <- data.table("total_passed" = test_file[1,1], "total_failed" = test_file[1,3],
           "secondary_passed" = test_file[2,1], "secondary_failed" = test_file[2,3],
          "mapped_passed" = test_file[5,1], "mapped_failed" = test_file[5,3],
          "perc_mapped_passed" = as.numeric(strsplit(perc_mapped, "-")[[1]][1]),
          "perc_mapped_failed" = as.numeric(strsplit(perc_mapped, "-")[[1]][2]),
          "species" = strsplit(gsub(processed_dir, "", filepath), "/")[[1]][2],
          "mapped_genome" = gsub("_concatinated", "", strsplit(gsub(processed_dir, "", filepath), "/")[[1]][4]))

cross_mapping_stats <- rbind(cross_mapping_stats, temp_df)
}

head(cross_mapping_stats)

cross_mapping_stats <- left_join(cross_mapping_stats, sp_df)

my_wt(cross_mapping_stats, "cross_mapping_stats.tsv")

ggplot(cross_mapping_stats[!is.na(color_class)], aes(x = perc_mapped_passed, fill = color_class)) + geom_histogram() + facet_wrap(~color_class, ncol = 4) + 
scale_fill_manual(values = class_colors)
ggsave("passed_QC_per_class.pdf", width = 12, height = 6 )


ggplot(cross_mapping_stats, aes(x = perc_mapped_passed)) + geom_histogram()
ggsave("passed_QC.pdf", width =4, height = 4 )

cross_mapping_stats_best <- cross_mapping_stats %>% arrange(desc(perc_mapped_passed)) %>%
group_by(species) %>%
filter(perc_mapped_passed == first(perc_mapped_passed))
setDT(cross_mapping_stats_best)

ggplot(cross_mapping_stats_best[!is.na(color_class)], aes(x = perc_mapped_passed,fill = color_class)) + geom_histogram() + facet_wrap(~color_class, ncol = 4) + 
scale_fill_manual(values = class_colors)
ggsave("passed_QC_per_class_best.pdf", width = 12, height = 6 )

