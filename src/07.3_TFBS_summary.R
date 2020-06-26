source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(TFBSTools)

wd=file.path(analysis_dir, "03_motifAnalysis/")
setwd(wd)
dir.create("summary")

jaspar_path = <your path to jaspar database>
motif_list <- names(readJASPARMatrix(fn =jaspar_path, 
          type = "all"))


###reading in the stats files (results of 03.2) - should be performed only once:
if(!file.exists("03_summary/TF_count_per_species_normalized_2020.csv")){
  stats_files = system("ls /scratch/lab_bock/shared/projects/compEpi/results_analysis/03_motifAnalysis/03.2_TF_freq_2020/*/counts.csv",intern=TRUE)
  ##list of species from filenames:
  sp_list <- unlist(lapply(stats_files, function(x) strsplit(x, "/")[[1]][10] ))
  
  TF_norm <- lapply(seq_along(stats_files), 
                  function(i) {df<-read.csv(stats_files[[i]])[,c("TF",  "n_norm")];  
                  colnames(df)<-c("TF", sp_list[[i]] ); 
                  return(df)})
  
  TF_norm_df <- Reduce(full_join, TF_norm)
  write.csv(TF_norm_df, "/03_motifAnalysis/03_summary/TF_count_per_species_normalized_2020.csv", quote = F)
}else{
  TF_norm_df <- read.csv( "/03_motifAnalysis/03_summary/TF_count_per_species_normalized.csv", row.names=1)
}

##assigning the TF ids as row.names 
row.names(TF_norm_df) <- TF_norm_df$TF
TF_norm_df <- TF_norm_df[motif_list, -1]

##NA indicates, that none of the TFBS where detected (i.e. no file created) -> replacing with 0s
TF_norm_df[is.na(TF_norm_df)] <- 0

TF_norm_df_upd <- TF_norm_df
TF_norm_df_upd[TF_norm_df_upd>1] <- 1


## TF count vs species
col_annot<-data.frame(sapply(sp_list, 
                             function(x) stats_annot[stats_annot$species==x, ]$color_class[[1]]))
colnames(col_annot)<-c("phyl.group")


pdf("summary/TF_frequency.pdf", width = 20, height = 10)
pheatmap(t(TF_norm_df_upd)[sp_df[sp_df$species %in% sp_list,]$species,],
         cluster_rows = F, show_rownames = F, fontsize_col = 5,
         annotation_row = col_annot, annotation_colors = list(phyl.group=class_colors),
         annotation_legend = F,  color = colorRampPalette(c("white", "blue"))(50)) 
dev.off()



##scaling each TF across species, so that they would be comparable
TF_norm_df_scaled<-scale(t(TF_norm_df), center = T, scale = T)

pdf("summary/TF_frequency_z_score.pdf", width = 20, height = 10)
pheatmap(TF_norm_df_scaled[sp_df[sp_df$species %in% sp_list,]$species,],
         cluster_rows = F, show_rownames = F, fontsize_col = 5,
         annotation_row = col_annot, annotation_colors = list(phyl.group=class_colors),
         annotation_legend = F,
         name = "freq. z-score") 
dev.off()


## how many species have heart and liver and more than one replica? 
replicas <- stats_annot[, c("species", "Tissue", "Sample_Name", "color_class", "conversion_type")] %>%
  filter(conversion_type == "converted") %>%
  distinct() %>%
  filter(Tissue == "Heart" | Tissue == "Liver") %>%
  dplyr::group_by(color_class, species, Tissue) %>%
  dplyr::summarise(n = n()) %>%
  spread(Tissue, n)

replicas <- as.data.frame(replicas)
replicas <- replicas[complete.cases(replicas),]
replicas$success <- replicas$species %in% sp_list


if(!file.exists("number_of_replicas.pdf")){
  pdf("number_of_replicas.pdf", height = 2, width = 10)
  ggplot(replicas, aes(x = Heart, y = Liver, color = success)) + geom_jitter() +
    geom_text_repel(data = replicas[!replicas$success & replicas$Heart>1 &replicas$Liver>1, ], 
                    aes(x = Heart, y = Liver, label = species)) + facet_wrap(~color_class, nrow = 1)
  dev.off()
}

## for further analysis we want only species that have at least two replicas in both heart & liver
head(replicas)
sp_two_repl <- replicas[replicas$Heart>1 & replicas$Liver>1,]$species
write.table(as.data.frame(sp_two_repl), 
            file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/meta/species_list_TF_analysis.txt"), 
            row.names = F, quote = F, col.names = F)

