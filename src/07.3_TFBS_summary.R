source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(TFBSTools)

wd=file.path(analysis_dir, "07_motifAnalysis/")
setwd(wd)
dir.create("07.3_freq_summary")

jaspar_path = file.path(data_dir, "resources", "JASPAR", "JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
motif_list <- names(readJASPARMatrix(fn =jaspar_path))
jaspar_annot <- read.csv(file.path(data_dir, "resources", "JASPAR","JASPAR.csv"), sep = ';')

###reading in the stats files (results of 03.2) - should be performed only once:
if(!file.exists("07.3_freq_summary/TF_count_per_species_normalized_2020.csv")){
    stats_files = list.files("07.2_TF_freq_2020", pattern = "counts.csv", recursive = T)
  ##list of species from filenames:
  sp_list <- unlist(lapply(stats_files, function(x) strsplit(x, "/")[[1]][1] ))
  
  TF_norm <- lapply(seq_along(stats_files), 
                  function(i) {df<-read.csv(file.path("07.2_TF_freq_2020",stats_files[[i]]))[,c("TF",  "n_norm")];  
                  colnames(df)<-c("TF", sp_list[[i]] ); 
                  return(df)})
  
  TF_norm_df <- Reduce(full_join, TF_norm)
  my_wt(TF_norm_df, "07.3_freq_summary/TF_count_per_species_normalized_2020.csv")
}else{
  TF_norm_df <- fread( "07.3_freq_summary/TF_count_per_species_normalized.csv")
}

                                 
##assigning the TF ids as row.names 
row.names(TF_norm_df) <- TF_norm_df$TF
TF_norm_df <- TF_norm_df[motif_list, -1]
                           
human_tfs <- jaspar_annot[sapply(jaspar_annot$Species, function(x) "Homo sapiens" %in% x),]$Name
human <- intersect(human_tfs, row.names(TF_norm_df))
TF_norm_df <- TF_norm_df[human, ]


##NA indicates, that none of the TFBS where detected (i.e. no file created) -> replacing with 0s
TF_norm_df[is.na(TF_norm_df)] <- 0

TF_norm_df_upd <- TF_norm_df
TF_norm_df_upd[TF_norm_df_upd>1] <- 1


## TF count vs species
col_annot<-data.frame(sapply(sp_list, 
                             function(x) stats_annot[stats_annot$species==x, ]$color_class[[1]]))
colnames(col_annot)<-c("phyl.group")


pdf("07.3_freq_summary/TF_frequency.pdf", width = 20, height = 10)
pheatmap(t(TF_norm_df_upd)[sp_df[sp_df$species %in% sp_list,]$species,],
         cluster_rows = F, show_rownames = F, fontsize_col = 5,
         annotation_row = col_annot, annotation_colors = list(phyl.group=class_colors),
         annotation_legend = F,  color = colorRampPalette(c("white", "blue"))(50)) 
dev.off()

##stats summary
TF_norm_df_melted <- left_join(TF_norm_df_melted, sp_df)
TF_norm_df_melted[TF_norm_df_melted$species == "BLC", ]$color_class <- "Actinopteri"
TF_norm_df_melted_summary <- TF_norm_df_melted %>% 
    group_by(color_class) %>%
    summarize(mean_fr = mean(freq), sd = sd(freq))
my_wt(TF_norm_df_melted_summary, "07.3_freq_summary/TF_frequency_stats_per_class_human.tsv")                             

##scaling each TF across species, so that they would be comparable
TF_norm_df_scaled<-scale(t(TF_norm_df), center = T, scale = T)

pdf("07.3_freq_summary/TF_frequency_z_score.pdf", width = 20, height = 10)
pheatmap(TF_norm_df_scaled[sp_df[sp_df$species %in% sp_list,]$species,],
         cluster_rows = F, show_rownames = F, fontsize_col = 5,
         annotation_row = col_annot, annotation_colors = list(phyl.group=class_colors),
         annotation_legend = F,
         name = "freq. z-score") 
dev.off()
                             
range <- max(abs(TF_norm_df_scaled))
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
#myBreaks <- c(seq(min(TF_norm_df_scaled), 0, length.out=ceiling(paletteLength/2) + 1), 
 #             seq(max(TF_norm_df_scaled)/paletteLength, max(TF_norm_df_scaled), length.out=floor(paletteLength/2)))

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(5/paletteLength, 5, length.out=floor(paletteLength/2)))


dev.off()
pdf("07.3_freq_summary/TF_frequency_z_score_withthr_white_human.pdf", width = 15, height = 10)
pheatmap(TF_norm_df_scaled[sp_df[sp_df$species %in% colnames(TF_norm_df),]$species,],
         cluster_rows = F, show_rownames = F, fontsize_col = 3,
         annotation_row = col_annot, annotation_colors = list(phyl.group=class_colors),
         annotation_legend = F,
         breaks = myBreaks,
         color = myColor,
         name = "freq. z-score") 
dev.off()
                

                             

TF_norm_df_scaled[TF_norm_df_scaled>5] <- 5
TF_norm_df_scaled[TF_norm_df_scaled<-5] <- -5
                             
try(dev.off())
pdf("07.3_freq_summary/TF_frequency_z_score_withthr.pdf", width = 20, height = 10)
pheatmap(TF_norm_df_scaled[sp_df[sp_df$species %in% colnames(TF_norm_df),]$species,],
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


#### only mouse TFs:
mouse_tfs <- jaspar_annot[sapply(jaspar_annot$Species, function(x) "Mus musculus" %in% x),]$Name
mouse <- intersect(mouse_tfs, row.names(TF_norm_df))
TF_norm_df_scaled_mouse <- scale(t(TF_norm_df[mouse, ]), center = T, scale = T)

                                 
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(5/paletteLength, 5, length.out=floor(paletteLength/2)))


dev.off()
pdf("07.3_freq_summary/TF_frequency_z_score_withthr_white_mouse.pdf", width = 15, height = 10)
pheatmap(TF_norm_df_scaled_mouse[sp_df[sp_df$species %in% row.names(TF_norm_df_scaled_mouse),]$species,],
         cluster_rows = F, show_rownames = F, fontsize_col = 10,
         annotation_row = col_annot, annotation_colors = list(phyl.group=class_colors),
         annotation_legend = F,
         breaks = myBreaks,
         color = myColor,
         name = "freq. z-score") 
dev.off()