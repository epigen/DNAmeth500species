## general overview and heatmaps for the inverted species detected

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
wd <- file.path(analysis_dir, "02_predict_meth/02.8_why_inv/summary")
dir.create(wd, recursive = T)
setwd(wd)
##uploading data
sp_full_order <- as.character(read.table(paste0(Sys.getenv("CODEBASE"), "compEpi/meta/species_ordered.txt"), header = 1)$sp_full_order)
sp_df <- as.data.frame(sp_full_order)

colnames(sp_df) <- c("species")
sp_df <- unique(left_join(sp_df, stats_annot[, c("species", "color_class")]))
row.names(sp_df) <- sp_df$species


#write.table(sp_df[sp_df$color_class == "Actinopteri","species"],
#          file.path(Sys.getenv("CODEBASE"),"compEpi/meta/Actinopteri.txt"), quote = F, row.names = F, col.names = F)
#inverted file 
invfile <- paste0(Sys.getenv("CODEBASE"), "compEpi/meta/inverted_species.txt")
auc_file <- paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/all_aucs_full.csv")
all_aucs<-read.csv(auc_file, row.names = 1)


full_auc_m <- read.csv(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/all_aucs_full_matrix.csv"), row.names = 1)
full_auc_m <- full_auc_m[, complete.cases(t(full_auc_m))]
  
##identifying the potential inverted species
full_auc_m <- as.matrix(full_auc_m)
  
s <- colMeans(full_auc_m)
s_rows <- rowMeans(full_auc_m)
inv <- union(names(s_rows[s_rows<0.5]), names(s[s<0.5]))
sd <- apply(full_auc_m, 1,sd)[inv]
  
#we want only those, that are not always aroud 0.5
s_rows[inv]-sd
  
## Then I manually check the ROC curves of the species in the folder
##FLA, NOP is mostly around 0.5 - not inverted, just not working
inv <- inv[!inv %in% c("FLA", "NOP", "YFS", "RBF") ]
  ##SBS should be treated carefully
if(!file.exists(invfile)){
  write.table(as.data.frame(inv), invfile, quote = F, row.names = F, col.names = F)
}

inv <- as.character(read.table(invfile)$V1)

sp_df$inv <- sp_df$species %in% inv


##number of inverted specis as a % of each phylo group
pdf("inverted_count.pdf", width = 6, height = 5)
ggplot(sp_df, aes(x = color_class, fill = inv)) + geom_bar(stat = "count", position = "stack") + 
  rotate_labels() + scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) + 
  labs(x = "phylogenetic group")
dev.off()

## heatmap inverted vs inverted
inv_ord <- sp_full_order[sp_full_order %in% inv]
auc_inv <- full_auc_m[inv_ord, inv_ord]

seq_count <- unique(read.csv(paste0(analysis_dir, "/02_predict_meth/02.1_within_species/summary/all_aucs.csv"), row.names = 1)[, c("species", "numSequences", "color_class")])
row.names(seq_count) <- seq_count$species
colnames(seq_count) <- c("species", "seqN", "class")
col_fun_num = colorRamp2(c(0, 2000), c("white", "darkgrey"))
col_annot<-columnAnnotation(class = sp_df[inv_ord, "color_class"], seqN = seq_count[inv_ord, "numSequences"], 
                            col = list(class = class_colors, seqN = col_fun_num), show_legend = c(F,T))
row_annot<-rowAnnotation(class = sp_df[inv_ord, "color_class"], col = list(class = class_colors),  show_legend = F)
col_fun = colorRamp2(c(0, 0.4,0.5,0.6, 1), c("#2b83ba","#abd9e9", "#ffffbf", "#fdae61","#d7191c"))

pdf("inv_vs_inv.pdf", width = 6, height = 5)
row_annot + Heatmap(auc_inv, top_annotation = col_annot, name = "ROC - AUC", cluster_rows = F, 
                    cluster_columns = F, col = col_fun) 
dev.off()



makeColorRampPalette <- function(colors, points, mat, num.colors.in.palette)
{
  stopifnot(length(colors) == 3)
  stopifnot(length(points) == 3)
  cutoff.fraction <- (points[[2]] - min(mat))/(max(mat) - min(mat))
  print(cutoff.fraction)
  bottom_color <- colorRampPalette(colors[1:2])(200)[round((min(mat) - points[[1]])/(points[[2]] - points[[1]])*200)]
  print(round((min(mat) - points[[1]])/(points[[2]] - points[[1]])*200))
  ramp1 <- colorRampPalette(c(bottom_color, colors[2]))(num.colors.in.palette * cutoff.fraction)
  
  top_color <- colorRampPalette(colors[2:3])(200)[round((max(mat) - points[[2]])/(points[[3]] - points[[2]])*200)]
  ramp2 <- colorRampPalette(c(colors[2], top_color))(num.colors.in.palette * (1 - cutoff.fraction))
  
  return(c(ramp1, ramp2))
}


cols <- makeColorRampPalette(c( "#abd9e9",    # distances 0 to 3 colored from white to red
                               "#ffffbf", "#d7191c"), # distances 3 to max(distmat) colored from green to black,
                             c(0.4, 0.5, 1), # vector of required color points
                             auc_inv, #matrix to draw the heatmap on
                             200)

my_color = list(class = class_colors, seqN = colorRampPalette(c("white", "darkgrey"))(200))

pdf("inv_vs_inv_upgraded.pdf", width = 8, height = 8)
pheatmap(auc_inv, cluster_rows = F, cluster_cols = F, color = cols, 
         annotation_colors = my_color,
         annotation_col = seq_count[, c("seqN", "class")],
         annotation_row = seq_count[, c("seqN", "class")],
         gaps_row = c(1,5, 10), 
         gaps_col = c(1,5,10))
dev.off()

###summary - for each species draw a distribution of AUCs in the species vs the random
pdf("inv_normal_aucs_distribution.pdf", width = 50, height = 100)
ggplot(all_aucs[all_aucs$ifRand=="noRandTest" | all_aucs$ifRand=="RandTest", ], aes(x = auc, fill = ifRand)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c("noRandTest" = "lightblue", "RandTest" = "grey", "noRandTrain" = "red", "RandTrain" = "black")) +
  geom_point(data = all_aucs[all_aucs$ifRand=="RandTrain" | all_aucs$ifRand=="noRandTrain", ],
             aes(y = 0, x = auc, fill = ifRand), shape = 21) +
  facet_wrap(~train_species, ncol=20)
dev.off()

sub_all_aucs <- all_aucs[all_aucs$train_species %in% c("SWF", "WHH"), ]
pdf("inv_normal_aucs_distribution_sub.pdf", width = 10, height = 5)
ggplot(sub_all_aucs[sub_all_aucs$ifRand=="noRandTest" | sub_all_aucs$ifRand=="RandTest", ], aes(x = auc, fill = ifRand)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c("noRandTest" = "lightblue", "RandTest" = "grey", "noRandTrain" = "red", "RandTrain" = "black")) +
  geom_point(data = sub_all_aucs[sub_all_aucs$ifRand=="RandTrain" | sub_all_aucs$ifRand=="noRandTrain", ],
             aes(y = 0, x = auc, fill = ifRand), shape = 21) +
  facet_wrap(~train_species, ncol=2) + labs(x = "ROC - AUC")
dev.off()


classes_with_inv <- unique(stats_annot[stats_annot$species %in% inv, ]$color_class)
species_for_analysis <- unique(stats_annot[stats_annot$color_class %in% classes_with_inv, ]$species)
if(!file.exists(paste0(Sys.getenv("CODEBASE"), "/compEpi/meta/classes_with_inv.txt"))){
  write.table(as.data.frame(species_for_analysis), paste0(Sys.getenv("CODEBASE"),
                                                          "/compEpi/meta/classes_with_inv.txt"), row.names = F, col.names = F, quote = F)
}
