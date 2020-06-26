source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(RColorBrewer)
wd = file.path(analysis_dir, "03_motifAnalysis/03.7_AME/SELEX")
dir.create(wd, recursive = T)
setwd(file.path(analysis_dir, "03_motifAnalysis/03.7_AME/"))

ans <- read.csv("summary/tf_presence.csv", sep = " ", row.names = 1)
head(ans)



##reading the selex data (Yin et al, S. table 3)
selex <- read.csv("../resources/Selex_data.csv", skip = 20, header = 21, sep = ";")
head(selex)

length(intersect(row.names(ans), selex$TF.name))

TF_meth <- data.frame(row.names(ans))
colnames(TF_meth) <- c("TF")

TF_meth$TF_upd <- unlist(sapply(as.character(TF_meth$TF), function(x) strsplit(x,"\\(")[[1]][1]))

##transferring the mouse gene names to human
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
musGenes <- TF_meth$TF_upd[sapply(TF_meth$TF_upd, function(x) x!=toupper(x))]
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", musGenes, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- setNames(genesV2[, 2], genesV2[, 1])

TF_meth$TF_upd <- unlist(sapply(TF_meth$TF_upd, function(x){if(x %in% names(humanx)) return(humanx[[x]]) 
                                                                        else return(x)}))
length(intersect(TF_meth$TF_upd, selex$TF.name))

TF_meth <- unique(left_join(TF_meth, selex[, c("TF.name","clone.type", "Call")], by = c("TF_upd" = "TF.name")))

ggplot(TF_meth, aes(x = Call)) + geom_bar() + rotate_labels()
ggsave("Meth_preferences.pdf", width = 5, height = 6)

##now lets assume that all the little effects and inconclusive can be united with undefined into grey area
TF_meth$Call_upd <- sapply(as.character(TF_meth$Call), 
                           function(x){if(!is.na(x) & (x=="MethylPlus"|x=="MethylMinus")) return(x)
                                                                      else return("undefined/unknown")})
                           
                           
##same heatmaps as before, but with methylation preferences
##species annot
annot_col = as.data.frame(sp_df$color_class)
row.names(annot_col) <- row.names(sp_df)
colnames(annot_col) <- c("class")

color=c(tissue_colors["L"], "#f7f7f7", tissue_colors["H"])
breaks=c(-2, -1,0,1)

##add TF family annotation
tf_annot <- read.csv("../JASPAR.csv", sep = ";")
all_tf <- row.names(ans)
head(tf_annot)
my_tf_annot <- unique(tf_annot[tf_annot$Name %in% all_tf, c("Name", "Class", "Family")])
row.names(my_tf_annot) <- my_tf_annot$Name
my_tf_annot$Family <- as.character(my_tf_annot$Family)

new_f <- sapply(my_tf_annot$Family, function(x) strsplit(x, "::")[[1]])
my_tf_annot$Family_short<-as.character(sapply(new_f, function(x) x[length(x)]))

new_c <- sapply(as.character(my_tf_annot$Class), function(x) strsplit(x, "::")[[1]])
my_tf_annot$Class_short<-as.character(sapply(new_c, function(x) x[length(x)]))


my_tf_annot <- my_tf_annot[order(my_tf_annot$Class_short),]
my_tf_annot$Family_short <- factor(my_tf_annot$Family_short, levels = unique(my_tf_annot$Family_short))

my_tf_annot_upd <- left_join(my_tf_annot, unique(TF_meth[, c("TF", "TF_upd", "Call_upd")]), by = c("Name" = "TF"))
#manually curating and dropping duplicates
##BCL6
id_to_drop <- c(row.names(my_tf_annot_upd[my_tf_annot_upd$Name=="BCL6" & my_tf_annot_upd$Call_upd=="undefined/unknown",]))

id_to_drop <- c(id_to_drop, row.names(my_tf_annot_upd[my_tf_annot_upd$Name=="KLF17" & 
                                          my_tf_annot_upd$Call_upd=="undefined/unknown",]))
id_to_drop <- c(id_to_drop, row.names(my_tf_annot_upd[my_tf_annot_upd$Name=="PAX7" & 
                                                        my_tf_annot_upd$Call_upd!="undefined/unknown",]))
id_to_drop

my_tf_annot_upd <- my_tf_annot_upd[-as.integer(id_to_drop),]
row.names(my_tf_annot_upd) <- my_tf_annot_upd$Name
#creating a distinct color palette
n <- length(unique(my_tf_annot$Class_short))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Class_idst_col<-setNames(sample(col_vector, n), unique(my_tf_annot$Class_short))
write.csv(as.data.frame(Class_idst_col), "classes_color_annot.csv", quote = F)
write.table(my_tf_annot_upd, "SELEX/annotation_TFs.csv",sep = ",", row.names = F )
my_colors <- list(class = class_colors, Class_short = Class_idst_col, 
                  Call_upd = c("MethylPlus" = "red", "MethylMinus" = "#006837", "undefined/unknown" = "grey"))
dev.off()
pdf("summary/final_heatmaps/TFBS_presence_turned_selex.pdf", width = 35, height = 12)
pheatmap(t(ans[,sp_df$species[sp_df$species %in% colnames(ans)]]), 
         cluster_rows = F,annotation_row = annot_col, annotation_col = my_tf_annot_upd[, c(5,7)],
         annotation_colors = my_colors, show_rownames = F,
         fontsize = 5, annotation_legend = T , breaks = breaks, color = color)
dev.off()

##add filtering
###and lets filter so that TF shold be detected at 10 species (5% of the species) and speceis, 
###that have at least 10 signals (10%)
tf_to_keep <- names(rowSums(abs(ans))[rowSums(abs(ans)) >= 10])
ans_f <- ans[tf_to_keep,]
sp_to_keep <- names(colSums(abs(ans_f))[colSums(abs(ans_f)) >= 5])
ans_f <- ans_f[,sp_to_keep]
f <- unique(my_tf_annot[my_tf_annot$Name %in% row.names(ans_f),]$Family_short)
my_colors <- list(class = class_colors, 
                  Call_upd = c("MethylPlus" = "red", "MethylMinus" = "#006837", "undefined/unknown" = "grey"),
                  Class_short = Class_idst_col[unique(my_tf_annot[my_tf_annot$Name %in% row.names(ans_f),]$Class_short)])

dev.off()
pdf("summary/final_heatmaps/TFBS_presence_turned_filtered_selex.pdf", width = 15, height = 10)
pheatmap(t(ans_f), 
         annotation_row = annot_col, annotation_col = my_tf_annot_upd[,  c(5,7)],
         annotation_colors = my_colors, show_rownames = F,
         fontsize = 10, breaks = breaks, color = color)
dev.off()

ph <- pheatmap(t(ans_f))

tf_tree <- ph$tree_col
tf_tree_dend <- as.dendrogram(tf_tree)

library(dendextend)

pdf("summary/TF_clusters_2.pdf",  height = 5, width = 12)
tf_tree_dend %>% set("labels_col", value = c("green", "blue", "red", "black", "orange"), k=2) %>% 
  plot(main = "Color labels per cluster in TFs")

dev.off()


cl <- cutree(tf_tree, k = 2)
df <- as.data.frame(cl)
df$TF <- row.names(df)
df <- left_join(df, TF_meth[, c("TF", "TF_upd", "Call_upd")])

write.table(df, "summary/clusters_and_selex.csv", row.names = F, quote = F)


###saving presence of FOXO4  and EGR1 nicely
ans_sub <- as.data.frame(t(ans[c("FOXO4", "EGR1"),]))
ans_sub_sub <- ans_sub[ans_sub$FOXO4!=0 | ans_sub$EGR1!=0, ]
annot <- sp_df[sp_df$species %in% row.names(ans_sub_sub), "color_class", drop = F]
col = list(color_class = class_colors)
pdf("summary/final_heatmaps/TFBS_selected.pdf", width = 10, height = 10)
pheatmap(ans_sub_sub[row.names(annot),], 
         cluster_rows = F,
          cluster_cols = F,
          annotation_row = annot,
         annotation_colors = col,
         fontsize = 10, annotation_legend = T , breaks = breaks, color = color)
dev.off()

ans_sub_sub <- cbind(ans_sub_sub, annot)
ans_sub_sub$species <- row.names(ans_sub_sub)
a <- melt(ans_sub_sub, id.vars = c("species", "color_class"))

a_summ <- a %>%
  group_by(variable, color_class, value) %>%
  filter(color_class %in% c("Aves","Marsupialia", "Mammalia" )) %>% 
  summarize(n = n())
a_summ[a_summ$color_class == "Marsupialia", ]$color_class <- "Mammalia"
a_summ <- as.data.frame(a_summ)
a_summ$Tissue <- "None"

a_summ[a_summ$value == 1,]$Tissue <- "Heart"
a_summ[a_summ$value == -1,]$Tissue <- "Liver"
a_summ$Tissue <- factor(a_summ$Tissue, levels = c("Liver", "None", "Heart"))
ggplot(a_summ, aes(x = variable, y = n, fill = Tissue)) + 
  geom_bar(stat = "identity", position="stack") + 
  scale_fill_manual(values = c("Liver" = "#e6ab02" ,"None" = "lightgrey", "Heart" = "#005fb1" )) + 
  theme(text = element_text(size = 15)) + labs(x="") +
  facet_wrap(~color_class, ncol =1, scales = "free") #+
ggsave("summary/final_heatmaps/TFBS_selected_count.pdf", width = 5, height = 10)
