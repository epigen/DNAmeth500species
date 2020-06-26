##!!! here everywhere we use the tissues the TFs where detected at, but the actual biological meaning is the opposite!
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(plyr)
library(RColorBrewer)


wd = file.path(analysis_dir, "03_motifAnalysis/03.7_AME")
setwd(wd)


## used to generate heatmaps
makeColorRampPalette_2 <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  
  ramp1<-rep(colors[[1]], num.colors.in.palette*cutoff.fraction)
  ramp2<-colorRampPalette(colors[2:3])(num.colors.in.palette*(1 - cutoff.fraction))
  
  return(c( ramp2, ramp1))
}


file_list_liver_top <- system(paste0("ls ", processed_dir, "/*/motifAnalysis/ded_top500_cov2_Liver.fa"), intern = TRUE)
file_list_heart_top <- system(paste0("ls ", processed_dir, "/*/motifAnalysis/ded_top500_cov2_Heart.fa"), intern = TRUE)

sp_list <- intersect(sapply(file_list_liver_top, function(x) strsplit(x, "/")[[1]][8]),
                           sapply(file_list_heart_top, function(x) strsplit(x, "/")[[1]][8]))


file_list_liver_add <- system(paste0("ls ", analysis_dir, "/03_motifAnalysis/diffMeth_additional_run/*/motifAnalysis/ded_top500_cov2_Liver.fa"), 
                              intern = T)

file_list_heart_add <- system(paste0("ls ", analysis_dir, "/03_motifAnalysis/diffMeth_additional_run/*/motifAnalysis/ded_top500_cov2_Heart.fa"), 
                              intern = T)

sp_list_add <- intersect(sapply(file_list_liver_add, function(x) strsplit(x, "/")[[1]][10]),
                     sapply(file_list_heart_add, function(x) strsplit(x, "/")[[1]][10]))
sp_list <- union(sp_list, sp_list_add)


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
sp_list <- sp_list[sp_list %in% sp_two_repl]

write.table(as.data.frame(sp_list), 
            file.path(Sys.getenv("CODEBASE"),"compEpi/meta/species_list_TF_analysis.txt"), 
            row.names = F, quote = F, col.names = F)


###############
#Functions#####
###############

##from the file list upload the AME scores and transform it into matrix TF/species
prepare_df <- function(sp_list, seq_set, bg_set){
  file_list <- paste0("screen/", sp_list, "/", seq_set, "_vs_", bg_set, "_top/ame.tsv")
  df_list <- sapply(file_list, read.csv, sep = "\t")
  print(length(file_list))
  col_n <-NCOL(df_list[[1]])
  col_names <- colnames(df_list[[1]])
  
  
  df_list_named <- sapply(seq_along(df_list), function(x) 
  {if(NROW(df_list[[x]])>3) {df <- df_list[[x]] %>% slice(1:(n()-3)) %>% mutate(species = sp_list[[x]],
                                                                                class = as.character(stats_annot[stats_annot$species == sp_list[[x]]]$color_class[[1]])) }
    else {df <- setNames(data.frame(matrix(ncol = col_n, nrow = 0)), col_names)}; 
    return(df)})
  
  df_full <- rbindlist(df_list_named[sapply(df_list_named, NROW)>0])
  df_mat <- df_full[, c("motif_alt_ID", "species", "adj_p.value")] %>% spread(species, adj_p.value) %>% replace(., is.na(.), 1)
  tfs <- df_mat$motif_alt_ID
  df_mat <- as.matrix(df_mat[, c(-1)])
  row.names(df_mat)<-tfs
  return(df_mat) 
}

## draw a small heatmap for top enriched sequences (for one class of species)
get_sub_heatmap<-function(df, name_hm, subdir){
  summ<-apply(df, 1, function(x) 1 - sum(x == 1)/NCOL(df))
  phm_list <-row.names(df)[summ>=0.1]
  return(phm_list)
}

## step - by - step analysis 
analyze_one <- function(seq_set, bg_set){
  subdir = file.path(wd, paste0("summary/", seq_set, "_vs_", bg_set))
  dir.create(subdir)
  
  df <- prepare_df(sp_list, seq_set = seq_set, bg_set = bg_set)
  
  ##drawing overall heatmap and creaing settings
  ####color palette
  cutoff.position <-  max(log(df[df!=1]))
  col<-makeColorRampPalette_2(c("#f0f0f0", "#e34a33", "#fdbb84"),
                              cutoff.position/min(log(df[df!=1])), 200)
  
  ##annotation
  annot_col = as.data.frame(sp_df$color_class)
  row.names(annot_col) <- row.names(sp_df)
  colnames(annot_col) <- c("class")
  my_colors <- list(class = class_colors)
  
  sp_list_ordered <- sp_full_order[sp_full_order %in% colnames(df)]
  print(length(sp_list_ordered))
  test_dev_off()
  pdf(paste0(subdir, "/Full_analysis.pdf"), width = 16, height = 10)
  pheatmap(log(df), color = col, annotation_col = annot_col,
           annotation_colors = my_colors, show_rownames = T, fontsize = 5, fontsize_row = 3)
  pheatmap(log(df[, sp_list_ordered]), color = col, annotation_col = annot_col,
           annotation_colors = my_colors, show_rownames = T, fontsize = 5, fontsize_row = 3, cluster_cols = F)
  dev.off()
  print("full saved")
  
  ## taking the subselection of TFs - the most enriched
  phm_all_names <- get_sub_heatmap(df[, sp_list_ordered], "all_species", subdir)
  print(phm_all_names)
  test_dev_off()
  pdf(paste0(subdir, "/top_tf_all.pdf"), width = 10, height = 5)
  pheatmap(log(df[phm_all_names, sp_list_ordered]), 
                  color = col, annotation_col = annot_col,
                  annotation_colors = my_colors, show_rownames = T, 
                  cluster_cols = F, fontsize_col = 5, main = "all species")
  dev.off()
  phm_all <- pheatmap(log(df[phm_all_names, sp_list_ordered]), 
                      cluster_cols = F, silent = T)
  
  phm_m_names <- get_sub_heatmap(df[, sp_list_ordered[annot_col[sp_list_ordered,] == "Mammalia"]],
                           "mammalia", subdir)
  print(phm_m_names)
  test_dev_off()
  pdf(paste0(subdir, "/top_tf_m.pdf"), width = 10, height = 5)
  pheatmap(log(df[phm_m_names, sp_list_ordered[annot_col[sp_list_ordered,] == "Mammalia"]]), 
           color = col, annotation_col = annot_col,
           annotation_colors = my_colors, show_rownames = T, 
           cluster_cols = F, fontsize_col = 5, main = "Mammalia")
  dev.off()
  
  
  phm_a_names <- get_sub_heatmap(df[, sp_list_ordered[annot_col[sp_list_ordered,] == "Aves"]],
                           "aves", subdir)
  if(length(phm_a_names) >0){
  print(phm_a_names)
  test_dev_off()
  pdf(paste0(subdir, "/top_tf_a.pdf"), width = 10, height = 5)
  pheatmap(log(df[phm_a_names, sp_list_ordered[annot_col[sp_list_ordered,] == "Aves"]]), 
           color = col, annotation_col = annot_col,
           annotation_colors = my_colors, show_rownames = T, 
           cluster_cols = F, fontsize_col = 5, main = "Aves")
  dev.off()
  }
  return(list(phm_all, df))
}

###no need to run, can just upload the files!
liver <- analyze_one(seq_set = "Liver", bg_set = "Heart")
liver_tree <- liver[[1]]$tree_rowliver_df <- liver[[2]]
rm(liver)
write.csv(liver_df, "summary/liver_df.csv", quote = F)

heart <- analyze_one(seq_set = "Heart", bg_set = "Liver")
heart_tree <- heart[[1]]$tree_row
heart_df <- heart[[2]]
rm(heart)
write.csv(heart_df, "summary/heart_df.csv", quote = F)

###CONTINUE HERE
liver_df <- read.csv("summary/liver_df.csv", row.names = 1)
heart_df <- read.csv("summary/heart_df.csv", row.names = 1)

b <- data.frame(TF = row.names(heart_df), as.data.frame(heart_df), row.names = NULL)
b <- gather(as.data.frame(b), species, adj.p.val, -1)
b$tissue <- "heart"

a <- data.frame(TF = row.names(liver_df), as.data.frame(liver_df), row.names = NULL)
a <- gather(as.data.frame(a), species, adj.p.val, -1)
a$tissue <- "liver"

ab <- rbind(a, b)
pdf("summary/pval_dist.pdf", width = 4, height =4)
ggplot(ab[ab$adj.p.val<1,], aes(x = log(adj.p.val))) + 
  geom_histogram(aes(fill = tissue), alpha = 0.5, position = "identity")
dev.off()
rm(ab)

##################################
##analysing the presence of TFBS##
##################################
## maybe modify threashold, f.e. p.val < 0.01 ?? <- doesn change a lot, checked

## as 1 are marked all the motifs, that are not discovered
heart_df_cp <- heart_df
heart_df_cp[heart_df_cp!=1] <- 0 

liver_df_cp <- liver_df

liver_df_cp[liver_df_cp!=1] <- 0 

tfs_common <- intersect(row.names(liver_df_cp), row.names(heart_df_cp))
species <- union(colnames(liver_df), colnames(heart_df))

ans <- as.data.frame(matrix(nrow = 0, ncol = length(species)))
colnames(ans) <- species

## pvalue == 1 -> non present in both -> "none"

for (tf in tfs_common){
  a <- liver_df_cp[tf,]
  b <- heart_df_cp[tf,]
  a <- c(a, setNames(rep(1, length(setdiff(species, names(a)))), setdiff(species, names(a))))
  b <- c(b, setNames(rep(1, length(setdiff(species, names(b)))), setdiff(species, names(b))))
  newrow <- list()
  
  for (sp in species){
    if (as.numeric(a[sp])){
      if(as.numeric(b[sp])) newrow[sp] <- "none"
      else newrow[sp] <- "heart"
    }else{
      if (as.numeric(b[sp])) newrow[sp] <- "liver"
      else newrow[sp] <- "both"
    }
  }
  ans <- rbind(ans, as.data.frame(t(newrow)))
}
ans <- apply(ans, 2, unlist)
row.names(ans) <- tfs_common

#tfs, detected only in heart
add_h <- heart_df_cp[setdiff(row.names(heart_df_cp),row.names(ans)), ]

add_h[add_h==0] <- "heart"
add_h[add_h==1] <- "none"

ans <- rbind.fill(as.data.frame(ans), as.data.frame(add_h))
row.names(ans) <- c(tfs_common, row.names(add_h))

##tfs, detected only in liver
add_l <- liver_df_cp[setdiff(row.names(liver_df_cp),row.names(ans)), ]
add_l[add_l==0] <- "liver"
add_l[add_l==1] <- "none"
ans <- rbind.fill(as.data.frame(ans), as.data.frame(add_l))
ans[is.na(ans)] <- "none"
row.names(ans) <- c(tfs_common, row.names(add_h), row.names(add_l))

ans_m <- melt(apply(ans, 2, unlist))
ans_m$Var2 <- factor(ans_m$Var2, levels = sp_full_order)


##for drawing annotations we need to switch to a heatmap with numerical values
##we assign to heart and liver values of 1 and -1, 
##so they are equally distanced from the 0 (being the center)

ans <- apply(ans, 2, unlist)
ans[ans=="none"] <- 0
ans[ans=="heart"] <- 1
ans[ans=="liver"] <- -1
ans <- apply(ans, 2, as.numeric)
ans <- ans[,colnames(ans) %in% sp_list ]
row.names(ans) <- c(tfs_common, row.names(add_h), row.names(add_l))
write.table(ans, "summary/tf_presence.csv", quote = F)

annot_col = as.data.frame(sp_df$color_class)

row.names(annot_col) <- row.names(sp_df)
colnames(annot_col) <- c("class")

ans[ans=="both"]
## we don´t have cases with "both", so we don´t have to consider it
#we need to be careful, which of the the colors indicate as which!
color=c(tissue_colors["L"], "lightgrey", tissue_colors["H"])
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

#creating a distinct color palette
n <- length(unique(my_tf_annot$Class_short))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Class_idst_col<-setNames(sample(col_vector, n), unique(my_tf_annot$Class_short))

dir.create("summary/final_heatmaps")

pdf("summary/TF_annot_by_Family.pdf", width = 15, height = 10)
ggplot(my_tf_annot, aes(x = Family_short, fill = Class_short)) + geom_bar(position="stack")+rotate_labels()+
  guides(fill = guide_legend(override.aes = list(size = 1), ncol = 4)) +
  theme(legend.title = element_text(size = 5), 
 legend.text  = element_text(size = 5),
 legend.key.size = unit(0.1, "lines"), legend.position = "bottom") + scale_fill_manual(values = Class_idst_col)
dev.off()

my_colors <- list(class = class_colors, Class_short = Class_idst_col)
test_dev_off()
pdf("summary/final_heatmaps/TFBS_presence.pdf", width = 15, height = 30)
pheatmap(ans[,sp_full_order[sp_full_order %in% colnames(ans)]], 
         cluster_cols = F,annotation_col = annot_col, annotation_colors = my_colors,
         fontsize = 5,  breaks = breaks[-5], color = color, 
         annotation_row = my_tf_annot[, 5, drop = F])
dev.off()

pdf("summary/final_heatmaps/TFBS_presence_turned.pdf", width = 28, height = 12)
pheatmap(t(ans[,sp_full_order[sp_full_order %in% colnames(ans)]]), 
         cluster_rows = F,annotation_row = annot_col, annotation_col = my_tf_annot[,  5, drop = F],
         annotation_colors = my_colors,
         fontsize = 5, annotation_legend = F , breaks = breaks[-5], color = color[-4])
dev.off()


##does the TFBS presence reconstruct the taxonomy?
test_dev_off()
pdf("summary/TFBS_presence_turned_taxonomy.pdf", width = 25, height = 10)
pheatmap(t(ans[,sp_full_order[sp_full_order %in% colnames(ans)]]), 
         annotation_row = annot_col, annotation_col = my_tf_annot[,  5, drop = F],
         annotation_colors = my_colors,
         fontsize = 7, annotation_legend = F , breaks = breaks[-5], color = color[-4])
dev.off()

###and lets filter so that TF shold be detected at 10 species (5% of the species) and speceis, 
###that have at least 10 signals (10%)
tf_to_keep <- names(rowSums(abs(ans))[rowSums(abs(ans)) >= 10])
ans_f <- ans[tf_to_keep,]
sp_to_keep <- names(colSums(abs(ans_f))[colSums(abs(ans_f)) >= 5])
ans_f <- ans_f[,sp_to_keep]
f <- unique(my_tf_annot[my_tf_annot$Name %in% row.names(ans_f),]$Family_short)
my_colors <- list(class = class_colors, Class_short = Class_idst_col[unique(my_tf_annot[my_tf_annot$Name %in% row.names(ans_f),]$Class_short)])
## filter species, that contain signals from at least 10 TFs
test_dev_off()
pdf("summary/final_heatmaps/TFBS_presence_turned_taxonomy_filtered.pdf", width = 15, height = 10)
pheatmap(t(ans_f), 
         annotation_row = annot_col, annotation_col = my_tf_annot[,  5, drop = F],
         annotation_colors = my_colors,
         fontsize = 7, breaks = breaks, color = color)
dev.off()

pdf("summary/final_heatmaps/TFBS_presence_taxonomy_filtered.pdf", width = 15, height = 10)
pheatmap(ans_f, 
         annotation_col = annot_col, annotation_row = my_tf_annot[,  5, drop = F],
         annotation_colors = my_colors,
         fontsize = 7, breaks = breaks, color = color)
dev.off()

ph <- pheatmap(t(ans_f))

tf_tree <- ph$tree_col
tf_tree_dend <- as.dendrogram(tf_tree)

library(dendextend)

pdf("summary/TF_clusters.pdf",  height = 5, width = 12)
tf_tree_dend %>% set("labels_col", value = c("green", "blue", "red", "black", "orange"), k=5) %>% 
  plot(main = "Color labels per cluster in TFs")
dev.off()

cl <- cutree(tf_tree, k = 5)
cl_df <- as.data.frame(sapply(unique(cl), function(x) paste(names(cl[cl == x]), collapse = ", ")))
colnames(cl_df) <- "clusters"
write.table(cl_df,"summary/TF_clusters_detailed.csv", sep = "\t", quote = F)


####is there correlation between TFBS count and the frequency of differentially meth. regions
ans <- read.csv("summary/tf_presence.csv", sep = " ", row.names = 1)
head(ans)

dir.create("diffRegFr/summary", recursive = T)

read_diff <- function(cov_tr){
  file_list_diffM <-system(paste0("ls ", analysis_dir,
                "/03_motifAnalysis/03.7_AME/diffRegFr/*/diffMethdedRefFreq_", cov_tr, ".csv"), intern = T)
  diffM_list <- sapply( file_list_diffM, read.csv, sep = " ", simplify = F)
  diffM <- rbindlist(diffM_list)
  diffM$cov_tr <- cov_tr
  return(diffM)
}

diffM_list <- lapply(c(2, 5, 10, 15, 30), read_diff)
diffM <- rbindlist(diffM_list)

diffM_2 <- read_diff(2)

diffM_total_m <- melt(diffM, measure.vars = c("Liver_all_cpgs", "Heart_all_cpgs"))
head(diffM_total_m)

#by coverage
ggplot(diffM_total_m, aes(x = value, fill = variable)) + geom_histogram(alpha = 0.5, position = "identity") +
  facet_wrap(~cov_tr, ncol = 1, scales = "free") + 
  scale_fill_manual(values = c("Heart_all_cpgs" = tissue_colors[["H"]], "Liver_all_cpgs" = tissue_colors[["L"]]))+ 
  geom_vline(xintercept = 500, linetype = "dashed", alpha = 0.5)
ggsave("diffRegFr/summary/count_of_diffMetFr_by_cov.pdf", width = 4, height = 15)

#by tissue
ggplot(diffM_total_m, aes(x = value, fill = as.factor(cov_tr))) + geom_histogram(alpha = 0.5, position = "identity") +
  facet_wrap(~variable, ncol = 1) + geom_vline(xintercept = 500, linetype = "dashed", alpha = 0.5)
ggsave("diffRegFr/summary/count_of_diffMetFr_by_tissue.pdf", width = 4, height = 5)

#ggplot(diffM, aes(x = cov_tr, y = Liver_all_cpgs, color = SP)) + geom_line() + 
#  theme(legend.position = "None")

N_of_TFs <- colSums(abs(ans))
diffM_2$N_of_TF <- sapply(as.character(diffM_2$SP), function(x) N_of_TFs[x])

diffM_2$Liver_fr <- diffM_2$Liver / diffM_2$Total
diffM_2$Heart_fr <- diffM_2$Heart / diffM_2$Total

diffM_2_m <- melt(diffM_2[complete.cases(diffM_2),c("SP", "Liver_fr", "Heart_fr", "N_of_TF")], 
                  id.vars= c("SP", "N_of_TF"))

ggplot(diffM_2_m, aes(x = value, y = N_of_TF, color = variable)) + geom_point(alpha = 0.5)+ 
      scale_color_manual(values = c("Heart_fr" = tissue_colors[["H"]], "Liver_fr" = tissue_colors[["L"]])) +
      xlab("count of diff.meth.ded.Ref")+
      geom_text_repel(data = diffM_2_m[diffM_2_m$value >  0.1, ], aes(x = value, y = N_of_TF, label = SP)) 
ggsave("diffRegFr/summary/frequency_check_count.pdf", width = 5, height = 5)      



###STOP THERE######
####all further analysis wasn´t included in the figure

##now let calculate how many times each TF Family was detected was detected for each species:
ans_df <- as.data.frame(ans)
ans_df$Name <- row.names(ans)
ans_df <- ans_df %>% gather("species", "presence",-Name)
ans_df <- left_join(ans_df, my_tf_annot[, c("Name", "Family_short")])
ans_df_by_class <- ans_df %>%
  filter(presence!=0)%>%
  group_by(species, Family_short) %>%
  dplyr::summarize(presence_by_class = n_distinct(presence), tissue = max(presence))

non_uniuqe_families <- unique(ans_df_by_class[ans_df_by_class$presence_by_class == 2,]$Family_short)

##first, TF families that don have ambigous tissue - specificity within one species!
ans_df_by_class_unique <- ans_df_by_class[!ans_df_by_class$Family_short %in% non_uniuqe_families,]%>%
  select(-presence_by_class) %>%
  spread(Family_short, tissue)%>%
  replace(., is.na(.), 0)

ans_df_by_class_unique <- as.data.frame(ans_df_by_class_unique)
row.names(ans_df_by_class_unique) <-ans_df_by_class_unique$species

col_annot <- unique(my_tf_annot[!my_tf_annot$Family_short %in% non_uniuqe_families, c(4,5)])
col_annot$Family_short <- as.character(col_annot$Family_short)
row.names(col_annot) <- col_annot$Family_short
my_colors <- list(class = class_colors, Class_short = Class_idst_col[unique(col_annot$Class_short)])
pdf("summary/TFBS_presence_by_family_taxonomy.pdf", width = 15, height = 25)
pheatmap(ans_df_by_class_unique[, c(-1)], 
         annotation_row = annot_col, 
         annotation_col = col_annot[, c(-1), drop = F],
         annotation_colors = my_colors,
         fontsize = 5, breaks = breaks[-5], color = color[-4], treeheight_row = 0)
dev.off()


test_dev_off()
pdf("summary/TFBS_presence_by_family.pdf", width = 15, height = 25)
pheatmap(ans_df_by_class_unique[sp_full_order[sp_full_order %in% row.names(ans_df_by_class_unique)], c(-1)], 
         annotation_row = annot_col, 
         annotation_col = col_annot[, c(-1), drop = F],
         annotation_colors = my_colors,
         cluster_rows = F,
         fontsize = 7, breaks = breaks[-5], color = color[-4], treeheight_row = 0)
dev.off()


###Those, that change the specificity within TF class for one species
ans_df_non_unique <- ans_df[ans_df$Family_short %in% non_uniuqe_families,] %>%
  select(-Family_short) %>%
  spread(Name, presence) %>%
  replace(., is.na(.), 0)

ans_df_non_unique <- as.data.frame(ans_df_non_unique)
row.names(ans_df_non_unique) <- ans_df_non_unique$species
ans_df_non_unique <- ans_df_non_unique[sp_full_order[sp_full_order %in% row.names(ans_df_non_unique)], c(-1)]
ans_df_non_unique[,!apply(ans_df_non_unique, 2, is.numeric)]
col_annot <- my_tf_annot[my_tf_annot$Family_short %in% non_uniuqe_families, c(4,5)]
my_colors <- list(class = class_colors, Class_short = Class_idst_col[unique(col_annot$Class_short)])
##dropping the species with no presence at all
ans_df_non_unique <- ans_df_non_unique[rowSums(ans_df_non_unique) > 0,]
test_dev_off()

pdf("summary/TFBS_presence_nonunique_in_family.pdf", width = 15, height = 10)
pheatmap(t(ans_df_non_unique), 
        annotation_col = annot_col, 
        annotation_row = col_annot,
        annotation_colors = my_colors,
        cluster_cols = F,
        fontsize = 7, breaks = breaks[-5], color = color[-4])
dev.off()

head(ans_m)
colnames(ans_m) <- c("TF", "species", "detected")
ans_m <- left_join(ans_m, sp_df, by = "species")
ans_m_class <- left_join(ans_m, my_tf_annot[, c("Name", "Family_short", "Class_short")], by = c("TF" = "Name"))

n_by_class <- sp_df[sp_df$species %in% sp_list,] %>%
  group_by(color_class) %>%
  dplyr::summarize(n = dplyr::n()) 
class_count <- setNames(n_by_class$n, n_by_class$color_class)

ans_m_class_count <- ans_m_class %>%
  filter(!Family_short %in% non_uniuqe_families) %>%
  group_by(Family_short, color_class, detected) %>%
  dplyr::summarize(n = n()) %>%
  spread(detected, n) %>%
  replace_na(list(liver=0, heart=0)) %>%
  select(-none)

#ans_m_class_count <- setDT(melt(ans_m_class_count, measure.vars = c("heart", "liver")))

ans_list <- data.frame(matrix(nrow= 0, ncol = 4))
colnames(ans_list) <- c("Family_short", "tissue", "presence", "color_class")
for(CLASS in names(class_count)){
  sub_df <- ans_m_class_count[ans_m_class_count$color_class==CLASS, -c(2)]
  
  sub_df <-sub_df %>%
    filter(liver >= 0.1*class_count[CLASS] | heart >= 0.1*class_count[CLASS]) %>%
    gather("tissue", "presence", -Family_short)
  sub_df$color_class <-CLASS
  #plot_class <-ggplot(sub_df, aes(x = Family_short, y = presence, fill = tissue)) + 
    #geom_bar(stat = "identity", position = "dodge") + rotate_labels() + ggtitle(CLASS) + 
    #scale_fill_manual(values = c("heart" = "#00a783", "liver" = "#D1B0FF" )) #+ theme(legend.position = "none")
  #ggsave(paste0("summary/TFBS_presence_heart_vs_liver_", CLASS, ".pdf" ), plot_class, width = 6, height = 5)
  ans_list <-rbind(ans_list, sub_df)
}

family_count <- my_tf_annot %>%
  group_by(Family_short) %>%
  dplyr::summarise(n_in_family = n())

ans_list <- left_join(ans_list, as.data.frame(family_count), by = "Family_short")

ans_list$color_class <- factor(ans_list$color_class, levels = names(class_colors))

ggplot(ans_list, aes(x = Family_short, y = presence, fill = tissue)) + 
  geom_bar(stat = "identity", position = "dodge") + rotate_labels() + 
  facet_wrap(~color_class, ncol = 1, scale="free_y") + 
  scale_fill_manual(values = c("heart" = "#00a783", "liver" = "#D1B0FF" )) +
  labs(y="count of TFBS detected")
ggsave(paste0("summary/TFBS_presence_heart_vs_liver.pdf" ), width = 8, height = 15)

sp_order <- sp_full_order[sp_full_order %in% colnames(ans)]
draw_class_hm<- function(CLASS, H, W){
  tf_list <- as.character(my_tf_annot[my_tf_annot$Class_short==CLASS,]$Name)
  print(tf_list)
  ans_sub <- ans[tf_list,sp_order]
  ans_sub <- ans_sub[,colSums(ans_sub!=0)>0]
  my_colors <- list(class = class_colors)
  test_dev_off()
  pdf(paste0("summary/presence_by_tf_class/", paste(strsplit(CLASS, " / ")[[1]], collapse = " "),".pdf"), width = W, height = H)
  pheatmap(ans_sub, main = CLASS,
           annotation_col = annot_col,
           annotation_row = my_tf_annot[, c(4), drop = F],
           annotation_colors = my_colors, fontsize = 5,
           annotation_legend = T , breaks = breaks[-5], color = color[-4])
  dev.off()
}

draw_class_hm("Basic helix-loop-helix factors (bHLH)", 5,5)
draw_class_hm("C2H2 zinc finger factors", 7, 8)
draw_class_hm("Fork head / winged helix factors", 7, 8)
draw_class_hm("Homeo domain factors", 7, 8)
draw_class_hm("Nuclear receptors with C4 zinc fingers", 6,8)
draw_class_hm("Tryptophan cluster factors", 7,8)

## now lets plot presence by TF in interesting families
tf_list <- as.character(my_tf_annot[my_tf_annot$Family_short=="PAS domain factors",]$Name)
sp_order <- sp_full_order[sp_full_order %in% colnames(ans)]
ans_sub <- ans[tf_list,sp_order]
ans_sub <- ans_sub[,colSums(ans_sub!=0)>0]
my_colors <- list(class = class_colors, 
      Class_short = Class_idst_col[unique(my_tf_annot[my_tf_annot$Family_short=="PAS domain factors",]$Class_short)])
dir.create("summary/presence_by_tf_class")
test_dev_off()
pdf("summary/presence_by_family/PAS_domain_factors.pdf", width = 5, height = 5)
pheatmap(t(ans_sub), cluster_rows = F, 
         annotation_row = annot_col,
         annotation_col = my_tf_annot[, c(5), drop = F],
         annotation_colors = my_colors, fontsize = 5,
         fontsize_row = 10, fontsize_col = 10,annotation_legend = T , breaks = breaks[-5], color = color[-4])
dev.off()
ans_summary <- ans_m %>%
  group_by(TF, color_class, detected) %>%
  dplyr::summarize(n = n()) %>%
  spread(detected, n) %>%
  replace_na(list(liver=0, heart=0))

ans_summary$n_of_species <- ans_summary$heart+ans_summary$liver + ans_summary$none
ans_summary$liver_f <- ans_summary$liver/ans_summary$n_of_species
ans_summary$heart_f <- ans_summary$heart/ans_summary$n_of_species
ans_summary <- as.data.frame(ans_summary)

pheatmap(ans_summary[ans_summary$color_class=="Actinopteri", c("liver_f", "heart_f")], cluster_cols = F)

ans_summary_total <- ans_summary %>%
  group_by(TF) %>%
  dplyr::summarise(total_liver = sum(liver), total_heart = sum(heart)) %>%
  filter(total_liver!=0 | total_heart != 0)

  ans_summary$freq <- ans_summary[, tissue]/ans_summary$n_of_species
  
  ans_only_by_class <- ans_summary[ans_summary$TF %in% specific_tfs, c("TF", "color_class", "freq")] %>%
    spread(color_class, freq)
  
  ans_only_by_class <- as.data.frame(ans_only_by_class)
  row.names(ans_only_by_class) <- ans_only_by_class$TF
  ## filtering
  ans_only_by_class_f <- ans_only_by_class[, c(-1)]
  ans_only_by_class_f[ans_only_by_class_f<0.05] <- 0
  ans_only_by_class_f <- ans_only_by_class_f[rowSums(ans_only_by_class_f) > 0,]
  cols = colorRampPalette(c("lightgrey", "red"))(200)
  test_dev_off()
  pdf(paste0("detected_only_in_", tissue, ".pdf"), width = 5, height = 10)
  pheatmap(ans_only_by_class_f, color = cols)
  dev.off()
  return(row.names(ans_only_by_class_f))
}

###ALWAYS DETECTED IN HEART###
liver_sp <- tissue_specific_tf("heart", ans_summary, ans_summary_total[ans_summary_total$total_liver==0,]$TF)
###ALWAYS DETECTED IN LIVER###
heart_sp <- tissue_specific_tf("liver", ans_summary, ans_summary_total[ans_summary_total$total_heart==0,]$TF)

##now lets take a look at the TFs that change their specificity:
non_sp_tfs <- ans_summary_total[ans_summary_total$total_heart!=0 & ans_summary_total$total_liver!=0,]$TF
length(non_sp_tfs)
ans_non_sp <- ans_summary[ans_summary$TF %in% non_sp_tfs, ]

ans_non_sp$presence <- 0 #in none
ans_non_sp$presence[ans_non_sp$heart>0 & ans_non_sp$liver==0 ] <- 1 #only in heart
ans_non_sp$presence[ans_non_sp$heart==0 & ans_non_sp$liver>0 ] <- 2 #only in liver
ans_non_sp$presence[ans_non_sp$heart>0 & ans_non_sp$liver>0 ] <- 3 #in both

pdf("TF_presence_by_group.pdf", width = 5, height = 5)
ggplot(ans_non_sp, aes(x = presence, fill = color_class)) + geom_bar(position = "dodge") + 
  scale_fill_manual(values = class_colors) + 
  scale_x_continuous(breaks=c(0,1,2,3), labels=c("none", "heart", "liver", "both")) +
  ggtitle("TF detection in tissue by taxonomy group")
dev.off()
  
ans_non_sp_m <- spread(ans_non_sp[, c("TF", "color_class", "presence")], color_class, presence)
rownames(ans_non_sp_m) <- ans_non_sp_m$TF
new_colnames<-c("TF", as.character(sapply(colnames(ans_non_sp_m[, c(-1)]), function(x) class_short[[x]])))
colnames(ans_non_sp_m) <- new_colnames


##full evoulutionary picture:
out <- pheatmap(ans[non_sp_tfs,sp_full_order[sp_full_order %in% colnames(ans)]], 
                cluster_cols = F,annotation_col = annot_col, annotation_colors = my_colors,
                fontsize = 5, annotation_legend = F , breaks = breaks[-5], color = color[-4])

pdf("summary/ambigous.pdf", width = 5, height = 15)
pheatmap(ans[non_sp_tfs,sp_full_order[sp_full_order %in% colnames(ans)]], 
         cluster_cols = F,annotation_col = annot_col, annotation_colors = my_colors,
         fontsize = 5, annotation_legend = F , breaks = breaks[-5], color = color[-4])
dev.off()

pdf("summary/ambigous_by_class.pdf", width = 4, height = 18)
pheatmap(ans_non_sp_m[non_sp_tfs, c(-1)], cluster_cols = F,
         cluster_rows = out$tree_row, breaks = breaks, color = color, fontsize_row = 5)
dev.off()

##now let??s integrate the known liver - diff. expressed genes from the Mouse Cell Atlas
df <- read.csv("../Liver_genes_MCA.csv", sep = ";")
check_overlap <- function(sub_df, sp){
  colnames(sub_df) <-c("Cell", "FACS", "10x", "MicrowellSeq")
  unique_genes <-union(union(sub_df$FACS, sub_df$`10x`), sub_df$microwellSeq)
  return(union(intersect(unique_genes, sp), intersect(toupper(unique_genes), sp)))
}

liver_sp_f <- ans_summary_total[ans_summary_total$total_liver==0,]$TF
liver_confirmed<-union(union(check_overlap(df[,c(1:4)], liver_sp), 
                             check_overlap(df[,c(10:14)], liver_sp)), 
      check_overlap(df[,c(19:22)], liver_sp))

##sanity - check - heart-sp
heart_sp <- ans_summary_total[ans_summary_total$total_heart==0,]$TF
check_overlap(df[,c(1:4)], heart_sp)
heart_sp



