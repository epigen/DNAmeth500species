#!/bin/env Rscript
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

library(ggdendro)
library(dendextend)
library(mclust)
library(circlize)

wd = file.path(analysis_dir, "04_phylogeny")
dir.create(wd)
setwd(wd)



##useful annotations
##from here we want the coveredCpGs
stats<-read.csv(file.path(analysis_dir, "01_basicStats", "01.10_all_seq_features","feature_summary_filtered.tsv"), sep = "\t")
row.names(stats)<-stats$species

##number of dedRef fragments
dedRefsize<-read.csv(file.path(analysis_dir, "02_vizStats", "02.3_stats_detail", "number_of_dedRef_fragments.csv"), row.names = 1)
##color mapping for both of them
col_fun = colorRamp2(c(min(dedRefsize$counts), max(dedRefsize$counts)),  c( "#efedf5", "#756bb1"))
col_fun_2 = colorRamp2(c(min(stats$coveredCpGs), max(stats$coveredCpGs)), c("#e5f5e0", "#31a354"))

annot_blast <- stats_annot %>%
  group_by(species) %>%
  summarize(rank_blast = mean(max_lowest_common))
annot_blast <- as.data.frame(annot_blast)
row.names(annot_blast) <- annot_blast$species
annot_blast[is.infinite(annot_blast$rank_blast),]$rank_blast <- 0

###3mers
k <- 3
kmer_table<-read.csv(paste0(analysis_dir, "/01_basicStats",
                            "/01.5_kmercount/summary/kmer0", k, ".csv"), row.names=1)
kmer_table<-t(kmer_table)
kmer_table <- kmer_table[unique(stats_annot$species),]
M <- matrix(nr = 0, nc = NROW(kmer_table))

#annotation (same order as in M) - not used in the final figure

#ha <- HeatmapAnnotation(class = sp_df[row.names(kmer_table),"color_class" ], 
#                        CpGs = stats[row.names(kmer_table), "coveredCpGs"],
#                        dedRef = dedRefsize[row.names(kmer_table), "counts"],
#                        blast = annot_blast[row.names(kmer_table), "rank_blast"],
#                        col = list(class = class_colors, 
#                                   dedRef = col_fun, CpGs = col_fun_2), 
#                        border = T)
sp_df <- as.data.frame(sp_df)
row.names(sp_df) <- sp_df$species
ha <- HeatmapAnnotation(class = sp_df[row.names(kmer_table),"color_class" ], col = list(class = class_colors), border = T)

clusters<-hclust(dist(kmer_table))
dd <- as.dendrogram(clusters)
print(Heatmap(M, cluster_columns = dd, top_annotation = ha, show_column_names = F))
which(sapply(labels(dd), 
             function(x) as.character(sp_df[x,]$color_class)) == "Actinopteri")
             


dd1 <- dendextend::rotate(dd, c(580:328,1:327))
Heatmap(M, cluster_columns = dd1, top_annotation = ha, show_column_names = F)

##swapping the all red and mixed branches:

which(sapply(labels(dd1), 
                 function(x) as.character(sp_df[x,]$color_class)) == "Actinopteri")
which(sapply(labels(dd1), 
             function(x) as.character(sp_df[x,]$color_class)) == "Amphibia")
##1st 40-101 2nd 102 - 252
dd2 <- dendextend::rotate(dd1, c(1:39, 255:105, 40:104, 256:580))
Heatmap(M, cluster_columns = dd2, top_annotation = ha, show_column_names = F)
             
##swapping the reptilia and aves, forming orange cluster
which(sapply(labels(dd2), 
             function(x) as.character(sp_df[x,]$color_class)) == "Reptilia")
which(sapply(labels(dd2), 
             function(x) as.character(sp_df[x,]$color_class)) == "Aves")
## start 312-387, 388-476
dd3 <- dendextend::rotate(dd2, c(1:311, 388:476,387:312,546:580,545:478,477))
Heatmap(M, cluster_columns = dd3, top_annotation = ha, show_column_names = F)

pdf(paste0("dendragram_rotated_", k, "_no_annot.pdf"), width = 15, height = 6)
print(Heatmap(M, cluster_columns = dd2, top_annotation = ha, show_column_names = F))
dev.off()

###for 6-mers
k <- 6
kmer_table<-read.csv(paste0(analysis_dir, 
                            "/99.4_kmer_count_filtered/summary/kmer0", k, ".csv"), row.names=1)
kmer_table<-t(kmer_table)
M <- matrix(nr = 0, nc = NROW(kmer_table))
ha <- HeatmapAnnotation(class = sp_df[row.names(kmer_table),"color_class" ], 
                        col = list(class = class_colors), border = T)
clusters<-hclust(dist(kmer_table))
dd <- as.dendrogram(clusters)
Heatmap(M, cluster_columns = dd, top_annotation = ha, show_column_names = F)
which(sapply(labels(dd), 
             function(x) as.character(sp_df[x,]$color_class)) == "Invertebrata")
##start 324 - rotating the invertebrata cluster
dd1 <- dendextend::rotate(dd, c(324:583,323:2,1))
Heatmap(M, cluster_columns = dd1, top_annotation = ha, show_column_names = F)
## rotating the middle Aves brunch, so that the actinopteri is closer to others
which(sapply(labels(dd1), 
             function(x) as.character(sp_df[x,]$color_class)) == "Aves")
##start 262 ends 373 (372 - act)
dd2 <- dendextend::rotate(dd1, c(1:261,373:262,374:583))
Heatmap(M, cluster_columns = dd2, top_annotation = ha, show_column_names = F)

## swapping the two ending branches: 
which(sapply(labels(dd2), 
             function(x) as.character(sp_df[x,]$color_class)) == "Aves")
which(sapply(labels(dd2), 
             function(x) as.character(sp_df[x,]$color_class)) == "Reptilia")
##2nd: 528-582
##1st: 374 - 527
dd3 <- dendextend::rotate(dd2, c(1:373,528:582,374:527,583))
Heatmap(M, cluster_columns = dd3, top_annotation = ha, show_column_names = F)
## roating the blue-purple-red branch in thr middle
which(sapply(labels(dd3), 
             function(x) as.character(sp_df[x,]$color_class)) == "Actinopteri")
which(sapply(labels(dd3), 
             function(x) as.character(sp_df[x,]$color_class)) == "Mammalia")
##ends: 260, starts 168
dd4 <- dendextend::rotate(dd3, c(1:167,260:168,261:583))
Heatmap(M, cluster_columns = dd4, top_annotation = ha, show_column_names = F)

pdf(paste0("dendragram_rotated_", k, "_no_annot.pdf"), width = 15, height = 6)
print(Heatmap(M, cluster_columns = dd4, top_annotation = ha, show_column_names = F))
dev.off()

## easier generation of the phylogeny, but without nice arragement
make_dendr<-function(k){
  #reading the data
  kmer_table<-read.csv(paste0(analysis_dir, "/01.5_kmercount/summary/kmer0", k, ".csv"), row.names=1)
  kmer_table<-t(kmer_table)
  
  #clustering based on frequency of kmers
  clusters<-hclust(dist(kmer_table))
  #turning to the dendrogramm (passed to ComplexHeatmap)
  dd.row <- as.dendrogram(clusters)
  
  #empty matrix, nameing of the columns should be same as in the object that was 
  #clustered (so row.names(kmer.table))
  M <- matrix(nr = 0, nc = NROW(kmer_table))
  colnames(M) <- row.names(kmer_table)
  
  #annotation (same order as in M)
  ha <- HeatmapAnnotation(class = sp_df[row.names(kmer_table),"color_class" ], 
                          CpGs = stats[row.names(kmer_table), "coveredCpGs"],
                          dedRef = dedRefsize[row.names(kmer_table), "counts"], 
                          col = list(class = class_colors, dedRef = col_fun, CpGs = col_fun_2), 
                          border = T)
  
  #saving
  pdf(paste0("dendragram_", k, ".pdf"), width = 15, height = 6)
  print(Heatmap(M, cluster_columns = dd1, top_annotation = ha, show_column_names = F))
  dev.off()  
}

make_dendr(6)
make_dendr(5)
make_dendr(4)
make_dendr(3)
make_dendr(2)
make_dendr(1)

