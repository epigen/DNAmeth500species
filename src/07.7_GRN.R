source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(biomaRt)

#directory setup
wd = file.path(analysis_dir, "03_motifAnalysis/03.7_AME")
setwd(wd)
outdir = file.path(wd, "GRN")
dir.create(outdir)


library(stringr)
library(topGO)

##TF data is downloaded from the TRRUST database
TF_Target <- read.table("trrust_rawdata.human.tsv")
head(TF_Target)
colnames(TF_Target) <- c("TF", "target", "action", "pmid")

gene_clusters <- read.table("summary/clusters_and_selex.csv", sep = " ", 
                            header = 1)

##clusters are coming from the 03.71, with cluster 1 responsible for genes, 
##hypomethylated in liver/hypermethylated in heart, 2 - other way around.
##In case the TF prefers the methylated BS (activity is MethylPlus), it is 
##active in the hypermethylated state, otherwise - in hypomethylated

gene_clusters$active_in <- c("L")
gene_clusters[gene_clusters$cl == 1 & 
                gene_clusters$Call_upd == "MethylPlus",]$active_in <-"H"
gene_clusters[gene_clusters$cl == 2 & 
                gene_clusters$Call_upd != "MethylPlus",]$active_in <- "H"
gene_clusters$TF_upd <- as.character(gene_clusters$TF_upd)

##identifying the targets of the TFs, that are active in heart/liber
H_targets <- unique(TF_Target[TF_Target$TF %in% 
                  gene_clusters[gene_clusters$active_in == "H",]$TF_upd, ]$target)
length(H_targets)

L_targets <- unique(TF_Target[TF_Target$TF %in% 
                  gene_clusters[gene_clusters$active_in == "L",]$TF_upd, ]$target)
length(L_targets)

#all targets combined
all_targets <- union(L_targets, H_targets)

#intersection of targets
common_targets <- intersect(H_targets, L_targets)

##all the genes, that are targeted by any of the transciption factors and their 
##regulators

#edges
sub_tf_t <- TF_Target[TF_Target$TF %in% gene_clusters$TF_upd, c(-4)]

colnames(sub_tf_t) <- c("source", "target", "interaction")
write.table(sub_tf_t, "GRN/edges_full.tsv", 
            sep = "\t", quote = F, row.names = F)

##nodes
nodes <- data.frame(name = unique(union(sub_tf_t$source, sub_tf_t$target)), 
                   action = c("none"))
nodes$action <- as.character(nodes$action)
nodes$name <- as.character(nodes$name)
nodes[nodes$name %in% 
        gene_clusters[gene_clusters$active_in == "H",]$TF_upd, ]$action <- "H"

nodes[nodes$name %in% 
        gene_clusters[gene_clusters$active_in == "L",]$TF_upd, ]$action <- "L"


nodes$role <- c("target")
nodes$role <- as.character(nodes$role)
nodes[nodes$name %in% TF_Target$TF, ]$role <- "TF"

write.table(nodes, "GRN/nodes_full.tsv", 
            sep = "\t", quote = F, row.names = F)


##subselecting and filtering - only the common actions and known interactions

sub_tf_s <- sub_tf_t[sub_tf_t$target %in% common_targets,]
sub_tf_s <- sub_tf_s[sub_tf_s$interaction != "Unknown", ]

#dropping ambigous interactions
to_drop_dupl <- as.integer(row.names(sub_tf_s[duplicated(sub_tf_s[, 
                                                  c("source", "target")]),]))
to_drop_dupl <- c(to_drop_dupl - 1, to_drop_dupl)
sub_tf_s <- sub_tf_s[!row.names(sub_tf_s) %in% as.character(to_drop_dupl),]

write.table(sub_tf_s, "GRN/edges_short.tsv", 
            sep = "\t", quote = F, row.names = F)
nodes_f <- nodes[nodes$name %in% unique(union(sub_tf_s$source, sub_tf_s$target)), ]
write.table(nodes_f, "GRN/nodes_short.tsv", 
            sep = "\t", quote = F, row.names = F)


features<- t(data.frame(c("TF", "circle"), c("target", "square"), 
             c("H", tissue_colors["H"]), c("L", tissue_colors["L"])))
colnames(features) <- c("feature", "value")
write.table(features, "GRN/feature_annot.tsv", 
            sep = "\t", quote = F, row.names = F)

#library(igraph)
#g <- graph_from_data_frame(sub_tf_t[sub_tf_t$target %in% common_targets,], directed = T, vertices = tfs)
#plot(g)
#color_action <- list("Activation" = "green", "Repression" = "red", "Unknown" = "grey")
#edge.color <- color_action[E(g)$action]
#E(g)$arrow.size <- .4

#plot(g)

#edge.col <- Reduce(c, color_action[E(g)$action])

#V(g)$size <- 5
#V(g)$vertex.label.cex <- 0.5
#V(g)$vertex.label.color <- "black"
#V(g)$vertex.label.dist <- 1
#E(g)$color <- edge.col
#plot(g, vertex.label.color = "black", label.dist = 1, label.cex = 0.5)

#pdf("GRN/igraph_GRN.pdf", height = 10, width = 15)
#print(plot(g, vertex.label.color = "black", vertex.label.dist = 0.65, 
  #         vertex.label.cex = 0.7, vertex.label.degree = - pi/2))
#dev.off()

