source(file.path(Sys.getenv("CODEBASE"),"DNAmeth300species/src/02.60_predicatbility_init.R"))
#for the further analysis we are making a table with the key data  from the sequences:
#important statistics
summary_conv<-stats_annot %>% 
  filter(conversion_type=="converted") %>%
  group_by(species) %>%
  summarize(mean_meth=mean(CpG_meth),coveredCpGs=mean(coveredCpGs), color_class=first(color_class), tier=mean(mean_qual_tier))

summary_unconv<-stats_annot %>% 
  filter(conversion_type=="unconverted") %>%
  group_by(species) %>%
  summarize(mean_meth_unconv=mean(CpG_meth),coveredCpGs_unconv=mean(coveredCpGs))

summary<-setDT(full_join(summary_conv, summary_unconv, by="species"))

#plotting converted vs unconverted
pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/coveredCpGs_conv_vs_unconv.pdf"))
draw_correlation(summary$coveredCpGs, summary$coveredCpGs_unconv, "covered_CpGs_converted", "covered_CpGs_unconverted")
dev.off()

pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/mean_meth_conv_vs_unconv.pdf"))
draw_correlation(summary$mean_meth, summary$mean_meth_unconv, "mean_meth_converted", "mean_meth_unconverted")
dev.off()

pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/mean_meth_quality.pdf"))
draw_correlation(summary$mean_meth, summary$tier, "mean_meth", "quality")
dev.off()
rm(summary_conv)                  
rm(summary_unconv)


kmer1<-t(read.csv(file.path(analysis_dir, "99.4_kmer_count_filtered/summary/kmer01.csv"), row.names=1))
kmer2<-t(read.csv(file.path(analysis_dir, "99.4_kmer_count_filtered/summary/kmer02.csv"), row.names=1))
kmer3<-t(read.csv(file.path(analysis_dir, "99.4_kmer_count_filtered/summary/kmer03.csv"), row.names=1))

frequencies<-merge(kmer1, kmer2,  by=0)
row.names(frequencies)<-frequencies$Row.names
frequencies<-merge(frequencies[, -1], kmer3, by=0)

rm(kmer1)
rm(kmer2)
rm(kmer3)

colnames(frequencies)[1]<-"species"

summary<-merge(summary, frequencies, by="species")

#calculating the sequence - derived features
colnames(summary)[which(colnames(summary)=="CG")]<-"CpG_O"
summary$CG_freq<-summary$G+summary$C
summary$CpG_OE<-summary$CpG_O/(summary$C*summary$G)


##uploading CpG - island frequencies:
stats_files=system(paste0("ls ", analysis_dir, "/99.5_CpG_islands_per_ref/*.csv"), intern=TRUE)
file_list <- lapply(stats_files, read.csv)
CpG_cr <- Reduce(full_join, file_list)

colnames(CpG_cr)<-c("species", "Gardiner_Garden", "Takai_Jones")

summary<-merge(summary, CpG_cr, by="species")

##adding the self AUC
auc_stats<-read.csv(file.path(analysis_dir, "02_predict_meth/02.1_within_species/summary/all_stats.csv"))
summary<-merge(summary, auc_stats[,c("Species", "k", "AUC")], by.x="species", by.y="Species", all.x=T)


pdf(file.path(analysis_dir, "02_predict_meth/02.6_predictability/feature_summary.pdf"), width=4, height=4*(length(summary)-2))
ggplot(melt(summary, id.vars = c("species", "color_class")), aes(x=value, fill=color_class)) + 
  geom_histogram() + facet_wrap(~variable, scales = "free", ncol = 1) + 
  scale_fill_manual(values=class_colors) + theme(legend.position = "None")
dev.off()

my_wt(summary, file.path(analysis_dir, "01_basicStats/feature_summary_filtered.tsv"))


## feautre groups

meth_features <- c("mean_meth", "coveredCpGs", "CpG_O", "CG_freq", "CpG_OE")
enzyme_features <- c("mspi_meth", "taq1_unmeth", "taq1_meth", "mspi_unmeth")
cpg_island_features <- c("Gardiner_Garden", "Takai_Jones")

kmers1 <- c("A", "C", "G", "T") # freq of T depends on others
kmers2 <- c("AA", "AC", "AG", "AT", "CA", "CC", "CpG_O", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"  )
kmers3 <- colnames(summary)[c(27:91)]

all_features <- list("meth" = meth_features, "enzyme_features" = enzyme_features,
                   "cpgs" = cpg_island_features, "kmers1" = kmers1, "kmers2" = kmers2, "kmers3" = kmers3)

## pairwise correlations within each group
summary=read.csv(file.path(analysis_dir, "01_basicStats/feature_summary_filtered.tsv"), sep="\t")

setDT(summary)

lapply(seq_along(all_features), function(x) draw_feature_correlations(all_features[[x]], 
                                                        names(all_features)[[x]], length(all_features[[x]])+1))


draw_feature_correlations(unique(unlist(all_features)), "ALL_F", length(unique(unlist(all_features))))
