source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/03.0_predict_meth_init.R"))

wd = file.path(analysis_dir, "03_prediction_lm/03.1_features")
dir.create(wd, recursive = TRUE)
setwd(wd)
#for the further analysis we are making a table with the key data  from the sequences:
#important statistics
summary_conv<-stats_annot %>% 
  filter(conversion_type=="converted") %>%
  group_by(species) %>% summarize(mean_meth=mean(CpG_meth),coveredCpGs=mean(coveredCpGs), color_class=first(color_class), 
            tier=mean(mean_qual_tier))

summary_unconv<-stats_annot %>% 
  filter(conversion_type=="unconverted") %>%
  group_by(species) %>% summarize(mean_meth_unconv=mean(CpG_meth),coveredCpGs_unconv=mean(coveredCpGs))

summary<-setDT(full_join(summary_conv, summary_unconv, by="species"))


#plotting converted vs unconverted
pdf("coveredCpGs_conv_vs_unconv.pdf")
draw_correlation(summary$coveredCpGs, summary$coveredCpGs_unconv, "covered_CpGs_converted", "covered_CpGs_unconverted")
dev.off()

pdf("mean_meth_conv_vs_unconv.pdf")
draw_correlation(summary$mean_meth, summary$mean_meth_unconv, "mean_meth_converted", "mean_meth_unconverted")
dev.off()

pdf("mean_meth_quality.pdf")
draw_correlation(summary$mean_meth, summary$tier, "mean_meth", "quality")
dev.off()


rm(summary_conv)                  
rm(summary_unconv)


kmer1<-t(read.csv(file.path(analysis_dir, "01_basicStats","01.5_kmercount", "summary","kmer01.csv"), row.names=1))

kmer2<-t(read.csv(file.path(analysis_dir, "01_basicStats","01.5_kmercount", "summary","kmer02.csv"), row.names=1))

kmer3<-t(read.csv(file.path(analysis_dir, "01_basicStats","01.5_kmercount", "summary","kmer03.csv"), row.names=1))

frequencies<-merge(kmer1, kmer2,  by=0)
row.names(frequencies)<-frequencies$Row.names

frequencies<-merge(frequencies[, -1], kmer3, by=0)

rm(kmer1)
rm(kmer2)
rm(kmer3)

colnames(frequencies)[1]<-"species"
frequencies$species <- as.character(frequencies$species)
                                    
summary<-left_join(summary, frequencies)

#calculating the sequence - derived features
colnames(summary)[which(colnames(summary)=="CG")] <- "CpG_O"

summary$CG_freq<-summary$G+summary$C
summary$CpG_OE<-summary$CpG_O/(summary$C*summary$G)


##uploading CpG - island frequencies:
if(!file.exists("CpG_island_per_ref.tsv")){
   stats_files=system(paste0("ls ", analysis_dir, "/01_basicStats/01.8_CpG_islands_per_ref/*.csv"), intern=TRUE)
file_list <- lapply(stats_files, read.csv)
CpG_cr <- rbindlist(file_list)

colnames(CpG_cr)<-c("species", "Gardiner_Garden", "Takai_Jones")
my_wt(CpG_cr, "CpG_island_per_ref.tsv") 
}else{
    CpG_cr <- fread("CpG_island_per_ref.tsv")
}

                                    
summary <- left_join(summary, CpG_cr)

##adding the self AUC
auc_stats<-read.csv(file.path(analysis_dir, "05_predict_meth/05.1_within_species/summary/all_aucs.csv"), row.names = 1)

summary<-left_join(summary, unique(auc_stats[,c("species", "k", "AUC")]))


pdf("feature_summary.pdf", width=4, height=4*(length(summary)-2))
ggplot(melt(summary, id.vars = c("species", "color_class")), aes(x=value, fill=color_class)) + 
  geom_histogram() + facet_wrap(~variable, scales = "free", ncol = 1) + 
  scale_fill_manual(values=class_colors) + theme(legend.position = "None")
dev.off()

my_wt(summary, file.path("feature_summary_filtered.tsv"))


## feautre groups

meth_features <- c("mean_meth", "coveredCpGs", "CpG_O", "CG_freq", "CpG_OE")
enzyme_features <- c("mspi_meth", "taq1_unmeth", "taq1_meth", "mspi_unmeth")
cpg_island_features <- c("Gardiner_Garden", "Takai_Jones")

kmers1 <- c("A", "C", "G", "T") # freq of T depends on others
kmers2 <- c("AA", "AC", "AG", "AT", "CA", "CC", "CpG_O", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"  )
kmers3 <- colnames(summary)[c(28:91)]

all_features <- list("meth" = meth_features,
                   "cpgs" = cpg_island_features, "kmers1" = kmers1, "kmers2" = kmers2, "kmers3" = kmers3)

## pairwise correlations within each group
summary=fread(file.path( "feature_summary_filtered.tsv"), sep="\t")

lapply(seq_along(all_features), function(x) draw_feature_correlations(all_features[[x]], 
                                                        names(all_features)[[x]], length(all_features[[x]])+1))


draw_feature_correlations(unique(unlist(all_features)), "ALL_F", length(unique(unlist(all_features))))
