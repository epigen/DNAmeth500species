source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd = file.path(analysis_dir,"00.2_QC")
setwd(wd)

##first - plot the frequencies of each nucleotide
kmer1 <- t(read.csv(file.path(analysis_dir, "99.4_kmer_count/summary/kmer01.csv"), row.names=1))
kmer1 <- merge(kmer1, unique(stats_annot[, c("species", "color_class")]), by.x=0, by.y="species")
colnames(kmer1)[[1]]<-"species"
kmer1m<-melt(kmer1, id.var=c("species", "color_class"))

pdf("QC_plots/Nucleotide_composition.pdf", width = 4, height = 16)
ggplot(kmer1m, aes(x = color_class, y = value, fill = color_class)) +  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values = class_colors) + facet_wrap(~variable, ncol = 1)+rotate_labels()+
  labs(x = "phylogenetic group", y = "Nucleotide composition")+theme(legend.position = "None")
dev.off()

pdf("QC_plots/Nucleotide_composition_compare.pdf", width = 8, height = 16)
ggplot(kmer1m, aes(x = variable, y = value, fill = color_class)) +  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values = class_colors) + facet_wrap(~color_class, ncol = 2)+
  labs(x = "phylogenetic group", y = "Nucleotide composition")+theme(legend.position = "None")
dev.off()

#second - plot the relationship between C/G and A/T ratio
kmer1["C_to_G"] <- kmer1$C/kmer1$G
kmer1["A_to_T"] <- kmer1$A/kmer1$T
kmer1["T_to_A"] <- kmer1$T/kmer1$A

pdf("QC_plots/ratio_compare.pdf", width = 4, height = 4)
ggplot(kmer1, aes(x = A_to_T, y = C_to_G)) +  
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+xlim(0, 1.1)+ylim(0,1.1)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "A/T", y = "C/G")+theme(legend.position = "None")

ggplot(kmer1, aes(x = T_to_A, y = C_to_G)) +  
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+ylim(0, 1.2)+xlim(0.9, 2.1)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "T/A", y = "C/G")+theme(legend.position = "None")
dev.off()

pdf("QC_plots/ratio_compare_by_class.pdf", width = 3, height = 16)

ggplot(kmer1, aes(x = T_to_A, y = C_to_G)) +  
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+ylim(0, 1.2)+xlim(0.9, 2.1)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "T/A", y = "C/G")+theme(legend.position = "None")+facet_wrap(~color_class, ncol = 1)

dev.off()

pdf("QC_plots/ratio_boxplot.pdf", width = 4, height = 4)

ggplot(kmer1, aes(x = color_class, y = C_to_G, fill = color_class)) +  
  geom_boxplot(outlier.shape = 21)+rotate_labels()+
  scale_fill_manual(values = class_colors) + ylim(0,1)+
  labs( y = "C/G")+theme(legend.position = "None")

ggplot(kmer1, aes(x = color_class, y = T_to_A, fill = color_class)) +  
  geom_boxplot(outlier.shape = 21)+rotate_labels()+ylim(0.9, 2.1)+
  scale_fill_manual(values = class_colors) + 
  labs( y = "T/A")+theme(legend.position = "None")
dev.off()


###than we upload the kmer counts for the filtered sequences
kmer1_filtered <- t(read.csv(file.path(analysis_dir, "99.4_kmer_count_filtered/summary/kmer01.csv"), row.names=1))
kmer1_filtered <- merge(kmer1_filtered, unique(stats_annot[, c("species", "color_class")]), by.x=0, by.y="species")
colnames(kmer1_filtered)[[1]]<-"species"

kmer1m_filtered <- melt(kmer1_filtered, id.var=c("species", "color_class"))
kmer1m_filtered$state <- "filtered"
kmer1m$state <- "initial"

kmer_merged <- rbind(kmer1m_filtered, kmer1m)
kmer_merged$state <- factor(kmer_merged$state, levels=c("initial", "filtered"))

pdf("QC_plots/Nucleotide_composition_filtered_vs_unfiltered.pdf", width = 6, height = 4)
ggplot(kmer_merged, aes(x = variable, y = value, fill = state)) + geom_boxplot()+
  labs(x="", y="% N")
dev.off()


pdf("QC_plots/Nucleotide_composition_filtered_vs_unfiltered_by_class.pdf", width = 10, height = 16)
ggplot(kmer_merged, aes(x = variable, y = value, fill = state)) + geom_boxplot()+
  labs(x="", y="% N")+facet_wrap(~color_class, ncol = 2)
dev.off()

## 
kmer1_filtered["C_to_G_f"] <- kmer1_filtered$C/kmer1_filtered$G
kmer1_filtered["T_to_A_f"] <- kmer1_filtered$T/kmer1_filtered$A


pdf("QC_plots/ratio_boxplot_filtered.pdf", width = 4, height = 4)
ggplot(kmer1_filtered, aes(x = color_class, y = T_to_A_f, fill = color_class)) +  
  geom_boxplot(outlier.shape = 21)+rotate_labels()+
  scale_fill_manual(values = class_colors) + 
  labs( y = "T/A")+theme(legend.position = "None")

ggplot(kmer1_filtered, aes(x = color_class, y = C_to_G_f, fill = color_class)) +  
  geom_boxplot(outlier.shape = 21)+rotate_labels()+
  scale_fill_manual(values = class_colors) + 
  labs( y = "C/G")+theme(legend.position = "None")
dev.off()

kmer1_and_f <- merge(kmer1, kmer1_filtered[, c("species", "T_to_A_f", "C_to_G_f")], by="species")

pdf("QC_plots/ratio_compare_full.pdf", width = 4, height = 4)

ggplot(kmer1, aes(x = T_to_A, y = C_to_G)) +  
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+ylim(0, 1.2)+xlim(0.9, 2.1)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "T/A", y = "C/G")+theme(legend.position = "None")


ggplot(kmer1_filtered, aes(x = T_to_A_f, y = C_to_G_f)) +  
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+ylim(0.8, 1.2)+xlim(0.8, 1.2)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "T/A_f", y = "C/G_f")+theme(legend.position = "None")

ggplot(kmer1_and_f, aes(x = T_to_A, y = T_to_A_f)) +  ylim(0.8, 2.1)+xlim(0.8, 2.1)+
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "T/A", y = "T/A_f")+theme(legend.position = "None")
ggplot(kmer1_and_f, aes(x = C_to_G, y = C_to_G_f)) +ylim(0.2, 1.2)+xlim(0.2, 1.2)+
  geom_point(shape = 21, aes(fill = color_class, alpha = 0.7, color = color_class))+
  geom_abline(linetype="dashed", alpha = 0.7)+
  scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors)+
  labs(x = "C/G", y = "C/G_f")+theme(legend.position = "None")
dev.off()

## filtering summary
stats_files = system("ls filtered_sequences/*/filter_stats.csv ",intern=TRUE)
sp_list <- unlist(lapply(stats_files, function(x) strsplit(x, "/")[[1]][2] ))
filtering_stats=data.table()

for (file in stats_files){
  f=fread(file)
  f$species<-strsplit(file, "/")[[1]][2]
  f$color_class<-stats_annot[stats_annot$species==f$species, ]$color_class[[1]]
  filtering_stats=rbindlist(list(filtering_stats,f))
}


filtering_stats$starts_with_cgg <- filtering_stats$only_cgg/filtering_stats$initial
filtering_stats$maps_uc_and_starts_cgg <- filtering_stats$cgg_and_match_uc/filtering_stats$initial

pdf("QC_plots/filtering_persentage.pdf", width = 6, height = 6)
ggplot(filtering_stats, aes(x=starts_with_cgg, fill = color_class))+geom_histogram()+
  scale_fill_manual(values = class_colors)

ggplot(filtering_stats, aes(x=maps_uc_and_starts_cgg, fill = color_class))+geom_histogram()+
  scale_fill_manual(values = class_colors)
dev.off()

fsm <- melt(filtering_stats[, -c("starts_with_cgg", "maps_uc_and_starts_cgg")], 
            id.vars = c("color_class", "species"))

pdf("QC_plots/dedRef_count_density.pdf", width = 10, height = 16)
ggplot(fsm, aes(x=value))+geom_density(aes(fill = variable, alpha = 0.1, color = variable))+
  facet_wrap(~color_class, scales="free_y", ncol=2)+rotate_labels()+labs(x="N of dedRefs")
dev.off()
