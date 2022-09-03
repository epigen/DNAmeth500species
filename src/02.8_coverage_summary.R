## ---------------------------
##
## Coverage analysis - summary
##
## Combining and visualizing the results of the 01.9 parallelized script
## For each species a coverage_stats.tsv file is generated
## Authors: Daria Romanovskaia
##
## Date Created: 2021-11-04
##
##
## ---------------------------

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd = file.path(analysis_dir, "01_basicStats", "01.9_coverage_analysis")
setwd(wd)

getwd()

library(patchwork)

dir.create("summary")

## collecting files with coverage information. We are mostly interested in columns N_rep and N_amp, that indicate the number of repeats

if(!file.exists( "summary/coverage_stats_combined.tsv")){
files = list.files(pattern = "coverage_stats.tsv", recursive = T)

coverage_stats_list <- sapply(files, fread, simplify = FALSE)

coverage_stats <- rbindlist(coverage_stats_list)

coverage_stats[, species:=strsplit(sample, "_")[[1]][1], by = row.names(coverage_stats)]

coverage_stats[, tissue:=strsplit(sample, "_")[[1]][3], by = row.names(coverage_stats)]

coverage_stats <- setDT(left_join(coverage_stats, sp_df))

coverage_stats[, individual := strsplit(sample, "_")[[1]][2], by = row.names(coverage_stats)]
    
coverage_stats[,N_rep_norm:=N_rep/N_frags,]
coverage_stats[,N_amp_norm:=N_amp/N_frags,]
coverage_stats[,N_ind_norm:=N_ind/N_frags,]
    
my_wt(coverage_stats, "summary/coverage_stats_combined.tsv")
}else coverage_stats <- fread( "summary/coverage_stats_combined.tsv")

#### filtering out the samples,excluded from 
stats_annot=stats_annot[!grepl("_uc$",Sample_Name)&species_check!="fail"&!grepl("Tumour|Cellline",Tissue)] 

coverage_stats <- coverage_stats[sample %in% stats_annot$Sample_Name,]

max(coverage_stats$individual)

coverage_stats[individual==10]

## normalizing the N_rep by the number of fragments

coverage_stats_mean <- coverage_stats %>% 
    group_by(species) %>% 
    summarise(n_tissue = length(unique(tissue)), n_samples = n(), n_indiv = max(individual), N_rep_norm_mean = mean(N_rep_norm),
              N_ind_norm_mean=mean(N_ind_norm), N_amp_norm_mean = mean(N_amp_norm), color_class = first(color_class), group = first(group))

setDT(coverage_stats_mean)

head(coverage_stats_mean)


wd = file.path(analysis_dir, "02_vizStats", "02.8_coverage_analysis")
dir.create(wd)
setwd(wd)

coverage_stats_m <- melt(coverage_stats_mean,
     measure.vars = c("N_rep_norm_mean", "N_amp_norm_mean",  "N_ind_norm_mean"))


p1 <- ggplot(coverage_stats_m, aes(x = as.factor(n_tissue), y = value, color = color_class)) + 
    geom_jitter() + facet_grid(~variable) + 
    scale_color_manual(values = class_colors) + theme(legend.position = "None")

p2 <- ggplot(coverage_stats_m, aes(x = as.factor(n_samples), y = value, color = color_class)) + 
    geom_jitter() + facet_grid(~variable) + 
    scale_color_manual(values = class_colors) + theme(legend.position = "None") + geom_vline(xintercept = 3.5, linetype = "dashed")

p3 <- ggplot(coverage_stats_m, aes(x = as.factor(n_indiv), y = value, color = color_class)) + 
    geom_jitter() + facet_grid(~variable) + 
    scale_color_manual(values = class_colors) + theme(legend.position = "None")+ geom_vline(xintercept = 1.5, linetype = "dashed")


g <- p1/p2/p3
g

ggsave("stats_biases_mean.pdf",g, width = 20, height = 10)

coverage_stats_filtered <- coverage_stats_mean[n_indiv > 1 & n_samples >3]

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


coverage_stats_filtered$group <- factor(coverage_stats_filtered$group, levels = class_short)

ggplot(coverage_stats_filtered, aes(x = group, y = N_rep_norm_mean,fill = color_class)) + 
    geom_boxplot(outlier.shape = 21) + 
    scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors) + 
    stat_summary(fun.data = give.n,fun.args = c(y=0.02), geom = "text",size=4)
ggsave("stats_N_rep_norm_mean.pdf", width = 8, height = 4)


ggplot(unique(coverage_stats_filtered), aes(x = group, y = N_ind_norm_mean,fill = color_class)) + 
    geom_boxplot(outlier.shape = 21) + 
    scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors) + 
    stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4)
ggsave("stats_N_ind_norm_mean.pdf", width = 8, height = 4)


ggplot(unique(coverage_stats_filtered), aes(x = group, y = N_amp_norm_mean,fill = color_class)) + 
    geom_boxplot(outlier.shape = 21) + 
    scale_fill_manual(values = class_colors) + scale_color_manual(values = class_colors) + 
    stat_summary(fun.data = give.n,fun.args = c(y=0.02), geom = "text",size=4)
ggsave("stats_N_amp_norm_mean.pdf", width = 8, height = 4)


ggplot(coverage_stats_filtered, aes(x = N_amp_norm_mean, y = N_rep_norm_mean, color = color_class)) + 
    geom_point() + scale_color_manual(values = class_colors) + 
    geom_label_repel(data = coverage_stats_filtered[N_amp_norm_mean>0.01 | N_rep_norm_mean > 0.015, ], aes(x = N_amp_norm_mean, y = N_rep_norm_mean, label = species))

ggsave("stats_N_rep_vs_amp_norm_mean.pdf", width = 6, height = 4)


coverage_stats <- left_join(coverage_stats, stats_annot[, c("Sample_Name", "Enrichment cycles")], by = c("sample" = "Sample_Name"))

ggplot(coverage_stats, aes(x = `Enrichment cycles`, fill = color_class)) + 
                geom_histogram() + 
                scale_fill_manual(values = class_colors)

setDT(coverage_stats)

coverage_stats <- coverage_stats[species %in%coverage_stats_filtered$species ]


ggplot(coverage_stats, aes(y = N_amp_norm, x = `Enrichment cycles`, color = color_class)) + 
    geom_jitter() + scale_color_manual(values = class_colors) + 
    geom_label_repel(data = coverage_stats[N_amp_norm>0.02, ], aes(y = N_amp_norm, x = `Enrichment cycles`, label = sample)) + 
    ggtitle(paste0("corr = ", round(cor(coverage_stats$N_amp_norm, coverage_stats$`Enrichment cycles`),3)))

ggsave("stats_PCR_vs_amp_norm.pdf", width = 6, height = 4)

coverage_stats$color_class <- factor(coverage_stats$color_class, levels = names(class_colors))

ggplot(coverage_stats, aes( x = as.factor(`Enrichment cycles`), fill = color_class)) + 
    geom_bar() + scale_fill_manual(values = class_colors) + xlab("Enrichment cycles") + ylab("")
ggsave("PCR_count.pdf", width = 6, height = 4)


ggplot(coverage_stats, aes(y = N_amp_norm, x = as.factor(`Enrichment cycles`))) + 
    geom_boxplot(aes(, color = color_class), outlier.size = 0.5, position = position_dodge(preserve = "single")) +
        scale_color_manual(values = class_colors) + 
        stat_summary(fun.data = give.n,fun.args = c(y=-0.005), geom = "text",size=3)+
    geom_label_repel(data = coverage_stats[N_amp_norm>0.04, ], 
                   aes(y = N_amp_norm, x = as.factor(`Enrichment cycles`),color = color_class, label = sample))
ggsave("PCR_n_amp.pdf", width = 9, height = 4)
