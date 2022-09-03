source(file.path("../../src/00.0_init.R"))

wd <- file.path(analysis_dir, "validation", "03_WGBS", "03.5_RRBS_test_in_WGBS")

setwd(wd)

summary_files = list.files(pattern = "auc_test", recursive = TRUE)



summary_auc_list <- sapply(summary_files, fread, simplify = F)

summary_auc <- rbindlist(summary_auc_list)



sp_list <- fread(file.path(meta_dir, "species_list.txt"), header = F)$V1

summary_auc[type == "Danio_rerio_GSE134055", type:="Danio_rerio",]

summary_auc$type <- gsub("_", " ", summary_auc$type)
summary_auc <- inner_join(summary_auc, sp_df, by=c("species_train" = "species"))

summary_auc[is.na(summary_auc$group)]

summary_auc_mean <- summary_auc %>% group_by(type, group) %>% summarize(mean_auc = mean(auc), sd_auc = sd(auc))


WGBS_species_order <- c('Branchiostoma lanceolatum','Danio rerio',
                        'Xenopus laevis', 'Chelydra serpentina',
                        'Gallus gallus', 'Phascolarctos cinereus', 'Mus musculus','Bos taurus')


summary_auc_mean$type <- factor(summary_auc_mean$type, levels = WGBS_species_order)

setDT(summary_auc_mean)

head(summary_auc_mean[is.na(group)])

p <- ggplot(summary_auc_mean, aes(x = group, y = type, fill = mean_auc)) + 
geom_tile() +  
scale_fill_gradient2(limits = c(0,1), low = "#2b83ba",mid ="#ffffbf",  high ="#d7191c", midpoint = 0.5 ) + 
geom_text(aes(label = paste0(round(mean_auc,2), "\n± ", round(sd_auc,2))), size = 2) + xlab("train") + ylab("test")
p

ggsave("WGBS_RRBS_summary_heatmap.pdf", p, width = 6, height = 4)

