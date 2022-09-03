source(file.path("../../src/00.0_init.R"))

wd <- file.path(analysis_dir, "validation", "03_WGBS", "03.6_WGBS_test_in_RRBS")

setwd(wd)



summary_files = list.files(pattern = "auc_test", recursive = TRUE)

summary_auc_list <- sapply(summary_files, fread, simplify = F)

summary_auc <- rbindlist(summary_auc_list)



summary_auc <- summary_auc[type %in% sp_df$species]

summary_auc[species_train == "Danio_rerio_GSE134055",species_train:="Danio_rerio",]

summary_auc$species_train <- gsub("_", " ", summary_auc$species_train)

summary_auc <- left_join(summary_auc, sp_df, by=c("type" = "species"))

summary_auc_mean <- summary_auc %>% group_by(species_train, group) %>% summarize(mean_auc = mean(auc), sd_auc = sd(auc))

WGBS_species_order <- c('Branchiostoma lanceolatum','Danio rerio','Xenopus laevis', 'Chelydra serpentina','Gallus gallus', 'Phascolarctos cinereus', 'Mus musculus','Bos taurus')

head(summary_auc_mean)

summary_auc_mean$species_train <- factor(summary_auc_mean$species_train, levels = WGBS_species_order)

p <- ggplot(summary_auc_mean, aes(x = group, y = species_train, fill = mean_auc)) + 
geom_tile() +  
scale_fill_gradient2(limits = c(0,1), low = "#2b83ba",mid ="#ffffbf",  high ="#d7191c", midpoint = 0.5 ) + 
geom_text(aes(label = paste0(round(mean_auc,2), "\n± ", round(sd_auc,2))), size = 2) + xlab("test") + ylab("train")
p

ggsave("RRBS_WGBS_summary_heatmap.pdf", p, width = 6, height = 4)



