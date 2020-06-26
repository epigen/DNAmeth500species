source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

mywd = file.path(analysis_dir, "/02_predict_meth/02.4_verify_inverted/")
setwd(mywd)

dir.create(file.path(mywd, "summary"))

## showing tissue-specific patterns:
tissue_inv <- fread(file.path(mywd,"tissue_check","WHH", "all_data.csv"))

tissue_norm <- fread(file.path(mywd,"tissue_check","AMP", "all_data.csv"))
tissues_common <- intersect(unique(tissue_norm$train_tissue), unique(tissue_inv$train_tissue))
tissue_inv <- tissue_inv[train_tissue %in% tissues_common]
tissue_norm <- tissue_norm[train_tissue %in% tissues_common]

tissue_inv[, c("fdr", "tpr"):=NULL]
tissue_inv <- unique(tissue_inv)
tissue_inv$species <- "WHH"

tissue_norm[, c("fdr", "tpr"):=NULL]
tissue_norm <- unique(tissue_norm)
tissue_norm$species <- "AMP"

tissue <- rbind(tissue_inv, tissue_norm)
tissue <- tissue[ifRand=="noRandTest"]

ggplot(tissue, aes(x = auc, fill = species)) + geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c("WHH" = "#ca0020", "AMP" = "#0571b0")) + 
  geom_point(data = tissue[species == type], 
             aes(x = auc, y = 0, fill = species), size = 3, alpha = 0.5, shape = 21) + 
  scale_color_manual(values = c("WHH" = "#ca0020", "AMP" = "#0571b0")) +
  facet_wrap(~train_tissue, nrow = 1) + 
  theme(text = element_text(size = 15))
ggsave("summary/tissue_check.pdf", width = 12, height = 3)

rm(tissue)
rm(tissue_inv)
rm(tissue_norm)
   
### individuals check:
samples_inv <- fread(file.path(mywd,"sample_check","WHH", "all_data.csv"))
samples_inv <- unique(samples_inv[, c("fdr", "tpr"):=NULL])
samples_inv$species <- "WHH"

samples_norm <- fread(file.path(mywd,"sample_check","AMP", "all_data.csv"))
samples_norm <- unique(samples_norm[, c("fdr", "tpr"):=NULL])
samples_norm$species <- "AMP"

samples <- rbind(samples_inv, samples_norm)
samples<-samples[ifRand=="noRandTest"]
samples$id <- paste(samples$species, samples$sample_N, sep = "_")
ggplot(samples, aes(x = auc,  fill = id)) + geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c("WHH_2"  = "#a50026", "WHH_1" ="#d73027", 
                               "AMP_1" = "#4575b4", "AMP_2" = "#74add1")) + 
  geom_point(data = samples[species == type], 
             aes(x = auc, y = 0, fill = id), size = 5,  shape = 21) +
  theme(text = element_text(size = 15))
ggsave("summary/sample_check.pdf", width = 6, height = 3)


###bootstraping experiment - here I just saved the AUCs
boot_inv_files  <- system("ls bootstrap/WHH/*_stats.tsv",intern=TRUE)
boot_inv_list <- sapply(boot_inv_files, read.csv, sep = "\t", simplify = F)
boot_inv <- rbindlist(boot_inv_list)
boot_inv$species <- "WHH"

boot_norm_files  <- system("ls bootstrap/AMP/*_stats.tsv",intern=TRUE)
boot_norm_list <- sapply(boot_norm_files, read.csv, sep = "\t", simplify = F)
boot_norm <- rbindlist(boot_norm_list)
boot_norm$species <- "AMP"

##############################

boot_df <- rbind(boot_inv, boot_norm)

ggplot(boot_df, aes(x = species, y = AUC, fill = species)) + geom_boxplot() + 
          geom_jitter() + scale_fill_manual(values = c("WHH" = "#ca0020", "AMP" = "#0571b0")) +
        theme(text = element_text(size = 15))
ggsave("summary/bootstrap.pdf", width = 5, height = 5)
