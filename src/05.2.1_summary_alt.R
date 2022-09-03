source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"05_predict_meth/05.1_within_species")
setwd(wd)

##collecting AUC results
if(!file.exists("summary/aucs_alternative.csv")){
  stats_files <- system("ls screen_alt/*/roc_res.csv",intern=TRUE)
  stats_files_list <- sapply(stats_files, function(x) 
    unique(read.csv(x, row.names = 1)[,c("auc", "f1","ifRand","type", "N_conv_samples")]), simplify = F)
  aucs_alt <- rbindlist(stats_files_list)
  write.table(aucs_alt,"summary/aucs_alternative.csv", quote = F, row.names = F)
}

aucs_alt <- read.table("summary/aucs_alternative.csv", header = 1)
head(aucs_alt)
colnames(aucs_alt)[[1]] <- "AUC_alt"
colnames(aucs_alt)[[2]] <- "f1_alt"

##uploading initial resuls
aucs <- read.csv("summary/all_aucs.csv", row.names = 1)
aucs <- right_join(aucs[!duplicated(aucs$species),], aucs_alt[aucs_alt$ifRand == "noRand",], by = c("species" = "type"))

##calculating difference
aucs$delta_AUC = aucs$AUC - aucs$AUC_alt

                             
ggplot(aucs,aes(x = N_conv_samples, fill=color_class)) + 
  geom_histogram() + scale_fill_manual(values = class_colors)
ggsave("summary/N_conv_samples.pdf", width = 6, height = 4)

pdf("summary/changes_in_AUC.pdf", width = 6, height = 4)
ggplot(aucs, aes(x = AUC, y = AUC_alt)) + 
  geom_point(shape = 21,aes(fill = color_class)) + 
  scale_fill_manual(values = class_colors) + 
  geom_text_repel(data = aucs[abs(aucs$delta_AUC) >= 0.015, ], aes(x = AUC, y = AUC_alt, label = species)) + 
    xlim(c(0,1)) + ylim(c(0,1)) + 
  ggtitle("changes in AUC")

ggplot(aucs, aes(x = AUC, y = AUC_alt)) + 
  geom_point(shape = 21,aes(fill = N_conv_samples)) +
  scale_fill_gradient2(midpoint = 20, low = "#e5f5f9", mid = "#99d8c9", 
                       high = "#2ca25f") + ggtitle("changes in AUC")+ 
  geom_text_repel(data = aucs[abs(aucs$delta_AUC) >= 0.015, ], aes(x = AUC, y = AUC_alt, label = species))

ggplot(aucs[aucs$N_conv_samples >= 2, ], aes(x = AUC, y = AUC_alt)) + 
  geom_point(shape = 21,aes(fill = color_class)) + 
  scale_fill_manual(values = class_colors)+ 
  ggtitle("changes in AUC (N samples > 1)")

ggplot(aucs, aes(x = delta_AUC, y = N_conv_samples)) + 
  geom_point(shape = 21,aes(fill = color_class)) + 
  scale_fill_manual(values = class_colors)+ 
  ggtitle("delta AUC")
dev.off()


pdf("summary/changes_in_AUC_full.pdf", width = 6, height = 4)

ggplot(aucs, aes(x = AUC, y = AUC_alt)) + 
  geom_point(shape = 21,size = 3, alpha = 0.5, aes(fill = color_class)) + 
  scale_fill_manual(values = class_colors) +coord_equal() + 
  ggtitle("changes in AUC")
  
  ggplot(aucs[aucs$N_conv_samples >= 2, ], aes(x = AUC, y = AUC_alt)) + 
    geom_point(shape = 21,size = 3, aes(fill = color_class)) + 
    scale_fill_manual(values = class_colors) +coord_equal() + 
    ggtitle("changes in AUC(N samples > 1)")
dev.off()
