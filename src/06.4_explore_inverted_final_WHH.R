##This is not a script to run, some analysis is involeved at some steps
library(msa)
library(tibble)
library(ggpubr)
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

SP <- "WHH"
wd = file.path(analysis_dir,"05_predict_meth/05.1_within_species")
setwd(wd)

subdir <- paste0(analysis_dir, "/06_inv/WHH_final/")
dir.create(subdir)

### uploading inverted species:
invfile <- paste0(Sys.getenv("CODEBASE"), "DNAmeth500species/meta/inverted_species.txt")
if(!file.exists(invfile)) print("inverted list does not exist, run 06.1 first")
inv <- read.table(invfile)$V1

sp_full_order <- read.table(paste0(Sys.getenv("CODEBASE"), "DNAmeth500species/meta/species_list_ordered_2021_2.txt"), header = 1)$V1

sp_df$inv <- sp_df$species %in% inv

##reading AUC summary
aucs <- fread(file.path(analysis_dir, "05_predict_meth","/05.2_test_on_other_species/summary/all_aucs_full.csv"))
self_aucs <- read.csv(file.path(analysis_dir,"05_predict_meth/05.1_within_species/summary/all_aucs.csv"), row.names = 1)
self_aucs<-unique(self_aucs[, c("species", "k", "AUC", "numSequences", "color_class", "group")])
row.names(self_aucs)<-self_aucs$species

##loading the feature weights
all_kmer3 <- read.table("../05.3_feature_weight_analysis/all_kmerWeights_3.csv", row.names = 1, header = 1)


##draw the heatmap of all the feature weights for the Actionpteri
sub_auc <- self_aucs[sp_df[sp_df$color_class == "Actinopteri",]$species,]

ha = columnAnnotation(inverted = sp_df[color_class == "Actinopteri",]$inv,  selfAUC = sub_auc$AUC,
 col=list(inverted=c("TRUE"="red","FALSE"="grey")))

pdf(paste0(subdir, "Actionpteri_f_weights.pdf"), width = 15, height = 10)
print(Heatmap(all_kmer3[sp_df[color_class=="Actinopteri",]$species], cluster_columns = F,
              top_annotation = ha, name="feature weights", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8)))
dev.off()

##analysing the feature weights
#identifying, which of the species in the same class showes the most "invertiness"
most_opp <- as.character(aucs[aucs$auc == min(aucs[train_species==SP & color_class_test == "Actinopteri" &  ifRand == "noRandTest", ]$auc),]$type)
sub <- all_kmer3[, c(SP, most_opp)]
sub$mean<-rowMeans(all_kmer3[, sp_df[sp_df$color_class=="Actinopteri" & !sp_df$inv, ]$species])

sd_vec<-apply(all_kmer3[, sp_df[sp_df$color_class=="Actinopteri" & !sp_df$inv, ]$species], 1, sd)
sub$kmer <- row.names(sub)
sub_m <- melt(sub, id.vars=c("kmer"))
sub_m$sd <- 0

sub_m[sub_m$variable=="mean", ]$sd <- sd_vec
sub_m$variable <- factor(sub_m$variable, levels=c( "WHH","PLF", "mean"))

##drawing barplot
pdf(paste0(subdir, "/kmer_weights_Actinopteri_WHH_PLF.pdf"),
    width = 20, heigh = 10)
ggplot(sub_m, aes(x = kmer, y = value )) + 
  geom_bar(stat = "identity", aes(fill = variable), alpha = 0.7, position = "identity") + rotate_labels()+
  scale_fill_manual(values = c( "mean" = "grey","WHH" = "red", "PLF" = "lightblue")) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.6, color = "grey", position=position_dodge())
dev.off()

sub$sd <- sd_vec
sub$top <- sub$mean + sub$sd

kmer_pos<-sub %>%
  filter(sign(WHH) == 1 & sign(mean) == -1 & sign(PLF) == -1) %>%
  filter(WHH > top) %>%
  .$kmer

sub$low <- sub$mean - sub$sd

kmer_neg <- sub %>%
  filter(sign(WHH) == -1 & sign(mean)== 1 & sign(PLF) == 1) %>%
  filter(WHH < low) %>%
  .$kmer

print(kmer_pos)
print(kmer_neg)

#sub_m_outliers <- sub_m[sub_m$kmer %in% union(kmer_pos, kmer_neg), ]
sub_m_outliers <- sub_m[sub_m$kmer %in% kmer_pos, ]
kmer_levels <- sub_m_outliers[sub_m_outliers$variable == "WHH", ][order(sub_m_outliers[sub_m_outliers$variable == "WHH", ]$value),"kmer" ]

sub_m_outliers$kmer <- factor(sub_m_outliers$kmer, levels = kmer_levels)
setDT(sub_m_outliers)

pdf(paste0(subdir, "/kmer_weights_Actinopteri_WHH_PLF_different.pdf"),
    width = 3, heigh = 3)
ggplot(sub_m_outliers[variable!="PLF" ], aes(x = kmer, y = value )) + 
  geom_bar(stat = "identity", aes(fill = variable), alpha = 0.7, position = "identity") + rotate_labels()+
  scale_fill_manual(values = c( "mean" = "grey","WHH" = "red", "PLF" = "lightblue")) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.6, color = "grey", position=position_dodge())
dev.off()
##selecting three top kmers
kmer_pos <- kmer_pos[2:4]
#first step is to analyze if there are repeats of the direct feature weights, that cause the shift of the svm model
rep <- paste0(kmer_pos, kmer_pos)

stats_high <- read.csv(paste0(analysis_dir,  "/06_inv/06.3_kmer_freq/repeat_frequencies_6/kmer_freq_high_6.csv"), row.names = 1)

stats_low <- read.csv(paste0(analysis_dir,
 "/06_inv/06.3_kmer_freq/repeat_frequencies_6/kmer_freq_low_6.csv"), row.names = 1)
##as a color I want to use the gradient of "invertiness of the species" in relation to the WHH
auc_d <- aucs %>%
  filter( ifRand == "noRand" | ifRand == "noRandTest" ) %>%
  filter( train_species == "WHH" & color_class_test == "Actinopteri" ) %>%
  mutate ( train_AUC = as.numeric(aucs[aucs$ifRand == "noRand" & aucs$train_species== "WHH", "auc"])) %>%
  mutate( auc_d = train_AUC - auc) %>%
  select(type, auc_d, auc)

make_df <-function(kmer){
  comp<-cbind(t(stats_high[stats_high$kmer==kmer, c(-1)]), t(stats_low[stats_low$kmer==kmer, c(-1)]))
  colnames(comp)<-c("high", "low")
  comp <- as.data.frame(comp) %>% rownames_to_column()
  comp<-full_join(comp, auc_d, by = c("rowname" = "type"))
  comp<- as.data.frame(comp)[complete.cases(comp),]
  comp$rep_kmer <- kmer
  return(comp)
}


####analysing the repeats of 3

rep3 <- paste0(kmer_pos, kmer_pos, kmer_pos)

stats_high <- read.csv(paste0(analysis_dir,                             "/06_inv/06.3_kmer_freq/repeat_frequencies_6/kmer_freq_high_9.csv"), row.names = 1)

stats_low <- read.csv(paste0(analysis_dir, "/06_inv/06.3_kmer_freq/repeat_frequencies_6/kmer_freq_low_9.csv"), row.names = 1)


comp_list3 <- lapply(rep3, make_df)

comp_all3<-rbindlist(comp_list3)
comp_all3$inv <- comp_all3$rowname %in% inv


##correlation coefficients:
comp_corr <- comp_all3 %>% 
  group_by(.dots = c("rep_kmer", "inv")) %>%
  summarise(corr = round(cor.test(high, low)$estimate,3), cor.pval = cor.test(high, low)$p.value, 
            x = 0.8*max(low), y = min(high))

comp_corr$corr_label = paste0("R= ", comp_corr$corr)


pdf(paste0(subdir,"repeat_frequencies_upgrade.pdf"), width = 20, height = 5)
ggplot(comp_all3, aes(x = low, y = high, color = auc)) + geom_point() + facet_wrap(~rep_kmer, scales = "free", nrow = 1)+ 
  geom_text_repel(data = comp_all3[comp_all3$inv, ], aes(x = low, y = high, label = rowname), color = "black") + 
  scale_color_gradient2(low = "#2b83ba", mid ="#ffffbf", high = "#d7191c", midpoint = 0.5) + 
  geom_smooth(method='lm', aes(group = inv), color = "black", size = 0.1, linetype = "dashed", se =F) + 
  labs(x = "low methylated dedRefs", y = "highly methylated dedRefs", color = "ROC-AUC") + 
  geom_text(data = comp_corr[!comp_corr$inv, ], aes(x = x, y = y, label = corr_label), color = "black")
dev.off()




pdf(paste0(subdir,"correlations_for_repeats.pdf"), width = 4, height = 4)
ggplot(comp_corr, aes(x = inv, y = rep_kmer)) + geom_point(aes(size = corr, color = log(cor.pval))) 
dev.off()


comp_all3$freq_d<-comp_all3$high - comp_all3$low

pdf(paste0(subdir, "auc_freq_corr.pdf"), width = 18, height = 5)
ggplot(comp_all3, aes(x = freq_d, y = auc, color = inv)) + geom_point() + 
  stat_cor(method = "pearson", label.x.npc = 0.5, label.y.npc = 0.5) + 
  facet_wrap(~rep_kmer, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TRUE"="red","FALSE"="darkgrey")) +
  geom_text_repel(data = comp_all3[comp_all3$inv, ], aes(x = freq_d, y = auc, label = rowname), color = "black") + 
    geom_vline(xintercept = 0.0, linetype = "dashed", size = 0.5)+ 
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.5)
dev.off()

