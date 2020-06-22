source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))


## reading tissue - specific prediction data
stats_files_ts=system("ls /scratch/lab_bock/shared/projects/compEpi/results_analysis/02_predict_meth/02.3_tissue_vs_sample/screen/*/all_aucs.csv",intern=TRUE)
aucs_ts<-lapply(stats_files_ts, read.csv, sep="\t")

sp_list_ts<-unlist(lapply(stats_files_ts, function(x) strsplit(x, "/")[[1]][11] ))

aucs_ts<-lapply(aucs_ts, function(auc)
  {auc$test_class<-unlist(lapply(as.character(auc$species), function(x) all_class_sp[x][[1]])); return(auc)})

aucs_ts<-lapply(aucs_ts, function(auc)  
  {auc$train_species<-unique(auc[auc$ifRand=="noRandTrain", "species"]); return (auc)})

aucs_ts<-lapply(aucs_ts, function(auc)  
  {auc$train_class<-unlist(lapply(as.character(auc$train_species), function(x) 
    all_class_sp[x][[1]])); return (auc)})

all_aucs_ts<-rbindlist(aucs_ts)

all_aucs_ts<-all_aucs_ts %>%
  filter(train_tissue!="Tumour_early" & train_tissue!="Tumour_mid" & train_tissue!="Tumour_late")

all_aucs_ts<-as.data.table(all_aucs_ts)

#summary statistics
n_t<-stats_annot %>%
  group_by(species) %>%
  summarize(n_tisssues=length(unique(Tissue)), mean_tier=mean(mean_qual_tier), 
            class_color=first(color_class))

##plotting number of tissues in a species:
ggplot(as.data.frame(n_t), aes(x=as.factor(n_tisssues), fill=class_color))+
  geom_bar(stat="count", position="stack")+scale_fill_manual(values=class_colors)


## median tissue specificity - how different is the predicition score is in the tissue we train on
## compared to median prediction in other species 
## in case when there is two tissues, it is just auc_train/auc_test
##if it is <1, eans that the prdiction is better in another tissue

MTS<-function(auc){
  ##self_auc
  auc_1<-auc %>%
    filter(ifRand=="noRandTrain") %>%
    group_by(train_tissue) %>%
    summarize(self_auc=median(auc), species=first(train_species))
  #auc in the same species, but in different tissues
  auc_2<-auc %>%
    filter(ifRand=="noRandTest") %>%
    filter(train_species==species) %>%
    group_by(.dots=c("train_tissue")) %>%
    summarize(cross_ts_auc=median(auc))
  
  
  auc<- full_join(auc_1, auc_2, by="train_tissue") %>%
    mutate(cross_ts=self_auc/cross_ts_auc)
  return(auc)
}

non_one_ids<-seq_along(sp_list_ts)[sp_list_ts %in% n_t[n_t$n_tisssues>1,]$species]
mts_list<-lapply(aucs_ts[non_one_ids], MTS)

mts_full<-Reduce(rbind, mts_list)

mts_full$class<-sapply(mts_full$species, function(x) stats_annot[stats_annot$species==x,]$color_class[[1]])

mts_full$group<-factor(unlist(lapply(mts_full$class, function(x) class_short[x])), levels=class_short)

tissue_sum<-mts_full %>%
  group_by(train_tissue) %>%
  summarize(c=n()) %>%
  filter(c>10)


mts_full<-mts_full[mts_full$train_tissue %in% tissue_sum$train_tissue, ]


mts_full_short<-mts_full[mts_full$train_tissue %in% c("Heart", "Liver"), ]

p_tissue_sp<-ggplot(mts_full_short[mts_full_short$class=="Mammalia",], aes(x=train_tissue, y=species, fill=cross_ts))+geom_tile()+
  scale_fill_gradient2(low = "red", high = "red", mid = "blue", midpoint = 1, 
                       limit = c(0,2), space = "Lab") 
p_tissue_sp 


p_mts<-ggplot(mts_full_short, aes(x=group, y=cross_ts, fill=class))+geom_boxplot()+
  scale_fill_manual(values=class_colors)+
  facet_wrap(~train_tissue, ncol=1, scales="free")+
  labs(x="Phylogenetic group", y="Median tissue specificity score")+
  theme(legend.position = "None")
p_mts

pdf(paste0("/scratch/lab_bock/shared/projects/compEpi/results_analysis/02_predict_meth/02.3_tissue_vs_sample/summary/mts_h_l.pdf"),width=4,height=8)
print(p_mts)
#print(hm)
dev.off()

p_liver_mts<-ggplot(mts_full[mts_full$train_tissue=="Liver",], aes(x=class, y=cross_ts, fill=class))+geom_boxplot()+
  scale_fill_manual(values=class_colors)
p_liver_mts

pdf(paste0("/scratch/lab_bock/shared/projects/compEpi/results_analysis/02_predict_meth/02.3_tissue_vs_sample/summary/cross_ts_hm.pdf"),width=4,height=8)
print(p_tissue_sp)
#print(hm)
dev.off()