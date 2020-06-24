source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"02_predict_meth/02.1_within_species")
setwd(wd)
dir.create("summary",recursive=TRUE)

stats_files=c(system("ls screen/*/??_stats.tsv",intern=TRUE), system("ls screen/*/???_stats.tsv",intern=TRUE))


##not included:
##AW, SRP, ISS - not enough data (cross must be bigger, than the number of samples)

all_stats=data.table()
for (file in stats_files){
  f=fread(file)
  all_stats=rbindlist(list(all_stats,f)) 
} 

head(all_stats)
setnames(all_stats,"Species","species")

all_stats_annot=merge(all_stats,stats_annot[,c("species","color_class","scientific_name"),],by="species")
all_stats_annot <- unique(all_stats_annot)

all_stats_annot$group<-factor(unlist(lapply(all_stats_annot$color_class, 
                                            function(x) class_short[x])), levels=class_short)


pdf("summary/model_summary.pdf",width=4, height=4)
#kmer distrib
ggplot(all_stats_annot,aes(x=k,fill=color_class,col=color_class))+geom_bar(alpha=1,size=1)+
  scale_fill_manual(values=class_colors)+scale_color_manual(values=class_colors)+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="", x="kmer length")

##C
ggplot(all_stats_annot,aes(x=as.factor(c),fill=color_class,col=color_class))+geom_bar(alpha=1,size=1)+
  scale_fill_manual(values=class_colors)+scale_color_manual(values=class_colors)+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="", x="C")

#AUC
ggplot(all_stats_annot,aes(x=group,y=AUC, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=class_colors)+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="ROC-AUC", x="")+
  stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4)

#f1
ggplot(all_stats_annot,aes(x=group,y=f1, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=class_colors)+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="f1-score", x="")+
  stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4)

#AUCvsF1
ggplot(all_stats_annot,aes(x=AUC,y=f1, fill=color_class, alpha = 0.7, color = color_class))+
  geom_point(shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=class_colors)+xlim(0.5, 1)+ylim(0.5, 1)+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="f1-score", x="ROC-AUC")

##nseq:
ggplot(all_stats_annot,aes(x=group, y=numSequences, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=unique(all_stats_annot[order(color_class)]$color))+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="numSequences", x="")
dev.off()

pdf("summary/kC.pdf", width=4, height=2)
ggplot(melt(all_stats_annot[, c("c", "k", "color_class")]), aes(x=as.factor(value), fill=color_class)) +
  geom_bar(size=1)+scale_fill_manual(values=class_colors)+
  facet_grid(~variable, scale="free", space="free")+theme(legend.position = "None")+
  labs(x="", y="Number of species")
dev.off()

pdf("summary/auc_for_pres.pdf", width=6, height=4)
ggplot(all_stats_annot,aes(x=color_class,y=AUC, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=unique(all_stats_annot[order(color_class)]$color))+rotate_labels()+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="ROC-AUC", x="")+
  stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4)
 dev.off()

#saving the table 
write.csv(all_stats_annot, "summary/all_aucs.csv", quote = F)

pdf("summary/auc_final.pdf", width=4, height=4)
ggplot(all_stats_annot,aes(x=group,y=AUC, fill=color_class))+
   geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
   scale_color_manual(values=class_colors)+
   theme(legend.position = "None", text=element_text(size=15))+labs(y="ROC-AUC", x="")+
   stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4)
dev.off()


##selected examples, that have AUC closest to the mean value:
selected <- all_stats_annot %>% 
  group_by(color_class) %>%
  mutate(mean_auc=mean(AUC), delta = abs(mean_auc-AUC)) %>%
  filter(delta==min(delta))
colnames(all_stats_annot)
selected[, c("scientific_name", "color_class")]

##draw the ROC-curves for them:
roc_files=system("ls /scratch/lab_bock/shared/projects/compEpi//results_analysis/02_predict_meth/02.1_within_species/screen1/*/roc_res.csv",intern=TRUE)

roc_tables <- lapply(roc_files, read.csv)
sp_list <- unlist(lapply(roc_files, function(x) strsplit(x, "/")[[1]][12] ))

roc_tables <- lapply(seq_along(roc_tables), function(x){roc <- roc_tables[[x]]; 
          roc$species <- sp_list[[x]]; roc$class <- stats_annot[stats_annot$species==sp_list[[x]]]$color_class[[1]];
          return(roc)} )
roc_tables <- lapply(seq_along(roc_tables), function(x){roc<-roc_tables[[x]]; 
roc$id <- paste0(roc_tables[[x]]$run,"_", roc_tables[[x]]$species); 
        roc$color_id<-paste0(roc_tables[[x]]$ifRand,"_", roc_tables[[x]]$class);  return(roc)})
all_rocs <- Reduce(rbind, roc_tables)

head(all_rocs)
new_names <- sapply(names(class_colors), function(x) paste0("noRand_", x))
new_names2 <- sapply(names(class_colors), function(x) paste0("rand_", x))
color_detailed <- c(setNames(unname(class_colors), unname(new_names)), 
                  setNames(rep("lightgray", length(new_names2)), unname(new_names2)))
p_auc_sel <- ggplot(all_rocs, aes(x=fdr, y=tpr, col=color_id))+
  geom_line(aes(group=id, alpha=ifRand))+scale_color_manual(values=color_detailed)+
  scale_alpha_manual(values=c("noRand"=1, "rand"=0.5))+
  labs(x="False Positive Rate", y="True Postitve Rate")

write.csv(all_rocs, "summary/roc_curves_selected.csv", quote = F)

pdf("summary/auc_sel_final.pdf", width=6, height = 4)
p_auc_sel
dev.off()
