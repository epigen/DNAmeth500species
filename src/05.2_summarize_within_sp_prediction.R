source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"05_predict_meth","05.1_within_species")
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

all_stats_annot=merge(all_stats,unique(stats_annot[,c("species","color_class","group"),]),by="species")

all_stats_annot <- unique(all_stats_annot)

pdf("summary/kC.pdf", width=4, height=2)
ggplot(melt(all_stats_annot[, c("c", "k", "color_class")]), aes(x=as.factor(value), fill=color_class)) +
  geom_bar(size=1)+scale_fill_manual(values=class_colors)+
  facet_grid(~variable, scale="free", space="free")+theme(legend.position = "None")+
  labs(x="", y="Number of species")
dev.off()


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

#AUCvsF1
ggplot(all_stats_annot,aes(x=AUC,y=f1, fill=color_class, alpha = 0.7, color = color_class))+
  geom_point(shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=class_colors)+xlim(0.5, 1)+ylim(0.5, 1) +
  theme(legend.position = "None", text=element_text(size=15))+labs(y="f1-score", x="ROC-AUC")

#AUC
all_stats_annot[group=="Jl.vb.",group:= "Inv.", ]
all_stats_annot[color_class=="Jawless_vertebrate",color_class:= "Invertebrata", ]

ggplot(all_stats_annot,aes(x=group,y=AUC, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=class_colors)+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="ROC-AUC", x="")+
  stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4) + 
  geom_point(data = all_stats_annot[species=="JL"],aes(x=group,y=AUC), fill=class_colors["Jawless_vertebrate"], shape = 21)

#f1
ggplot(all_stats_annot,aes(x=group,y=f1, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=class_colors)+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="f1-score", x="")+
  stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4) + 
  geom_point(data = all_stats_annot[species=="JL"],aes(x=group,y=f1), fill=class_colors["Jawless_vertebrate"], shape = 21)


##nseq:
ggplot(all_stats_annot,aes(x=group, y=numSequences, fill=color_class))+
  geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
  scale_color_manual(values=unique(all_stats_annot[order(color_class)]$color))+
  theme(legend.position = "None", text=element_text(size=15))+labs(y="numSequences", x="") + 
  geom_point(data = all_stats_annot[species=="JL"],aes(x=group,y=numSequences ), fill=class_colors["Jawless_vertebrate"], shape = 21)
dev.off()


#saving the table 
write.csv(all_stats_annot, "summary/all_aucs.csv", quote = F)
all_stats_annot <- fread("summary/all_aucs.csv")

pdf("summary/auc_final.pdf", width=4, height=4)
ggplot(all_stats_annot,aes(x=group,y=AUC, fill=color_class))+
   geom_boxplot(outlier.shape=21)+scale_fill_manual(values=class_colors)+
   scale_color_manual(values=class_colors)+
   theme(legend.position = "None", text=element_text(size=15))+labs(y="ROC-AUC", x="")+
   stat_summary(fun.data = give.n,fun.args = c(y=0.4), geom = "text",size=4) + 
  geom_point(data = all_stats_annot[species=="JL"],aes(x=group,y=AUC), fill=class_colors["Jawless_vertebrate"], shape = 21)
dev.off()


###AUC plot for JL
JL_stats <- all_stats_annot
JL_stats[,plot_group:=ifelse(species=="JL","Lamprey",ifelse(color_class=="Invertebrata","Invertebrata","Vertebrata")),]


#cap outliers for better plotting
#JL_stats[,AUC:=ifelse(value>quantile(value,0.99),quantile(value,0.99),value),by=c("variable")]

pdf("summary/Lamprey_auc.pdf",2,2)
ggplot(JL_stats,aes(x=AUC,fill=plot_group,col=plot_group))+geom_density(lwd=1,alpha=0.3)+
  geom_point(data=JL_stats[species=="JL"],y=0, size=3) + xlim(c(0.4, 1)) + theme(legend.position = "None")
dev.off()

##selected examples, that have AUC closest to the mean value:
selected <- all_stats_annot %>% 
  group_by(color_class) %>%
  mutate(mean_auc=mean(AUC), delta = abs(mean_auc-AUC)) %>%
  filter(delta==min(delta))
colnames(all_stats_annot)
left_join(selected[, c("species", "color_class")], unique(stats_annot[, c("species", "scientific_name")]))

getwd()
roc_file <- fread("screen/BD/roc_res.csv")

##draw the ROC-curves for them:
roc_tables <- lapply(selected$species, function(x) fread(file.path("screen", x,"roc_res.csv")))
                     

roc_tables <- lapply(seq_along(roc_tables), function(x){roc <- roc_tables[[x]]; 
          roc$species <- selected$species[[x]]; roc$class <- stats_annot[stats_annot$species==selected$species[[x]]]$color_class[[1]];
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
                     
                     
### plotting k-AUC summary
                     
kauc_df=data.table()
                     
files_list=c(system("ls screen/*/??_AUC_k_grid.tsv",intern=TRUE), system("ls screen/*/???_AUC_k_grid.tsv",intern=TRUE))

for(file in files_list){
    temp_df <- fread(file)
    temp_df$species <- strsplit(file,"/")[[1]][2]
    kauc_df = rbind(kauc_df, temp_df)
}
my_wt(kauc_df, "k_AUC_grid.tsv")

kauc_mt <- kauc_df %>% pivot_wider(names_from = species,values_from = GridCol_1)
kauc_mt <- as.data.frame(kauc_mt)
row.names(kauc_mt) <- kauc_mt$k
kauc_mt <- kauc_mt[, -1]

sp_order <- sp_df[species %in%colnames(kauc_mt)]$species
kauc_mt <- kauc_mt[,sp_order]

annot_row = sp_df[species %in% sp_order, c("color_class"), drop = F] 
row.names(annot_row) <- sp_order

                     
p <- pheatmap(kauc_mt, cluster_rows = F, scale = "column",
         cluster_cols = F, show_colnames = F, main="AUC grid from k", 
        annotation_col = annot_row, annotation_color = list(color_class = class_colors))
                     
                     
                     
ggsave("summary/Heatmap_k_AUC_grid_nonscaled.pdf",plot = p, width = 7, height = 3)
