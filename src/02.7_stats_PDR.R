source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
#where to save the output
wd=file.path(analysis_dir,"02_vizStats/02.7_stats_PDR")
dir.create(wd)
setwd(wd)



## mean PDR 

ggplot(stats_annot,aes(x=group,y=PDR_mean,fill=color_class)) +
  geom_boxplot(outlier.shape = 21)+scale_fill_manual(values = class_colors)+
#  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
#  rotate_labels(angle = 60,vjust = 1)+
  ylab("mean PDR")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 15))
ggsave("PDR_boxplots.pdf", height=4,width=6)