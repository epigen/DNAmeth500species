source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
#where to save the output
wd=file.path(analysis_dir,"02_vizStats/02.8_stats_nonCpG_meth")
dir.create(wd)
setwd(wd)


#remove unconverted and other "bad" samples
stats_annot_unconv=stats_annot[grepl("_uc$",Sample_Name)]
#stats_annot=stats_annot[!grepl("_uc$",Sample_Name)&species_check!="fail"&!grepl("Tumour|Cellline",Tissue)&!is.na(Tissue)] 
stats_annot=stats_annot[!grepl("_uc$",Sample_Name)&!grepl("Tumour|Cellline",Tissue)&!is.na(Tissue)] 


## all boxplot cphph
ggplot(stats_annot,aes(x=group,y=avg_meth_cphph,fill=color_class)) +
  geom_boxplot(outlier.shape = 21)+scale_fill_manual(values = class_colors)+
#  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
#  rotate_labels(angle = 60,vjust = 1)+
  ylab("av_meth_cphph")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 15))
ggsave("av_meth_cphph_boxplots.pdf", height=4,width=6)

## all boxplot cphpg
ggplot(stats_annot,aes(x=group,y=avg_meth_cphpg,fill=color_class)) +
  geom_boxplot(outlier.shape = 21)+scale_fill_manual(values = class_colors)+
#  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
#  rotate_labels(angle = 60,vjust = 1)+
  ylab("av_meth_cphpg")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 15))
ggsave("av_meth_cphpg_boxplots.pdf", height=4,width=6)

## correlations:

#cphpg-cphph
ggplot(stats_annot,aes(x=avg_meth_cphph,y=avg_meth_cphpg,fill=color_class)) +
  geom_point(shape = 21)+scale_fill_manual(values = class_colors)+
#  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
#  rotate_labels(angle = 60,vjust = 1)+
 # ylab("av_meth_cphpg")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 15))
#ggsave("av_meth_cphpg_boxplots.pdf", height=4,width=6)


### per tissue:

#core tissues
core_inv=c("Muscle","Tube_feet","Tentacle","Arm","Gills","Gonad","Pharynx")
core_vert=c("Liver","Heart","Brain","Spleen","Muscle","Gills","Fin")

ggplot(stats_annot[Tissue %in% core_vert],aes(x=group,y=avg_meth_cphph,fill=color_class)) +
  geom_boxplot(outlier.shape = 21)+scale_fill_manual(values = class_colors)+
#  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
#  rotate_labels(angle = 60,vjust = 1)+
  ylab("av_meth_cphph")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 15)) + 
    facet_wrap(~Tissue, ncol = 2)
ggsave("av_meth_cphph_boxplots_pertissue.pdf", height=12,width=6)


### brain vs others

stats_annot[, brain_tissue := ifelse(Tissue=="Brain", "Brain", "all other"),]

ggplot(stats_annot[Tissue %in% core_vert],aes(x=group,y=avg_meth_cphph, fill=brain_tissue)) +
  geom_boxplot(outlier.shape = 21)+
#  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
#  rotate_labels(angle = 60,vjust = 1)+
  ylab("av_meth_cphph")+xlab("")+
  theme(text = element_text(size = 15))