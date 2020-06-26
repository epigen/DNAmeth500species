source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(maps)

wd=file.path(analysis_dir,"01_basicStats/01.6_stats_detail")
dir.create(wd,recursive=TRUE)
setwd(wd)

#######setup######
#class order
stats_annot[ncbi_order%in%c("Diprotodontia","Dasyuromorphia"),color_class:="Marsupialia",]
stats_annot[,color_class:=factor(color_class,levels=names(class_colors)),] 

#adding short names, rhat we will use for figure panels
stats_annot$group<-factor(unlist(lapply(stats_annot$color_class, function(x) class_short[x])),
                          levels=class_short)

#core tissues
core_inv=c("Muscle","Tube_feet","Tentacle","Arm","Gills","Gonad","Pharynx")
core_vert=c("Liver","Heart","Brain","Spleen","Muscle","Gills","Fin")


#remove unconverted and other "bad" samples
stats_annot_unconv=stats_annot[grepl("_uc$",Sample_Name)]
stats_annot=stats_annot[!grepl("_uc$",Sample_Name)&species_check!="fail"&!grepl("Tumour|Cellline",Tissue)&!is.na(Tissue)] 

#fix sex and age
stats_annot[sex=="male",sex:="m",]
stats_annot[is.na(sex), sex:="not specified",]
stats_annot$sex=factor(stats_annot$sex, levels = c("not specified", "f", "m"))
stats_annot[age%in%c("juv","juvenil","juvenile"),age:="young"]
stats_annot[age%in%c("1 j","2 j"),age:=NA]
stats_annot[is.na(age), age:="not specified"]
stats_annot$age=factor(stats_annot$age, levels=c("not specified", "young", "adult", "old"))

#fix duplicate entry for on species (ACT), this is due to duplicte entry in patho annotation --> fix there
stats_annot=stats_annot[!duplicated(Sample_Name,fromLast=TRUE)]

#check deduced reference for biases
ref_stats=fread(file.path(analysis_dir,"01_basicStats/ref_stats_long.tsv"))
classes=unique(stats_annot[,c("species","color_class","group"),])
ref_stats=merge(ref_stats,classes,by="species")

pdf("ref_stats.pdf",height=7,width=11)
ggplot(ref_stats,aes(x=group,y=value,fill=color_class))+geom_boxplot()+
  facet_wrap(~variable,scale="free")+scale_fill_manual(values = class_colors)
dev.off()


#preparations
stats_annot[,c("N_cs","med_meth_cs"):=list(.N,median(CpG_meth)),by=color_class]
stats_annot[,c("N_os","med_meth_os"):=list(.N,median(CpG_meth)),by=ncbi_order]
stats_annot[,c("N_tissue"):=list(.N),by=Tissue]
stats_annot[,c("N_ss"):=list(.N),by="species"]
stats_annot[,c("N_rep", "N_tissue_per_sp"):=list(length(unique(replicate)),length(unique(Tissue))),by="species"]

stats_summary_data <- do.call(cbind, lapply(stats_annot[, c("N_ss", "N_rep", "N_tissue_per_sp")], summary))

write.csv(stats_summary_data,  "stats_summary_per_speices.csv", quote = F)

#for each class find species with most samples
stats_annot[,unique(species[order(N_ss,decreasing=TRUE)])[1:10],by=color_class]

##Bubble heatmap, showing per class per tissue stats:
tissue_count <- unique(stats_annot[, c("Tissue", "N_tissue")])
tissue_count <- tissue_count[order(tissue_count$N_tissue, decreasing = T),]
tissue_count_per_class <-stats_annot %>% 
  group_by(color_class,Tissue) %>%
  filter(Tissue %in% unique(stats_annot[N_tissue>=50]$Tissue)) %>%
  summarize(n=n())

other_per_class <- stats_annot %>% 
  filter(Tissue %in% unique(stats_annot[N_tissue<50]$Tissue)) %>%
  group_by(color_class) %>%
  summarize(n=n()) %>%
  mutate(Tissue = "other")

tissue_count_per_class <- as.data.frame(bind_rows(tissue_count_per_class, other_per_class))
tissue_count_per_class$Tissue <- factor(tissue_count_per_class$Tissue, 
                                        levels = c(tissue_count$Tissue, "other"))
my_wt(tissue_count_per_class, "tissue_stats.tsv")

ggplot(tissue_count_per_class, aes(y = Tissue, x = color_class)) + geom_point(aes(color = n, size = n)) + 
  rotate_labels() + coord_equal() + scale_color_gradient(low = "#6baed6", high = "#08306b")
ggsave("bubble_heatmap.pdf", width = 6, height = 6)

#basic stats
##plot count of sanoles per tissue, stacked by color class
pdf("basic_tissues_class.pdf",height=4,width=8)
stats_annot[,Tissue:=factor(Tissue,levels=unique(Tissue[order(N_tissue,decreasing=TRUE)])),]
ggplot(stats_annot,aes(x=Tissue,fill=color_class))+geom_bar()+scale_fill_manual(values = class_colors)+
  rotate_labels(angle = 60,vjust = 1)+ylab("Number of samples") + theme(legend.position = "None")
dev.off()

##plot count of samples per tissue, stacked by color class
pdf("basic_tissues_class_core.pdf",height=4,width=5)
ggplot(stats_annot[N_tissue>=10],aes(x=Tissue,fill=color_class)) +
  geom_bar()+scale_fill_manual(values = class_colors) +
  rotate_labels(angle = 60,vjust = 1)+ylab("Number of samples") + 
  theme(legend.position = "None",  text = element_text(size=15))
dev.off()

##number of samples per class
pdf("basic_samples_class.pdf",height=4,width=4)
ggplot(stats_annot,aes(x=group,fill=color_class)) + geom_bar() +
  scale_fill_manual(values = class_colors)+
  ylab("Number of samples")+theme(legend.position = "None", text = element_text(size=15))
dev.off()


pdf("basic_samples_class_pres.pdf",height=4,width=4)
ggplot(stats_annot,aes(x=color_class,fill=color_class)) +
  geom_bar()+scale_fill_manual(values = class_colors) +
  ylab("Number of samples") + theme(legend.position = "None") + xlab("") + rotate_labels()
dev.off()


sex_colors=c("f"="#f28e00","m"="#22b900", "not specified"="lightgrey")
age_colors=c("old"="#3182bd","adult"="#9ecae1","young"="#deebf7","not specified"="lightgrey")

pdf("basic_samples_class_sex_age.pdf",height=4,width=5)
ggplot(stats_annot,aes(x=group,fill=sex))+geom_bar(position = "fill") +
  ylab("Fraction of samples")+xlab("")+theme(legend.position = "bottom", text = element_text(size=15)) +
  scale_fill_manual(values=sex_colors)
ggplot(stats_annot,aes(x=group,fill=age))+geom_bar(position = "fill") +
  ylab("Fraction of samples")+xlab("")+theme(legend.position = "bottom", text = element_text(size=15)) +
  scale_fill_manual(values = age_colors)
dev.off()



orders = stats_annot[,.(N_orders=length(unique(ncbi_order)),N_species=length(unique(ncbi_name)),
                        N_tissues=length(unique(Tissue))),by=c("color_class")]
orders$group <- factor(unlist(lapply(orders$color_class, 
                                     function(x) class_short[x])), levels=class_short)


pdf("basic_phylo_class.pdf",height=4,width=4)
ggplot(orders,aes(x=group,y=N_orders,fill=color_class))+geom_bar(stat="identity") + 
  scale_fill_manual(values = class_colors)+
  xlab("")+ylab("Number of different orders")+theme(legend.position = "None",  text = element_text(size=15))
ggplot(orders,aes(x=group,y=N_species,fill=color_class))+geom_bar(stat="identity") + 
  scale_fill_manual(values = class_colors)+
  xlab("")+ylab("Number of different species")+theme(legend.position = "None",  text = element_text(size=15))
ggplot(orders,aes(x=group,y=N_tissues,fill=color_class)) + geom_bar(stat="identity") + 
  scale_fill_manual(values = class_colors)+
  xlab("")+ylab("Number of different tissues")+theme(legend.position = "None",  text = element_text(size=15))
dev.off()




#CpG methylation
pdf("CpG_meth_all_for_pres.pdf",height=4,width=6)
ggplot(stats_annot,aes(x=color_class,y=CpG_meth,fill=color_class)) +
  geom_boxplot(outlier.shape = 21)+scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
  rotate_labels(angle = 60,vjust = 1)+
  ylim(c(0,100))+ylab("% DNA methylation (CpG)")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 15))
dev.off()

#CpG methylation

ggplot(stats_annot,aes(x=group,y=CpG_meth,fill=color_class))+geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=2)+
  ylim(c(0,100))+ylab("% DNA methylation (CpG)")+xlab("")+theme(legend.position = "None")+
  theme(text = element_text(size = 10))
ggsave("CpG_meth_all.pdf",height=7,width=7, units = "cm")


pdf("CpG_meth_tissue.pdf",height=4,width=12)
ggplot(stats_annot[Tissue%in%core_vert],aes(x=Tissue,y=CpG_meth,fill=color_class))+
  geom_boxplot()+scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=2.5)+
  rotate_labels(angle = 60,vjust = 1)+ylim(c(0,100))+facet_wrap(~color_class,ncol=4)+
  ylab("% DNA methylation (CpG)")+xlab("")
dev.off()

#COV across different levels of abstraction
stats_annot_coreTissues=stats_annot[Tissue%in%c(core_inv,core_vert)]

COV_class=stats_annot_coreTissues[,.(COV=sd(CpG_meth)/mean(CpG_meth),N=.N,type="groups"),by=c("color_class")]
COV_spec=stats_annot_coreTissues[,.(COV=sd(CpG_meth)/mean(CpG_meth),N=.N,type="species"),by=c("abbreviation_sp")]
COV_tissues=stats_annot_coreTissues[,.(COV=sd(CpG_meth)/mean(CpG_meth),N=.N,type="tissues"),by=c("Tissue")]
COV_repl=stats_annot_coreTissues[,.(COV=sd(CpG_meth)/mean(CpG_meth),N=.N,type="species+\nreplicate"),by=c("abbreviation_sp","replicate")]
COV_st=stats_annot_coreTissues[,.(COV=sd(CpG_meth)/mean(CpG_meth),N=.N,type="species+\ntissue"),by=c("abbreviation_sp","Tissue")]

COV=rbindlist(list(COV_class,COV_spec,COV_tissues,COV_repl[,-c("abbreviation_sp")],COV_st[,-c("abbreviation_sp")]))
pdf("CpG_meth_COV.pdf",height=2.5,width=4)
ggplot(COV,aes(x=type,y=COV))+geom_boxplot()+stat_summary(fun.data = give.n,fun.args = c(y=-0.03), geom = "text",size=2.5)+ylab("COV % DNA methylation")+xlab("")
dev.off()


#focus groups
stats_annot[,ncbi_order:=factor(ncbi_order,levels=unique(ncbi_order[order(med_meth_os)])),]

pdf("CpG_meth_mammals.pdf",height=5,width=5)
ggplot(stats_annot[color_class%in%c("Mammalia","Marsupialia")&N_os>5&!is.na(ncbi_order)],aes(x=ncbi_order,y=CpG_meth,col=color_class))+
  geom_boxplot(fill="transparent",outlier.shape = NA)+
  geom_point(alpha=0.5,pch=21,position=position_jitter())+
  stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",angle=90,hjust=0)+rotate_labels(angle = 60,vjust = 1)+ylim(c(0,100))+scale_color_manual(values=c("Mammalia"="black","Marsupialia"="red"))+ylab("% DNA methylation (CpG)")+xlab("")
dev.off()

pdf("CpG_meth_invertebrata.pdf",height=5,width=7)
sub_inv=stats_annot[color_class%in%c("Invertebrata")&!is.na(ncbi_order)&N_os>1]
sub_inv[,Tissue:=ifelse(Tissue%in%core_inv,Tissue,"other"),]
ggplot(sub_inv,aes(x=ncbi_order,y=CpG_meth))+
  geom_boxplot(fill="transparent",outlier.shape = NA)+
  geom_point(alpha=0.5,size=2.5,pch=21,aes(fill=Tissue),position=position_jitter())+
  stat_summary(fun.data = give.n,fun.args = c(y=-10), geom = "text",angle=90,hjust=0)+rotate_labels(angle = 60,vjust = 1)+ylim(c(-10,100))+ylab("% DNA methylation (CpG)")+xlab("")+scale_fill_manual(values = col_vector[-c(1,2)])
dev.off()


#base composition
base_compo=melt(stats_annot[,c("perc_As", "perc_Ts", "perc_Cs", "perc_Gs", "Sample_Name","ncbi_order","color_class","N_os", "group"),],id.vars = c("Sample_Name","ncbi_order","color_class","N_os", "group"),variable.name = "Base")

base_compo_mean=base_compo[,.(value=mean(value),N_os=unique(N_os)),by=c("Base","color_class","ncbi_order")]


pdf("bases_orders.pdf",height=4,width=14)
ggplot(base_compo_mean,aes(x=ncbi_order,y=value,fill=Base))+geom_bar(stat="identity")+rotate_labels(angle = 60,vjust = 1)+geom_text(aes(label=paste0("N=",N_os)),angle=90,y=0,hjust=0,vjust=0.5,size=2)+facet_grid(~color_class,scale="free_x",space = "free_x")+ylim(c(0,100))+ylab("%")
dev.off()

pdf("bases_classes.pdf",height=5,width=6)
ggplot(base_compo,aes(x=color_class,y=value,col=Base))+geom_boxplot()+rotate_labels(angle = 60,vjust = 1)+stat_summary(fun.data = give.n,fun.args = c(y=75), geom = "text",angle=90,hjust=0,size=2,position=position_dodge(width=1))+ylim(c(0,100))+ylab("%")
dev.off()

pdf("bases_mammalia.pdf",height=5,width=6)
ggplot(base_compo[color_class%in%c("Mammalia","Marsupialia")&!is.na(ncbi_order)&N_os>2],aes(x=ncbi_order,y=value,color=color_class))+
  geom_boxplot()+rotate_labels(angle = 60,vjust = 1)+
  facet_wrap(~Base,scale="free_y")+
  stat_summary(fun.data = give.n, geom = "text",angle=90,hjust=-0.5,size=3)+
  rotate_labels(angle = 60,vjust = 1)+ylim(c(0,100))+scale_color_manual(values=c("Mammalia"="black","Marsupialia"="red"))+ylab("%")
dev.off()

pdf("base_composition_by_tax.pdf", width = 5, height = 3.5)
ggplot(base_compo,aes(x=group,y=value))+
  geom_boxplot(aes(fill=Base), outlier.shape = 21)+
  stat_summary(fun.data = give.n, 
               fun.args = c(y=60), 
               geom = "text", angle = 90,
               hjust=0,size=3, position=position_dodge(width=1))+
  ylim(c(0,70))+ rotate_labels()+
  ylab("%") + xlab("") + 
  theme( text = element_text(size = 13))
dev.off()

##nuber of ded_ref_fragments (i.e. genome size)
if(file.exists("dedRef_CpG_count.csv")){
  dedRef_count <- read.csv("dedRef_CpG_count.csv", header = 1, sep = ";")
  dedRef_count$species <- row.names(dedRef_count)
  dedRef_count <- unique(inner_join(dedRef_count, stats_annot[,c("species","color_class", "group")]))

###comparing parameters
pdf("full_vs_filtered_CpG_count.pdf", width = 10, height = 5)
ggplot(dedRef_count, aes(x = dedRef_count_full, y = dedRef_count_filtered, fill = color_class)) + 
  geom_point(shape = 21, alpha = 0.5) + scale_fill_manual(values = class_colors)+
  ylab("filtered number of dedRef fragments")+xlab("full number of dedRef fragments")+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(label=scientific_10) + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.5) + 
  facet_wrap(~color_class, ncol = 4) + coord_fixed()

ggplot(dedRef_count, aes(x = CpG_count, y = CpG_count_filtered, fill = color_class)) + 
  geom_point(shape = 21, alpha = 0.5) + scale_fill_manual(values = class_colors) +
  ylab("filtered number of CpGs")+xlab("full number of CpGs") +
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(label=scientific_10) + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.5) + 
  facet_wrap(~color_class, ncol = 4)

ggplot(dedRef_count, aes(x = dedRef_count_filtered, y = CpG_count, fill = color_class)) + 
  geom_point(shape = 21, alpha = 0.5) + scale_fill_manual(values = class_colors)+
  ylab("number of CpGs") + xlab("filtered number of dedRef fragments") +
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(label=scientific_10) + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.5) + 
  facet_wrap(~color_class, ncol = 4)
dev.off()

###comparing parameters
pdf("full_vs_filtered_CpG_count_merged.pdf", width = 3.5, height = 3.5)
ggplot(dedRef_count, aes(x = dedRef_count_full, y = dedRef_count_filtered, fill = color_class)) + 
  geom_point(shape = 21, alpha = 0.5) + scale_fill_manual(values = class_colors)+
  ylab("filtered number of dedRef fragments")+xlab("full number of dedRef fragments")+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(label=scientific_10) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.5) 

ggplot(dedRef_count, aes(x = CpG_count, y = CpG_count_filtered, fill = color_class)) + 
  geom_point(shape = 21, alpha = 0.5) + scale_fill_manual(values = class_colors)+
  ylab("filtered number of CpGs")+xlab("full number of CpGs")+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(label=scientific_10) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.5) 

ggplot(dedRef_count, aes(x = dedRef_count_filtered, y = CpG_count, fill = color_class)) + 
  geom_point(shape = 21, alpha = 0.5) + scale_fill_manual(values = class_colors)+
  ylab("number of CpGs")+xlab("filtered number of dedRef fragm.")+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10) + 
  scale_x_continuous(label=scientific_10) + coord_fixed() + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', alpha = 0.5) 

dev.off()


##genome stats (for supplementary)
pdf("genome_stats.pdf",height=3.5,width=3.5)

ggplot(dedRef_count, aes(x = group, y = dedRef_count_full, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0.85*max(dedRef_count$dedRef_count_full)), 
               geom = "text",size=3, angle = 90, hjust=0,
               position=position_dodge(width=1))+
  ylab("number of dedRef fragm.")+xlab("")+theme(legend.position = "None", text = element_text(size = 13)) + 
  rotate_labels() + scale_y_continuous(label=scientific_10) 

ggplot(dedRef_count, aes(x = group, y = dedRef_count_filtered, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0.85*max(dedRef_count$dedRef_count_filtered)), 
               geom = "text",size=3, angle = 90, hjust=0, position=position_dodge(width=1))+
  ylab("number of dedRef fragm.(filtered)")+xlab("") + 
  theme(legend.position = "None", text = element_text(size = 13)) + 
  rotate_labels() + scale_y_continuous(label=scientific_10)

ggplot(dedRef_count, aes(x = group, y = CpG_count, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0.8*max(dedRef_count$CpG_count)), 
               geom = "text",size=3, angle = 90, hjust=0, position=position_dodge(width=1))+
  ylab("count of GpGs")+xlab("")+theme(legend.position = "None", text = element_text(size = 13)) + 
  rotate_labels() + scale_y_continuous(label=scientific_10)

ggplot(dedRef_count, aes(x = group, y = CpG_count_filtered, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=0.8*max(dedRef_count$CpG_count_filtered)), 
               geom = "text",size=3, angle = 90, hjust=0, position=position_dodge(width=1))+
  ylab("count of GpGs (filtered)")+xlab("")+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  rotate_labels() + scale_y_continuous(label=scientific_10)

dev.off()
}else{
  print("dedRef_CpG_count file not found. Please, run the script 01.61 to generate it")
}

###sequencing stats (for supplementary)
pdf("sequencing_stats_conv.pdf",height=3.5,width=3)
#number of covered CpGs
ggplot(stats_annot[stats_annot$conversion_type == "converted",], 
       aes(x = group, y = coveredCpGs, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,fun.args = c(y=1.01*max(stats_annot$coveredCpGs)),
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+ 
  ylab("number of covered CpGs") + xlab("") + rotate_labels() +
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(label=scientific_10, limits = c(0, 1.15*max(stats_annot$coveredCpGs)))  

#mapping rate
ggplot(stats_annot[stats_annot$conversion_type == "converted",], 
       aes(x = group, y = mapping_efficiency, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,
               fun.args = c(y=1.01*max(stats_annot$mapping_efficiency)),
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+ 
  ylab("mapping efficiency")+xlab("")+rotate_labels()+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(
          limits = c(0, 1.15*max(stats_annot$mapping_efficiency)))  

#pre-fragmentation
ggplot(stats_annot[stats_annot$conversion_type == "converted",], 
       aes(x = group, y = others, fill = color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors) +
  stat_summary(fun.data = give.n,
               fun.args = c(y=1.01*max(stats_annot$others)), 
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+
  ylab("% prefragnemtation")+xlab("")+rotate_labels()+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(
    limits = c(0, 1.15*max(stats_annot$others)))  

#contamination rate
ggplot(stats_annot[stats_annot$conversion_type == "converted",],
       aes(x=group,y=cont_rat,fill=color_class))+
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors)+
  stat_summary(fun.data = give.n,
               fun.args = c(y=1.01*max(stats_annot$cont_rat)), 
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+
  ylab("% contamination")+xlab("")+rotate_labels()+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(
    limits = c(0, 1.15*max(stats_annot$cont_rat)))  

#enrichment cycles
ggplot(stats_annot[stats_annot$conversion_type == "converted",],
       aes(x=group,y=`Enrichment cycles`,fill=color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors) + 
  stat_summary(fun.data = give.n,
               fun.args = c(y=1.01*max(stats_annot$`Enrichment cycles`)), 
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+
  ylab("PCR enrichment cycles") + xlab("") + rotate_labels()+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(
    limits = c(5, 1.1*max(stats_annot$`Enrichment cycles`)))

#conversion rate
ggplot(stats_annot[stats_annot$conversion_type == "converted",],
       aes(x=group,y=conversionRate,fill=color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors) + 
  stat_summary(fun.data = give.n,
               fun.args = c(y=1.01*max(stats_annot$conversionRate)), 
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+
  ylab("conversion rate") + xlab("") + rotate_labels()+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(
    limits = c(70, 1.05*max(stats_annot$conversionRate)))

dev.off()

stats_conv <- stats_annot[, c("conversionRate", "k1_unmeth","color_class")]
stats_conv$k1_unmeth <- 100*stats_conv$k1_unmeth
stats_conv$k1_conv_rate <- 100 - stats_conv$k1_unmeth
stats_conv[,best_conv_rate:=pmax(conversionRate,k1_conv_rate),]
                                     
stats_conv$risk <- c( "FALSE")
stats_conv[stats_conv$conversionRate > 98, ]$risk <- c( "TRUE")
stats_conv[stats_conv$k1_unmeth < 2, ]$risk <- c( "TRUE")

##additional convertion plots:
ggplot(stats_conv, aes( y = k1_unmeth, x = conversionRate, fill = risk) )+ 
  geom_point(shape = 21) + xlim(c(0, 100)) + ylim(c(0, 100)) + 
  facet_wrap(~color_class, nrow = 1) + 
  geom_vline(xintercept = 98, linetype = "dashed", color = "grey", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey", alpha = 0.5) + 
  scale_fill_manual(values = c("TRUE" = "grey", "FALSE" = "red")) + coord_fixed() +
  labs(x = "% converted Cs", y = "% spike-in methylation", 
       fill = "Conversion rate > 98% or\nUnmethylated spike-in < 2%") + 
  theme(legend.position = "bottom", text = element_text(size = 13))
ggsave("conv_stats.pdf", width = 15, height = 4)
                                     
#conversion plot showing best of non-CpG C conversion rate and k1_unmeth
pdf("../max_conv_rate.pdf",height=3.5,width=3)      #put into top level folder because I don't have permission for the actual folder                               
ggplot(stats_conv,aes(x=color_class,y=best_conv_rate,fill=color_class)) + 
  geom_boxplot(outlier.shape = 21) + scale_fill_manual(values = class_colors) + 
  stat_summary(fun.data = give.n,
               fun.args = c(y=100.3), 
               geom = "text",size=3, angle = 90, hjust=0, 
               position=position_dodge(width=1))+                                    
  ylab("Conversion efficiency (%)") + xlab("") + rotate_labels()+
  theme(legend.position = "None", text = element_text(size = 13)) + 
  scale_y_continuous(
    limits = c(95, 101))+scale_x_discrete(labels=class_short)
dev.off()                                    
                                     
                                     
                                     
###########plotting onto tree################################

library(ape)
library(ggtree)
library(ggstance)
#library(rphylopic)
#library(ggimage)


mammalia=stats_annot[color_class%in%c("Mammalia","Marsupialia")&!is.na(ncbi_order)]
#write(paste0(as.character(unique(mammalia$ncbi_order)),"\n"),file = "mammalia.txt") # not actually need (would be needed for NCBI taxonomy)
#get ded_ref base frequencies
ded_ref_bases=fread(paste0(analysis_dir,"/01_basicStats/feature_summary_filtered.tsv"))

mammalia=merge(mammalia,ded_ref_bases,by.y="species",by.x="abbreviation_sp")

#get Mammalia tree from http://www.timetree.org/ in section "BUILD A TIMETREE" an save as mammalia_order.nwk
tree=read.tree(paste0(meta_dir,"/mammalia_order.nwk")) 

#build tree
p=ggtree(tree)+theme_tree2()+xlim_expand(c(0,280), panel = "Tree")+geom_tiplab()
#plot tree with node labels (to select nodes to plot pcs)
p+geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

#create pic data
pic=data.table(id=c("Diprotodontia",
                    "Dasyuromorphi",
                    "Cingulata",
                    "Proboscidea",
                    "Rodontia",
                    "Primates",
                    "Carnivora"),
               species=c("Macropus eugenii",
                         "Sarcophilus harrisii",
                         "Dasypus novemcinctus",
                         "Loxodonta",
                         "Mus musculus",
                         "Pongo",
                         "Canis lupus familiaris"),
               uid=c("821f487d-9d47-4a36-ac51-8782cf2c7c93",
                     "3d7479e2-0708-4336-9a3c-136ed9a0b581",
                     "61a5cea9-e8e6-419a-8185-68669b6abafb",
                     "6e47a104-e62a-4901-81f7-fca2ecacdb0e",
                     "3b89954e-b28e-4062-89cd-1ad5df0eb431",
                     "63c557ce-d82c-42e6-a26a-a9f0f05c2c18",
                     "e4e306cd-73b6-4ca3-a08c-753a856f7f12"),
               node=c(31,34,47,48,43,45,41))
#pic[,uid:=search_text(text=species)$uid[1],by=1:nrow(pic)] #good idea but yiels many invalid uids --> one needs to find UID by clicking on "copy image address"

#create label annotation data
d1 <- data.table(id=tree$tip.label)
d1[,present:=ifelse(id%in%unique(mammalia$ncbi_order),"True","False"),]

#create data object for base frequencies
d2 <- melt(mammalia[,c("perc_As", "perc_Ts", "perc_Cs", "perc_Gs", "ncbi_order"),],id.vars = c("ncbi_order"),variable.name = "Base")
d2=d2[,.(N=.N,value=mean(value)),by=c("ncbi_order", "Base")]
d2[,Base:=gsub("perc_|s","",Base),]
setnames(d2,c("ncbi_order"),c("id"))

#create data object for base frequencies (deduced reference)
d2.1 <- melt(mammalia[,c("A", "T", "C", "G", "ncbi_order"),],id.vars = c("ncbi_order"),variable.name = "Base")
d2.1=d2.1[,.(N=.N,value=mean(value)),by=c("ncbi_order", "Base")]
d2.1[,Base:=as.character(Base),]
setnames(d2.1,c("ncbi_order"),c("id"))

#create data object for CpG methylation
d3 <-mammalia[,c("ncbi_order","CpG_meth"),]
setnames(d3,"ncbi_order","id")

#add everything to the tree
p1 <- p %<+% d1 + geom_tiplab(aes(color=present))+scale_color_manual(values = c("True"="black","False"="grey"))

p1 <- facet_plot(p1, panel="% DNA methylation", data=d3, geom=geom_boxploth, aes(x=CpG_meth,group=label),width=0.5)+xlim_expand(c(0,100), panel = "% DNA methylation")

p1 <- facet_plot(p1, panel="% DNA methylation", data=unique(d2[,c("id","N"),]), geom=geom_text, hjust= 0,mapping=aes(x=0, label=paste0("N=",N)))

p1 <- facet_plot(p1, panel="Bases (library)", data=d2, geom=geom_barh, mapping=aes(x=value,fill=Base),stat="identity")
p1 <- facet_plot(p1, panel="Bases (reference)", data=d2.1, geom=geom_barh, mapping=aes(x=value,fill=Base),stat="identity")


#p1 <- p1 %<+% pic[1:3] + geom_nodelab(aes(image=uid),size=0.1, geom="phylopic", alpha=.5, color='steelblue')

pdf("../mammalia_meth_bases.pdf",height=5,width=15)
p1  +xlim_expand(c(0,100), panel = "% DNA methylation")+ theme(legend.position="right")
dev.off()


#### size of the genome (dedRef fragments)

files <- system(paste0("ls ",processed_dir,"/*/reduced/consensus/toSelf_filtered_0.08mm_final.fa"), intern=TRUE)
counts <- sapply(files, function(x) as.numeric(system(paste0("wc -l <", x), intern = T)))
species <- sapply(names(counts), function(x) strsplit(x, "/")[[1]][8])

counts <- setNames(counts, species)
counts <- counts/2

write.csv(as.data.frame(counts), "number_of_dedRef_fragments.csv", quote = F)
