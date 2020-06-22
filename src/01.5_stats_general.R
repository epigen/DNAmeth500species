source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(tableHTML)

wd = file.path(analysis_dir,"01_basicStats")
setwd(wd)


stats_annot=fread("all_stats.tsv")
stats_long_annot=fread("all_stats_long.tsv")
ref_long=fread("ref_stats_long.tsv")
#remove unconverted samples
stats_annot=stats_annot[grep("_uc",Sample_Name,invert=TRUE)]
stats_long_annot=stats_long_annot[grep("_uc",Sample_Name,invert=TRUE)]


#plot all sample stats
pdf("sampleStats.pdf",width=180,height=48)
ggplot(stats_long_annot,aes(x=species,y=value,col=sub_average))+
  geom_point(alpha=0.5,position=position_jitter(width = .5))+facet_wrap(~variable,scale="free",ncol=1)+scale_color_manual(values=c("FALSE"="black","TRUE"="red"))
dev.off()

#plot all ref stats
pdf("refStats.pdf",width=180,height=25)
ggplot(ref_long,aes(x=species,y=value,fill=sub_average))+geom_bar(stat="identity")+facet_wrap(~variable,scale="free",ncol=1)+scale_fill_manual(values=c("FALSE"="black","TRUE"="red"))
dev.off()

#plot correlation between mapping rate and non motif reads
pdf("mappingRateVSnonMotif.pdf",width=30,height=12)
ggplot(stats_annot,aes(x=mapping_efficiency,y=others,col=species))+
  geom_point(alpha=0.6)+geom_abline(intercept=100,slope=-1)+ylab("% non-motif reads")+
  xlim(c(0,100))+ylim(c(0,100))
dev.off()

#plot specific stats as boxplot
pdf("stats_boxpl.pdf",height=15,width=15)
ggplot(stats_long_annot,aes(x=variable,y=value))+geom_violin()+
  geom_boxplot(outlier.shape=NA,fill="white",col="black")+xlab("")+facet_wrap(~variable,scales="free")
dev.off()

#focus on conversion rate
conversion_problem=stats_annot[conversionRate<98|k1_unmeth>0.02]
stats_annot[conversionRate<98&k1_unmeth>0.02]

pdf("stats_underconv.pdf",height=4,width=7)
ggplot(conversion_problem,aes(x=conversionRate,y=k1_unmeth))+geom_point(alpha=0.5,aes(col=conversionRate>98|k1_unmeth<0.02))+geom_hline(yintercept=0.02,lty=20)+geom_vline(xintercept=98,lty=20)+scale_color_manual(values=c("TRUE"="grey","FALSE"="red"))+ylim(0,1)+xlim(c(0,100))+theme(aspect.ratio=1)
dev.off()

ggplot(conversion_problem,aes(x=conversionRate,y=cont_rat))+geom_point()


#######################################################################################################
###plot subgroups/selected
sub_name="NWU"
sub_stats_long=stats_long_annot[Box=="NWU1"]
sub_ref=ref_long[species%in%sub_stats_long$abbreviation_sp]
sub_stats=stats_annot[Sample_Name%in%sub_stats_long$Sample_Name]

pdf(paste0("sampleStats_",sub_name,".pdf"),width=30,height=48)
ggplot(sub_stats_long,aes(x=species,y=value,col=sub_average))+geom_point(alpha=0.5,position=position_jitter(width = .5))+facet_wrap(~variable,scale="free",ncol=1)+scale_color_manual(values=c("FALSE"="black","TRUE"="red"))
dev.off()

#plot all ref stats
pdf(paste0("refStats_",sub_name,".pdf"),width=30,height=25)
ggplot(sub_ref,aes(x=species,y=value,fill=sub_average))+geom_bar(stat="identity")+facet_wrap(~variable,scale="free",ncol=1)+scale_fill_manual(values=c("FALSE"="black","TRUE"="red"))
dev.off()

#plot correlation between mapping rate and non motif reads
pdf(paste0("mappingRateVSnonMotif_",sub_name,".pdf"),width=9,height=6)
ggplot(sub_stats,aes(x=mapping_efficiency,y=others,col=species))+
  geom_point(alpha=0.6,size=1.2)+geom_abline(intercept=100,slope=-1)+ylab("% non-motif reads")+xlim(c(0,100))+ylim(c(0,100))
dev.off()


