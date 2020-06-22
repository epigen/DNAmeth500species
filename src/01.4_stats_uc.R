source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"01_basicStats")
setwd(wd)


uc_samples=stats_annot[conversion_type=="unconverted"]$Sample_Name_unif

stats_uc=stats_annot[Sample_Name_unif%in%uc_samples,]
stats_uc=stats_uc[order(Sample_Name)]

pdf("uc_coverage_boxpl.pdf",height=3,width=3)
ggplot(stats_uc,aes(x=conversion_type,y=fragments_uncovered_perc))+geom_point(col="grey",alpha=0.3,position=position_jitter(width=0.3))+geom_boxplot(outlier.shape=NA,fill="transparent",col="black")+xlab("")+ylab("Uncovered fragments (%)")
dev.off()

pdf("uc_mapeff_boxpl.pdf",height=3,width=3)
ggplot(stats_uc,aes(x=conversion_type,y=mapping_efficiency))+geom_point(col="grey",alpha=0.3,position=position_jitter(width=0.3))+geom_boxplot(outlier.shape=NA,fill="transparent",col="black")+xlab("")+ylab("Mapping efficiency (%)")
dev.off()

pdf("uc_lowest_common_boxpl.pdf",height=3,width=3)
ggplot(stats_uc,aes(x=conversion_type,y=max_lowest_common))+geom_point(col="grey",alpha=0.3,position=position_jitter(width=0.3))+geom_boxplot(outlier.shape=NA,fill="transparent",col="black")+geom_hline(yintercept=9,lty=20,col="red")+xlab("")+ylab("Lowest common rank")
dev.off()


stats_uc_calc=stats_uc[,list(total_reads=total_reads[2]/total_reads[1],mapping_efficiency=mapping_efficiency[2]/mapping_efficiency[1],coveredCpGs=coveredCpGs[2]/coveredCpGs[1],CTGG=(CGG[2]+TGG[2])/(CGG[1]+TGG[1]),fragments_uncovered_perc=fragments_uncovered_perc[2]/fragments_uncovered_perc[1],lowest_common_rank=max_lowest_common[2]/max_lowest_common[1]),by="Sample_Name_unif"]

stats_uc_calc_long=melt(stats_uc_calc,id.vars="Sample_Name_unif")
stats_uc_calc_long[,log2_value:=log2(value)]
stats_uc_calc_long[,species:=sub("_.*","",Sample_Name_unif),]

stats_uc_calc_long[,mean_value:=mean(log2_value,na.rm=TRUE),by="variable"]
stats_uc_calc_long[,sd_value:=sd(log2_value,na.rm=TRUE),by="variable"]
stats_uc_calc_long[,abs_diff_mean:=abs(log2_value-mean_value),]
stats_uc_calc_long[,mark:=ifelse(abs_diff_mean>sd_value,"diff","ok"),]

stats_uc_calc_long[,mark_count:=sum(mark=="diff"),by="Sample_Name_unif"]
stats_uc_calc_long[,species:=factor(species,levels=unique(species[order(log2_value[variable=="fragments_uncovered_perc"])]))]

pdf("uc_stats.pdf",height=10,width=90)
ggplot(stats_uc_calc_long,aes(x=species,y=log2_value,col=mark))+geom_point()+facet_wrap(~variable,ncol=1,scale="free_y")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+ylab("log2(unconv/conv)")
dev.off()

pdf("uc_stats_sel.pdf",height=10,width=25)
ggplot(stats_uc_calc_long[mark_count>=2],aes(x=species,y=log2_value,col=mark))+geom_point()+facet_wrap(~variable,ncol=1,scale="free_y")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+ylab("log2(unconv/conv)")
dev.off()

pdf("uc_stats_boxpl.pdf",height=6,width=8)
ggplot(stats_uc_calc_long,aes(x=variable,y=log2_value))+geom_point(aes(col=abs(log2_value)>1),alpha=0.3,position=position_jitter(width=0.3))+geom_boxplot(outlier.shape=NA,fill="transparent",col="black")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+ylab("log2(unconv/conv)")+xlab("")+scale_color_manual(values=c("TRUE"="red","FALSE"="grey"))
dev.off()


#check most important parameters to assess quality of unconverted
#too few reads (probably pooling problem) SRP_1_H_uc rerun failed
summary(stats_uc[grepl("_uc",Sample_Name)]$total_reads)
low_reads_samples=stats_uc[Sample_Name_unif%in%Sample_Name_unif[total_reads<8000000&grepl("_uc",Sample_Name)]][order(total_reads)][order(Sample_Name)]
write.table(low_reads_samples,"low_reads_samples.tsv",quote=FALSE,sep="\t",row.names=FALSE)


#too low mapping rate (might indicate sample mixup) --> conclusion: all seem ok (rel. similar mapping rates und uncovererd in uc and not uc) except for CH_2_LU_uc --> investigate with Amelie
low_map_samples=stats_uc[Sample_Name_unif%in%Sample_Name_unif[mapping_efficiency<40&others<60]][,c("Sample_Name","species","total_reads","mapping_efficiency","fragments_uncovered_perc"),with=FALSE]
write.table(low_map_samples,"low_map_samples.tsv",quote=FALSE,sep="\t",row.names=FALSE)

