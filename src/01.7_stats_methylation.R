source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"01_basicStats/01.7_stats_methylation")
dir.create(wd)
setwd(wd)


files=system(paste0("ls ",data_dir,"/results_pipeline/*/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/*_mean_meth.tsv"),intern = TRUE)


get_histo=function(m,c){
  h=hist(m,breaks=c(0:100),plot=FALSE)
  res=data.table(bin=1:100,counts=h$counts,frequency=h$density)
  res[,cov_thres:=c,]
  return(res)
}


cov_thresholds=c(0,1,2,3,4,5,10,15,20,30)

all_means=data.table()
all_histo=data.table()

for(file in files){
  print(file)
  meth_data_mean=fread(file)
  
  meth_data_mean_long=melt(meth_data_mean,id.vars=c("meta","score","consN","mean_diffNts","max_diffNts",
                                "min_diffNts"),measure=patterns(".cov", ".meth"),variable.factor=FALSE,
                                  variable.name="sample_index",value.name=c("cov","meth"),na.rm=TRUE)
  
  sample_names=data.table(sample=unique(sub(".meth","",grep(".meth",names(meth_data_mean),value=TRUE))))
  sample_names[,sample_index:=as.character(1:nrow(sample_names)),]
  
  meth_data_mean_long=merge(meth_data_mean_long,sample_names,by="sample_index")
  
  for (cov_thres in cov_thresholds){
    means=meth_data_mean_long[cov>=cov_thres,.(mean=mean(meth),sd=sd(meth),cov_thres=cov_thres),by=sample]
    histo=meth_data_mean_long[cov>=cov_thres,get_histo(meth,cov_thres),by=sample]
    
    all_means=rbindlist(list(all_means,means))
    all_histo=rbindlist(list(all_histo,histo))
    }
}

write.table(all_means,"all_meth_means.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(all_histo,"all_meth_histo.tsv",sep="\t",quote=FALSE,row.names=FALSE)


all_means_annot=merge(all_means,stats_annot,by.x="sample",by.y="Sample_Name")
all_histo_annot=merge(all_histo,stats_annot,by.x="sample",by.y="Sample_Name")


pdf("meth_means_boxpl.pdf",height=5,width=15)
ggplot(all_means_annot,aes(y=mean,x=color_class,fill=color_class))+geom_boxplot()+rotate_labels()+
facet_wrap(~cov_thres,ncol=5,scale="free")+ylab("% DNA methylation (CpG)")+
#stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)+
scale_fill_manual(values = class_colors)+scale_x_discrete(labels=class_short)
dev.off()



pdf("meth_hist_boxpl.pdf",height=5,width=20)
ggplot(all_histo_annot,aes(y=frequency,x=bin,col=color_class))+stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),position=position_dodge(width=5), geom="errorbar", width=0.1)+stat_summary(fun.y=mean, geom="point",position=position_dodge(width=5),size=1)+scale_color_manual(values = class_colors)+facet_wrap(~cov_thres,ncol=5,scale="free")
dev.off()


all_histo_annot[,large_bins:=factor(ifelse(bin<21,"low",ifelse(bin>79,"high","medium")),levels=c("low","medium","high")),]
all_histo_annot_red=all_histo_annot[,.(frequency=sum(frequency)),by=c("large_bins","color_class","sample","species","cov_thres")]


pdf("meth_hist_boxpl_grouped.pdf",height=5,width=20)
ggplot(all_histo_annot_red,aes(y=frequency,x=large_bins,col=color_class))+geom_boxplot()+
scale_color_manual(values = class_colors)+facet_wrap(~cov_thres,ncol=5,scale="free")+rotate_labels()
dev.off()
