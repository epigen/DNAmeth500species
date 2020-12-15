source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(MASS)

wd=file.path(analysis_dir,"04_tissue")
dir.create(wd)
setwd(wd)

files=system(paste0("ls ",data_dir,"/results_pipeline/*/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/*_mean_meth.tsv"),intern = TRUE)

stats_annot[,mean_qual_tier_relaxed:=mean(c(mapping_efficiency_qual_tier, others_qual_tier,fragments_uncovered_perc_qual_tier,cont_rat_qual_tier)),by=1:nrow(stats_annot)]
stats_annot[,mean_overlap_perc:=(mean_overlap/coveredCpGs)*100,]


#-----set filter TRUE/FALSE and filter threshold for mean_overlap_perc
filter=FALSE
thres=50
#---------------------#


remove_samples=stats_annot[conversion_type=="converted"][mean_overlap_perc<thres|species_check=="fail"]$Sample_Name
remove_samples=c(paste0(remove_samples,".meth"),paste0(remove_samples,".cov"))


all_dist=data.table()
for (file in files){
  print(file)
  mean_meth=fread(file)
  
  if(filter==TRUE){
    remove_samples_sel=remove_samples[remove_samples%in%colnames(mean_meth)]
    print(paste0("Removing ",length(remove_samples_sel)," samples."))
    if(length(remove_samples_sel)==(ncol(mean_meth)-6)){
      next
    }
    if (length(remove_samples_sel)>0){
    mean_meth=mean_meth[,(remove_samples_sel):=NULL,]
    }
  }
  

  meth.dist=as.data.table(as.matrix(dist(t(as.data.frame(mean_meth)[,seq(7,ncol(mean_meth),by=2)]),method="eucl")),keep.rownames="sample1")
  meth.dist_long=melt(meth.dist,id.vars = "sample1",value.name = "dist",variable.name = "sample2")
  meth.dist_long[,sample2:=as.character(sample2),]
  
  meth.dist_long[,spec1:=unlist(strsplit(sample1,"_|\\."))[1],by=1:nrow(meth.dist_long)]
  meth.dist_long[,repl1:=unlist(strsplit(sample1,"_|\\."))[2],by=1:nrow(meth.dist_long)]
  meth.dist_long[,tis1:=unlist(strsplit(sample1,"_|\\."))[3],by=1:nrow(meth.dist_long)]
  meth.dist_long[,spec2:=unlist(strsplit(sample2,"_|\\."))[1],by=1:nrow(meth.dist_long)]
  meth.dist_long[,repl2:=unlist(strsplit(sample2,"_|\\."))[2],by=1:nrow(meth.dist_long)]
  meth.dist_long[,tis2:=unlist(strsplit(sample2,"_|\\."))[3],by=1:nrow(meth.dist_long)]
  
  meth.dist_long_red=meth.dist_long[tis1!=tis2|repl1!=repl2]
  if(nrow(meth.dist_long_red)==0){next}
  meth.dist_long_red[,combi:=paste0(sort(c(sample1,sample2)),collapse="__"),by=1:nrow(meth.dist_long_red)]
  meth.dist_long_red=meth.dist_long_red[!duplicated(combi)]
  
  mean_repl=meth.dist_long_red[tis1==tis2,.(mode="replicate",species=spec1[1],N=.N,mean=mean(dist),median=median(dist),sd=sd(dist)),]
  mean_tis=meth.dist_long_red[repl1==repl2,.(mode="tissue",species=spec1[1],N=.N,mean=mean(dist),median=median(dist),sd=sd(dist)),]
  
  all_dist=rbindlist(list(all_dist,mean_repl,mean_tis))
}

if(filter==TRUE){
  write.table(all_dist,paste0("all_dist_intfilt_",thres,".tsv"),sep="\t",quote=FALSE,row.names = FALSE)
}else if (filter==FALSE){
  write.table(all_dist,"all_dist.tsv",sep="\t",quote=FALSE,row.names = FALSE)
}


all_dist_rat=all_dist[,.(ratio_N=N[mode=="tissue"]/N[mode=="replicate"],N_tissue=N[mode=="tissue"],N_replicate=N[mode=="replicate"],ratio_sd=sd[mode=="tissue"]/sd[mode=="replicate"],ratio_mean=mean[mode=="tissue"]/mean[mode=="replicate"]),by=species]

sp_annot=stats_annot[conversion_type=="converted",.(scientific_name=scientific_name[1],color_class=color_class[1],ncbi_order=ncbi_order[1],species_check=all(species_check=="pass"),min_overlap=mean(min_overlap_perc),mean_overlap=min(mean_overlap_perc),mean_qual_tier=mean(mean_qual_tier,na.rm=TRUE),max_qual_tier=max(max_qual_tier)),by=species]

all_dist_rat_annot=merge(all_dist_rat,sp_annot,by="species")

all_dist_rat_annot$class_short <- sapply(as.character(all_dist_rat_annot$color_class),
                                         function(x) class_short[[x]])
all_dist_rat_annot$class_short <- factor(all_dist_rat_annot$class_short, levels = class_short)


if (filter==TRUE){
pdf(paste0("dist_ratio_class_boxpl_intfilt_",thres,"_short_names.pdf"),height=3.5,width=4.5)
pl=ggplot(all_dist_rat_annot,aes(y=log2(ratio_mean),x=class_short,col=color_class))+geom_boxplot()+scale_color_manual(values = class_colors)+geom_hline(yintercept = 0,lty=21)+ylab("log2(tissue dist./replicate dist.")+xlab("")
pl2=ggplot(all_dist_rat_annot[N_replicate>1],aes(y=log2(ratio_mean),x=class_short,col=color_class))+geom_boxplot()+scale_color_manual(values = class_colors)+geom_hline(yintercept = 0,lty=21)+ylab("log2(tissue dist./replicate dist.")+xlab("")
print(pl)
print(pl2)
dev.off()

}else if (filter==FALSE){

pdf("dist_ratio_class_boxpl_short_names.pdf",height=3.5,width=4.5)
pl=ggplot(all_dist_rat_annot,aes(y=log2(ratio_mean),x=class_short,col=color_class))+geom_boxplot()+scale_color_manual(values = class_colors)+geom_hline(yintercept = 0,lty=21)+rotate_labels()+ylab("log2(tissue dist./replicate dist.")+xlab("")

pl2=ggplot(all_dist_rat_annot[N_replicate>1],aes(y=log2(ratio_mean),x=class_short,col=color_class))+geom_boxplot()+scale_color_manual(values = class_colors)+geom_hline(yintercept = 0,lty=21)+ylab("log2(tissue dist./replicate dist.")+xlab("")
print(pl)
print(pl2)
dev.off()

pdf("dist_ratio_class_boxpl_filt.pdf",height=3.5,width=4.5)
pl2=ggplot(all_dist_rat_annot[species_check==TRUE&mean_qual_tier<2.5],aes(y=log2(ratio_mean),x=color_class,col=color_class))+geom_boxplot()+scale_color_manual(values = class_colors)+geom_hline(yintercept = 0,lty=21)+rotate_labels()+ylab("log2(tissue dist./replicate dist.")+xlab("")
print(pl2)
dev.off()

}
