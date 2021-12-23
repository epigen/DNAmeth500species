source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(MASS)

wd=file.path(analysis_dir,"04_tissue")
dir.create(wd)
setwd(wd)

files=system(paste0("ls ",data_dir,"/results_pipeline/*/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/*_mean_meth.tsv"),intern = TRUE)

stats_annot[,mean_qual_tier_relaxed:=mean(c(mapping_efficiency_qual_tier, others_qual_tier,fragments_uncovered_perc_qual_tier,cont_rat_qual_tier)),by=1:nrow(stats_annot)]
stats_annot[,mean_overlap_perc:=(mean_overlap/coveredCpGs)*100,]


stats_annot$class_short <- sapply(as.character(stats_annot$color_class),function(x) class_short[[x]])
stats_annot$class_short <- factor(stats_annot$class_short, levels = class_short)



#-----set filter TRUE/FALSE and filter threshold for mean_overlap_perc
filter=TRUE
thres=50
#---------------------#

remove_samples=stats_annot[conversion_type=="converted"][mean_overlap_perc<thres|species_check=="fail"]$Sample_Name
remove_samples=c(paste0(remove_samples,".meth"),paste0(remove_samples,".cov"))

all_dist=data.table()
all_dist_orig=data.table()
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
  
  if(ncol(mean_meth)<10){next}
  meth_df=as.data.frame(mean_meth)[,seq(7,ncol(mean_meth),by=2)]
  cov_df=as.data.frame(mean_meth)[,seq(8,ncol(mean_meth),by=2)]
  meth_df[cov_df<4]=NA
  
  meth.dist=as.data.table(as.matrix(cor(meth_df,use="pairwise.complete.obs",method="pear")),keep.rownames="sample1")
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
  
  mean_repl=meth.dist_long_red[tis1==tis2,.(mode="replicate",species=spec1[1],N=.N,mean=mean(dist),meanR2=mean(dist^2),median=median(dist),sd=sd(dist)),]
  mean_tis=meth.dist_long_red[repl1==repl2,.(mode="tissue",species=spec1[1],N=.N,mean=mean(dist),meanR2=mean(dist^2),median=median(dist),sd=sd(dist)),]
  
  all_dist_orig=rbindlist(list(all_dist_orig,meth.dist_long_red))
  all_dist=rbindlist(list(all_dist,mean_repl,mean_tis))
}


if(filter==TRUE){
  write.table(all_dist,paste0("all_cor_intfilt_",thres,".tsv"),sep="\t",quote=FALSE,row.names = FALSE)
  write.table(all_dist_orig,paste0("all_cor_orig_intfilt_",thres,".tsv"),sep="\t",quote=FALSE,row.names = FALSE)
}else if (filter==FALSE){
  write.table(all_dist,"all_cor.tsv",sep="\t",quote=FALSE,row.names = FALSE)
}

#-----------------------
#load data if not yet there
#all_dist=fread(paste0("all_cor_intfilt_",thres,".tsv"))
#all_dist_orig=fread(paste0("all_cor_orig_intfilt_",thres,".tsv"))
#-------------------------

#annotate the species
sp_annot=stats_annot[conversion_type=="converted",.(scientific_name=scientific_name[1],color_class=color_class[1],ncbi_order=ncbi_order[1],
                                                    species_check=all(species_check=="pass"),min_overlap=mean(min_overlap_perc),
                                                    mean_overlap=min(mean_overlap_perc),mean_qual_tier=mean(mean_qual_tier,na.rm=TRUE),
                                                    max_qual_tier=max(max_qual_tier),max_prefrag=max(others),mean_prefrag=mean(others),
                                                    mean_PDR=sum(PDR_summ)/sum(N_sites),mean_meth=mean(CpG_meth)),by=species]
all_dist_annot=merge(all_dist,sp_annot,by="species")
all_dist_annot[,rt:=ifelse(.N==2,TRUE,FALSE),by="species"]
                                  
#relationship with mean methylation
pdf("meanR2_mean_meth.pdf",height=2.5,width=7)
ggplot(all_dist_annot[rt==TRUE],aes(x=mean_meth,y=meanR2,col=color_class))+geom_point(alpha=0.6)+geom_smooth(method="lm",se = F)+facet_wrap(~mode)+scale_color_manual(values = class_colors)+ylab("Mean variance explained")+xlab("CpG methylation (%)")
ggplot(all_dist_annot[rt==TRUE],aes(x=mean_meth,y=meanR2,col=color_class))+geom_point(alpha=0.6)+geom_smooth(method="lm",se = F,color="black")+facet_wrap(~mode)+scale_color_manual(values = class_colors)+ylab("Mean variance explained")+xlab("CpG methylation (%)")
dev.off()
                                  
#check relationship between variance explained and prefragmentation/PDR
get_cor=function(x,y){
  if (!length(x)>2){return(list('cor'=0,'p'=1,x=min(x),y=max(y)))}
  r=cor.test(x,y)
  return(list('cor'=r$estimate,'p'=r$p.value,x=min(x),y=max(y)))
}

cors_pdr=all_dist_annot[rt==TRUE&color_class!="Chondrichthyes",get_cor(mean_PDR,meanR2),by=c("color_class","mode")]
cors_pdr[,p_adj:=p.adjust(p,"fdr"),]

cors_frag=all_dist_annot[rt==TRUE&color_class!="Chondrichthyes",get_cor(mean_prefrag,meanR2),by=c("color_class","mode")]
cors_frag[,p_adj:=p.adjust(p,"fdr"),]


pdf("meanR2_prefragmentation.pdf",height=3,width=11)
ggplot(all_dist_annot[rt==TRUE&color_class!="Chondrichthyes"],aes(x=mean_prefrag,y=meanR2,col=color_class))+geom_point(alpha=0.6)+
  geom_smooth(method="lm",se = T,color="black")+facet_grid(mode~color_class)+
  geom_text(data=cors_frag,size=3,vjust=0.5,hjust=0,aes(x=0,y=0.25,label=paste0('p=',signif(p_adj,3),'\nr=',signif(cor,3))))+
  scale_color_manual(values = class_colors)+ylab("Mean variance explained")+xlab("DNA prefragmentation (%)")
dev.off()

pdf("meanR2_pdr.pdf",height=3,width=11)
ggplot(all_dist_annot[rt==TRUE&color_class!="Chondrichthyes"],aes(x=mean_PDR,y=meanR2,col=color_class))+geom_point(alpha=0.6)+
  geom_smooth(method="lm",se = T,color="black")+facet_grid(mode~color_class)+
  geom_text(data=cors_pdr,size=3,vjust=0.5,hjust=0,aes(x=0,y=0.25,label=paste0('p=',signif(p_adj,3),'\nr=',signif(cor,3))))+
  scale_color_manual(values = class_colors)+ylab("Mean variance explained")+xlab("PDR")
dev.off()

#correlate difference in variance explained with PDR/prefrag --> no consistent/significant signal
deltas=all_dist_annot[rt==TRUE&color_class!="Chondrichthyes",.(delta_var=meanR2[mode=="replicate"]-meanR2[mode=="tissue"]),by=c("species","color_class","mean_PDR","mean_prefrag")]

cors_delta_pdr=deltas[,get_cor(delta_var,mean_PDR),by="color_class"]
cors_delta_pdr[,p_adj:=p.adjust(p,"fdr"),]

cors_delta_frag=deltas[,get_cor(delta_var,mean_prefrag),by="color_class"]
cors_delta_frag[,p_adj:=p.adjust(p,"fdr"),]


ggplot(deltas,aes(x=mean_prefrag,y=delta_var,col=color_class))+geom_point(alpha=0.6)+
  geom_smooth(method="lm",se = T,color="black")+facet_grid(~color_class)+
  geom_text(data=cors_delta_frag,size=3,vjust=0.5,hjust=0,aes(x=0,y=-0.1,label=paste0('p=',signif(p_adj,3),'\nr=',signif(cor,3))))+
  scale_color_manual(values = class_colors)+ylab("Delta variance explained")+xlab("DNA prefragmentation (%)")

ggplot(deltas,aes(x=mean_PDR,y=delta_var,col=color_class))+geom_point(alpha=0.6)+
  geom_smooth(method="lm",se = T,color="black")+facet_grid(~color_class)+
  geom_text(data=cors_delta_pdr,size=3,vjust=0.5,hjust=0,aes(x=0,y=-0.1,label=paste0('p=',signif(p_adj,3),'\nr=',signif(cor,3))))+
  scale_color_manual(values = class_colors)+ylab("Delta variance explained")+xlab("PDR")

                                  
#transform to wide format
all_dist_annot_wide=dcast(all_dist_annot,species+color_class~mode,value.var="meanR2")

#remove incomplete species for difference test
all_dist_annot_wide_filt=all_dist_annot_wide[!is.na(replicate)&!is.na(tissue)]

#perform paired t-test on classes with >1 species, to test for difference between 
all_dist_annot_wide_filt[,N:=.N,by="color_class"]
#p.values=all_dist_annot_wide_filt[N>1,signif(t.test(replicate,tissue,paired=TRUE)$p.value,digits = 3),by="color_class"]
p.values=all_dist_annot_wide_filt[N>1,signif(wilcox.test(replicate,tissue,paired=TRUE)$p.value,digits = 3),by="color_class"]


pdf("tissues_cor.pdf",height=6,width=12)
ggplot(all_dist_annot_wide,aes(x=replicate,y=tissue,fill=color_class))+geom_text(data = p.values,aes(x=0.1,y=0.9,label=paste0("p=",V1)),hjust=0)+geom_point(pch=21,size=2,col="black",alpha=0.6)+geom_abline(slope=1,size=1,lty=2,intercept=c(0,0))+scale_fill_manual(values = class_colors)+coord_fixed()+xlim(c(0,1))+ylim(c(0,1))+facet_wrap(~color_class,ncol=4)+xlab("Variance explained by tissue")+ylab("Variance explained by individual")
dev.off()

#make wordcloud of tissues used in each class
library(wordcloud)

used_samples=all_dist_orig[(tis1==tis2|repl1==repl2)&spec1%in%all_dist_annot_wide_filt$species,gsub(".meth","",unique(c(sample1,sample2))),]
tissues=stats_annot[Sample_Name%in%used_samples,.N,by=c("color_class","Tissue")]

#calculate relative tissue frequencies per class
tissues[,rel_freq:=N/sum(N),by="color_class"]

#order by class
tissues[,caption:=factor(color_class,levels=classes),]
tissues=tissues[order(caption)]


wordcoud_title = function(x){
  wordcloud(x$Tissue,x$rel_freq,colors=brewer.pal(8, "Dark2"))
  title(unique(x$caption))
}


pdf("tissues_wordcloud.pdf",height=7,width=15)
layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow = TRUE))
margin(0,0,0,0)
tissues[,wordcoud_title(.SD),by="color_class"]
dev.off()


#check "genetic similarity" of used samples 
#--> no striking difference between groups that show more tissue specific meth and those that show less
stats_annot_used=stats_annot[Sample_Name%in%used_samples,]
pdf("tissues_CoG_overlap.pdf",height=2.5,width=4)
ggplot(stats_annot_used,aes(y=mean_overlap_perc,x=class_short,fill=color_class))+geom_boxplot()+scale_fill_manual(values = class_colors)+ylab("Mean CpG overlap (%)")+rotate_labels(angle = 45,vjust = 1)
dev.off()
                                  
# plot correlations of all the sample pairs per species grouped by class as qc
all_dist_orig=fread(paste0("all_cor_orig_intfilt_",thres,".tsv"))

#annotate the species
sp_annot=stats_annot[conversion_type=="converted",.(scientific_name=scientific_name[1],color_class=color_class[1],ncbi_order=ncbi_order[1],species_check=all(species_check=="pass")),by=species]
all_dist_orig_annot=merge(all_dist_orig,sp_annot,by.x="spec1",by.y="species")

pdf("methylation_correlation.pdf",height=4,width=6)
ggplot(all_dist_orig_annot,aes(x=color_class,y=dist))+geom_boxplot(aes(fill=color_class))+scale_fill_manual(values = class_colors)+ylab("Pairwise pearson correlation")+rotate_labels(angle = 45,vjust = 1)+ylim(c(0,1))+geom_hline(yintercept = 0.75,lty=2)+stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=2)
dev.off()

                                  
