source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

library(MASS)

wd=file.path(analysis_dir,"08_tissue")
dir.create(wd)
setwd(wd)

files=system(paste0("ls ",data_dir,"/results_pipeline/*/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/*_mean_meth.tsv"),intern = TRUE)
sel_species=c("NSS_","WIS_","NOP_","AX_","HE_","EAB_","NY_","KAN_")

unique(stats_annot[grepl(paste0(sel_species,collapse = "|"),Sample_Name),c("scientific_name","English","abbreviation_sp")])
unique(stats_annot[grepl(paste0(sel_species,collapse = "|"),Sample_Name),c("Tissue")])

files=grep(paste0(sel_species,collapse = "|"),files,value = TRUE)
all_mds=data.table()
for (file in files){
  mean_meth=fread(file)
  #MDS on fragments
  meth.dist=dist(t(as.data.frame(mean_meth)[,seq(7,ncol(mean_meth),by=2)]),method="eucl")

  meth.dist[meth.dist<=0]=0.00001
  meth.mds=isoMDS(meth.dist)
  MDS=data.table(Sample_Name=as.character(lapply(row.names(meth.mds$points),function(x){unlist(strsplit(x,"\\."))[1]})),MDS1=meth.mds$points[,1],MDS2=meth.mds$points[,2])

  all_mds=rbindlist(list(all_mds,MDS))
}

all_mds[,species:=unlist(lapply(strsplit(Sample_Name,"_"),"[[",1)),]
all_mds[,repl:=unlist(lapply(strsplit(Sample_Name,"_"),"[[",2)),]
all_mds[,tissue:=unlist(lapply(strsplit(Sample_Name,"_"),"[[",3)),]

all_mds[,tissue_simpl:=ifelse(tissue%in%names(tissue_colors),tissue,"X"),]


all_mds=merge(all_mds, stats_annot[,c("Enrichment cycles","Sample_Name"),],by="Sample_Name")

all_mds[,species:=factor(species,levels=c("NSS","WIS","NOP","AX","HE","EAB","KAN","NY")),]

pdf("tissueMDS_repl.pdf",height=6,width=4.5)
ggplot(all_mds,aes(x=MDS1,y=MDS2,col=tissue_simpl),)+geom_text(aes(label=repl),fontface = "bold")+facet_wrap(~species,ncol=2,scale="free")+scale_color_manual(values = tissue_colors)+theme(axis.text=element_blank(), axis.ticks=element_blank())
dev.off()

pdf("tissueMDS_shape.pdf",height=6,width=4.5)
ggplot(all_mds,aes(x=MDS1,y=MDS2,col=tissue_simpl,shape=repl),)+geom_point()+facet_wrap(~species,ncol=2,scale="free")+scale_color_manual(values = tissue_colors)+theme(axis.text=element_blank(), axis.ticks=element_blank())
dev.off()


pdf("tissueMDS_cycles.pdf",height=3.5,width=9)
ggplot(all_mds,aes(x=MDS1,y=MDS2,col=tissue_simpl))+geom_text(aes(label=`Enrichment cycles`),fontface = "bold")+facet_wrap(~species,ncol=4,scale="free")+
  scale_color_manual(values = tissue_colors)+theme(axis.text=element_blank(), axis.ticks=element_blank())
dev.off()

