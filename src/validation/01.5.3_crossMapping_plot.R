source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(uwot)
library(igraph)
library(biomaRt)
library(rpart)
library("ROCR")
library("caret")

wd=file.path(analysis_dir,"SET_PROPER_PATH/01.5_crossMapping")
dir.create(wd,recursive=TRUE)
setwd(wd)


genomes=fread("species_ucsc_matches.tsv")
genomes[ucsc_db=="canFam4",ucsc_db:="canFam5"]
profiles=stats=fread("Methylation_profiles_merged.tsv.gz")
stats=fread("Methylation_profile_stats_summary.tsv")
ref=fread("tracks_present.tsv")
ref=unique(ref)

#make sure only one genome per species is used
profiles[,sel_genome:=sort(mapped_genome)[1],by=c("sample")]
profiles=profiles[sel_genome==mapped_genome]
stats[,sel_genome:=sort(mapped_genome)[1],by=c("label")]
stats=stats[sel_genome==mapped_genome]


profiles[,N:=.N,by="sample"]
Ns=unique(profiles[,c("sample","N")])
profiles[,species:=sub("_.*","",sample),by="sample"]
profiles[,tissue:=sub(".*_","",sample),by="sample"]
profiles[,individual:=as.numeric(sub("_[A-Za-z]*","",sub("[A-Z]*_","",sample))),by="sample"]
profiles[,N_spec:=.N,by="species"]
profiles[species=="JL",color_class:="Invertebrata",]
profiles[species=="SQM",color_class:="Mammalia",]


########profile analysis#######################
transcript_factor=200
flank_factor=100

mean_profiles=profiles[N_spec>4000,.(meth=mean(meth)),by=c("species","color_class","trans_type","pos","mapped_genome")]
mean_profiles[,color_class:=factor(color_class,levels=names(class_colors)),]
mean_profiles[,color_class2:=ifelse(species=="JL","JL","other"),]


pdf("profiles.pdf",h=3,w=8)
ggplot(mean_profiles,aes(x=pos,y=meth,color=color_class2))+
  stat_smooth(geom="line",aes(group=species),method='loess',span=0.3,se=FALSE,alpha=0.2)+
  stat_smooth(geom="line",aes(group=color_class2),method='loess',span=0.3,se=FALSE,alpha=1,lwd=1,color="black")+
  geom_vline(xintercept = c(-flank_factor,0,transcript_factor,flank_factor+transcript_factor),lty=20)+
  xlab('')+ylab('DNA methylation at CpGs (%)')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_wrap(~color_class,ncol=4)
dev.off()

#Special examples 
#Chondrichthyes: ESH (Elephant shark) WIS (winter skate)
#Invertebrata: GSU/PSU (Green/Purple Sea Urchin)
#Amphibia: ABU (Asian Bullfrog)

sub_profiles=profiles[species%in%c("ESH","WIS","GSU","ABU")]
sub_profiles=merge(sub_profiles,stats[,c("label","perc_mapped_passed")],by.x="sample", by.y="label")

pdf("profiles_samples.pdf",h=3,w=4)
ggplot(sub_profiles,aes(x=pos,y=meth,color=species))+
  stat_smooth(geom="line",aes(group=sample),method='loess',span=0.3,se=FALSE,alpha=1,lwd=0.5)+
  geom_vline(xintercept = c(-flank_factor,0,transcript_factor,flank_factor+transcript_factor),lty=20,lwd=0.5)+
  xlab('')+ylab('DNA methylation at CpGs (%)')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_wrap(~color_class,ncol=2)
dev.off()

#####Axolotl analysis##################
ax_profile=fread("AX_methylation_profile.tsv")
ax_profile[,N_samp:=.N,by="sample"]
ax_profile_mean=ax_profile[,.(meth=mean(meth)),by=c("sample","pos","trans_type","N_samp")]

pdf("profiles_ax2.pdf",h=3,w=4.5)
ggplot(ax_profile,aes(x=pos,y=meth))+
  stat_smooth(geom="line",aes(group=sample),color="green4",method='loess',span=0.03,se=FALSE,alpha=0.5)+
  stat_smooth(geom="line",method='loess',span=0.03,se=FALSE,alpha=1,lwd=1,color="black")+
  geom_vline(xintercept = c(-flank_factor,0,transcript_factor,flank_factor+transcript_factor),lty=20)+
  xlab('')+ylab('DNA methylation at CpGs (%)')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylim(0,100)
dev.off()


#####mapping stats analysis##################
ref[refGene==1|refSeqComposite==1,transc:="RefSeq"]
ref[xenoRefGene==1,transc:="XeneoRefSeq"]

stats=merge(stats,ref[,c("ucsc_db","transc")],by.x="mapped_genome",by.y="ucsc_db")
stats=merge(stats,Ns, by.x="label",by.y="sample")
stats[,mapping:=ifelse(perc_mapped_passed<20,"low",ifelse(perc_mapped_passed<70,"mid","high")),]
stats[is.na(color_class),color_class:="Invertebrata",]
stats[,color_class:=factor(color_class,levels=names(class_colors)),]
stats[,mapped_genome:=factor(mapped_genome,levels=unique(mapped_genome[order(color_class)])),]
stats[ymin_Upstream<0,ymin_Upstream:=0,] # fix smoothing artefacts (methylation values can't be <0)


stats_red=stats[,.(y95_Transcript=mean(y95_Transcript[N>1000]),ymin_Upstream=mean(ymin_Upstream[N>1000]) ,Nspec=sum(N[N>1000])),by=c("color_class","mapped_genome","perc_mapped_passed","species","transc")]
stats_red[,median_mapped:=median(perc_mapped_passed),by="mapped_genome"]

#annotate with evolutionary distance
stats_red_annot=merge(stats_red,genomes[,c("species","ncbi_name","ucsc_species","ucsc_ncbi_name","common_level","ucsc_db")],
                      by.x="mapped_genome",by.y="ucsc_db",allow.cartesian=TRUE,suffixes=c("",".g"))
stats_red_annot=stats_red_annot[species==species.g]
#simplify levels
stats_red_annot[,common_level:=gsub("below_|sub|parv|infra|super","",gsub("[0-9]","",common_level))]
stats_red_annot[common_level=="cohort",common_level:="class"]
stats_red_annot[common_level=="tribe",common_level:="family"]
table(stats_red_annot$common_level)

stats_red_annot[,common_level:=factor(common_level,levels=c("class","order","family","genus","species")),]
stats_red_annot[,annot_genome:=paste0(ucsc_species," (",mapped_genome,")"),]
stats_red_annot[,annot_genome:=factor(annot_genome,levels=unique(annot_genome[order(median_mapped)])),]


pdf("mapping_evol_dist.pdf",h=5,w=4.4)
ggplot(stats_red_annot,aes(x=common_level,y=perc_mapped_passed,color=color_class))+
  geom_text(aes(label=species),position=position_jitter(width=0.2),size=1.5)+
  geom_boxplot(color="black",outlier.color = NA,fill="transparent")+
  rotate_labels()+
  stat_summary(fun.data = give.n,fun.args = c(y=-4), geom = "text",size=3,color="black")+
  scale_color_manual(values=class_colors)+xlab("")+ylab("Mapped consensus fragments (%)")
dev.off()


pdf("mapping_genome_spec.pdf",h=7,w=12)
ggplot(stats_red_annot, aes(x=annot_genome,y=perc_mapped_passed))+
  #geom_point(position=position_jitter(0.3),alpha=0.5,aes(color=color_class))+
  geom_text(aes(label=species,color=color_class),size=2)+
  geom_boxplot(width=0.6,outlier.color = NA,fill="transparent",lwd=0.4)+
  scale_color_manual(values=class_colors)+
  rotate_labels()+ylab("Mapped consensus fragments (%)")

ggplot(stats_red_annot, aes(x=annot_genome,y=perc_mapped_passed))+
  geom_point(position=position_jitter(0.3),alpha=0.5,aes(color=color_class))+
  geom_boxplot(width=0.6,outlier.color = NA,fill="transparent")+
  scale_color_manual(values=class_colors)+
  rotate_labels()+ylab("Mapped consensus fragments (%)")
dev.off()


pdf("mapping_vs_dip.pdf",h=3,w=5)
ggplot(stats_red[Nspec>0], aes(y=log(y95_Transcript/ymin_Upstream),x=perc_mapped_passed))+
  geom_point(alpha=0.7,aes(color=color_class))+
  geom_smooth()+
  scale_color_manual(values=class_colors)+
  xlab("Mapped consensus fragments (%)")+ylab("Dip-score")
dev.off()



######comparative promoter methylation analysis in human ortholog gene space#######################
#prepare orthologXsample matrix
profiles_prom=profiles[trans_type=="upstream"&pos>-50,.(meth=mean(meth)),by=c("trans_id","sample","mapped_genome","color_class","species","tissue","individual")]

#try gene body methylation --> works really well for class and tissue prediction (all of the below can be run for gene body instead of promoter) but too much for now
#profiles_prom=profiles[trans_type=="transcript",.(meth=mean(meth)),by=c("trans_id","sample","mapped_genome","color_class","species","tissue","individual")]

profiles_prom[,rank_high:=rank(-meth),by="sample"]
profiles_prom[,trans_id_simpl:=sub("\\.[0-9]+","",trans_id)]

##match geneids to gene2refseq db using command line grep
##make sure refsqs.txt has unix end of line characters
#write.table(unique(profiles_prom$trans_id_simpl),"refesqs.txt",row.names = FALSE,col.names = FALSE, quote=FALSE)
##get gene2refseq here:#get here: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
##run in bash (takes 5 sec per loop)
##for i in {1..7}
##do
##time head -n $((100000*$i)) refesqs.txt|tail -n 100000 |grep -f - gene2refseq|awk '{print $1"\t"$2"\t"$4"\t"$14"\t"$19"\t"$13"\t"$14}'>>mappings.tsv
##done
##sort -u mappings.tsv > mappings_uniq.tsv


orthologs=fread("gene_orthologs.gz") #get here: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
orthologs=orthologs[orthologs$`#tax_id`=="9606"] #only keep human
mapping=fread("mappings_uniq.tsv",col.names = c("taxid","geneid","refseq","genome","gene_name","genomic_start","genomic_end"))
mapping[,refseq_simpl:=sub("\\.[0-9]+","",refseq)]
mapping=mapping[!duplicated(refseq)] # remove duplications due to different genome versions
mapping[,genomic_start:=as.numeric(ifelse(genomic_start=="-",0,genomic_start)),] #calculate lengths of transcripts
mapping[,genomic_end:=as.numeric(ifelse(genomic_end=="-",0,genomic_end)),]
mapping[,genomic_length:=genomic_end-genomic_start,]

#merge profiles,stats per species, transcript ids, and orthologs together --> where possible transcripts across species get a human ortholog assigned
profiles_prom=merge(profiles_prom,mapping,by.x="trans_id_simpl",by.y="refseq_simpl")
profiles_prom=merge(profiles_prom,orthologs,by.x="geneid",by.y="Other_GeneID")
profiles_prom=merge(profiles_prom,stats[,c("label","perc_mapped_passed","y95_Downstream", "y95_Transcript", "y95_Upstream", "ymin_Upstream", "xmin_Upstream")],by.x="sample",by.y="label")

#per ortholog and sample only keep the longest transcript
profiles_prom[,keep:=trans_id[which.max(genomic_length)],by=c("sample","GeneID")]
profiles_prom=profiles_prom[keep==trans_id]

#add some stats per ortholog
profiles_prom[,c("otho_sample","ortho_species"):=list(length(unique(sample)),length(unique(species))),by="GeneID"]

#add dip score for each sample
profiles_prom[,dip:=log(y95_Transcript/ymin_Upstream)]

#add gene names for human orthologs
mart<- biomaRt::useMart(biomart = 'ensembl',dataset = "hsapiens_gene_ensembl")
names <- as.data.table(biomaRt::getBM(attributes = c("entrezgene_id","hgnc_symbol"), 
                                 filters=c("entrezgene_id"),
                                 values=unique(profiles_prom$GeneID),
                                 mart=mart))
names=names[hgnc_symbol!=""]
names=names[!duplicated(entrezgene_id)]
profiles_prom=merge(profiles_prom,names,by.x = "GeneID",by.y="entrezgene_id",all.x=TRUE)

# save orthologXsample matrix
write.table(profiles_prom,"promoter_meth_orthologs.tsv",sep="\t",quote=F,row.names=F)#for external use
save(profiles_prom,file = "promoter_meth_orthologs.RData")
# load orthologXsample matrix to start actual analysis
load("promoter_meth_orthologs.RData")

##### start analysis in ortholog space
# stats per ortholog
profiles_prom_red=unique(profiles_prom[,c("GeneID","species","mapped_genome","color_class")])
profiles_prom_red[,ortho_species:=.N,by="GeneID"]
profiles_prom_red[,GeneID:=factor(as.character(GeneID),levels=unique(as.character(GeneID)[order(ortho_species,decreasing=TRUE)])),]
profiles_prom_red[,color_class:=factor(color_class,levels=names(class_colors)),]

pdf("ortholog_representation.pdf", 6,3)
ggplot(profiles_prom_red[GeneID%in%sample(unique(GeneID),size = 500),],aes(x=GeneID,fill=color_class))+
  geom_bar()+scale_fill_manual(values=class_colors)
dev.off()

#make matrix sampleXOrtholog matrix
mat=dcast(profiles_prom,sample+color_class+tissue+species+mapped_genome+dip+perc_mapped_passed~hgnc_symbol+GeneID,value.var = "meth")

#prep annotations
sample_annot=data.frame(mat[,1:7],row.names = "sample")
gene_annot=profiles_prom[,.(N_samples=length(unique(sample)),N_species=length(unique(species)),N_classes=length(unique(color_class))),by=c('hgnc_symbol','GeneID')]
gene_annot[,gene_id:=paste0(hgnc_symbol,'_',GeneID),]
gene_annot=data.frame(gene_annot[order(match(gene_id, colnames(mat[,-c(1:7)])))],row.names='gene_id')
stopifnot(gene_annot$gene_id==colnames(mat[,-c(1:7)]))

#prep numerical matrix
mat_num=as.matrix(mat[,-c(1:7)])
row.names(mat_num)=mat$sample

cs=colSums(!is.na(mat_num))
rs=rowSums(!is.na(mat_num))


####sample-wise clustering
mat_red=mat_num[rs>400,cs>100]
gene_annot_red=gene_annot[colnames(mat_red),]
sample_annot_red=sample_annot[rownames(mat_red),]

dim(mat_red)
#[1]  1524 14339
sum(is.na(mat_red))/prod(dim(mat_red))
#[1] 0.8140666

as.data.table(sample_annot_red)[,.(species=length(unique(species)),samples=.N),]
#species samples
#    382    1524


#scrambled (mix up methylation values, but retain NA structure)
#set.seed(42)
#mat_red[!is.na(mat_red)]=sample(mat_red[!is.na(mat_red)],sum(!is.na(mat_red)))

cor=cor(t(mat_red),use = "pairwise.complete.obs")
row.names(cor)=row.names(mat_red)
colnames(cor)=row.names(mat_red)

sum(is.na(cor))


#plot umap for all
set.seed(43)
#umap <- umap(cor, n_neighbors = 10, learning_rate = 0.5, init = "random", n_epochs = 200,min_dist=0.5,spread=1.5) #makes very tight clusters --> bad readability
umap <- umap(cor, n_neighbors = 20, learning_rate = 0.5, init = "random", n_epochs = 200,min_dist=2,spread=1)

umap_dt=as.data.table(sample_annot_red,keep.rownames = TRUE)
umap_dt[,umap_1:=umap[,1],]
umap_dt[,umap_2:=umap[,2],]
umap_dt[,color_class:=factor(color_class,levels=names(class_colors)),]

umap_dt_genomes=umap_dt[,.(umap_1=mean(umap_1),umap_2=mean(umap_2)),by="mapped_genome"]

#for scrambled (uncomment scrambled above)
#pdf("umap_scrambled.pdf", 10,4)
#ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=color_class))+
#  geom_text(aes(label=rn),size=0.7)+
#  scale_color_manual(values=class_colors)+coord_fixed()
#dev.off()

pdf("umap.pdf", 10,4)
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=color_class))+
  geom_text(aes(label=rn),size=0.7)+
  geom_text(data=umap_dt_genomes,aes(label=mapped_genome),color="black",size=2)+
  scale_color_manual(values=class_colors)+coord_fixed()
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=tissue))+
  geom_text(aes(label=rn),size=0.7)+coord_fixed()+
geom_text(data=umap_dt_genomes,aes(label=mapped_genome),color="black",size=2)
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=mapped_genome))+
  geom_text(aes(label=rn),size=0.7)+coord_fixed()+
geom_text(data=umap_dt_genomes,aes(label=mapped_genome),color="black",size=2)
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=dip))+
  geom_text(aes(label=rn),size=0.7)+scale_color_gradient(high='red',low='blue')+coord_fixed()
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=perc_mapped_passed))+
  geom_text(aes(label=rn),size=0.7)+scale_color_gradient(high='red',low='blue')+coord_fixed()
dev.off()


####gene-wise clustering (interpretation of leiden clusters i.e. gene set enrichment, clustering of tissues/classes is hard )
mat_red=mat_num[rs>50,cs>400] 
gene_annot_red=gene_annot[colnames(mat_red),]
sample_annot_red=sample_annot[rownames(mat_red),]

#one more filtering step to frequent tissues
tissue_table=table(sample_annot_red$tissue)
sample_annot_red=sample_annot_red[sample_annot_red$tissue%in%names(tissue_table[tissue_table>=8]),]
#sample_annot_red=sample_annot_red[sample_annot_red$tissue%in%c("H","L"),]

mat_red=mat_red[row.names(sample_annot_red),]

dim(mat_red)
sum(is.na(mat_red))/prod(dim(mat_red))

cor=cor(mat_red,use = "pairwise.complete.obs")
row.names(cor)=colnames(mat_red)
colnames(cor)=colnames(mat_red)
sum(is.na(cor))

#plot umap
set.seed(5)
umap <- umap(cor, n_neighbors = 15, learning_rate = 0.5, init = "random", n_epochs = 500,min_dist=0.05,spread=1.5)
umap_dt=as.data.table(gene_annot_red,keep.rownames = TRUE)
umap_dt[,umap_1:=umap[,1],]
umap_dt[,umap_2:=umap[,2],]

#leiden clustering
adj=dist(data.frame(umap_dt[,c("rn","umap_1","umap_2")],row.names="rn"))

graph <- graph.adjacency(exp(-as.matrix(adj)), weighted=TRUE, mode="undirected")
leiden=cluster_leiden(graph,resolution_parameter = 0.06,n_iterations=10)
sort(table(leiden$membership))
stopifnot(leiden$names==umap_dt$rn)
umap_dt$leiden=leiden$membership

pdf("umap_genes.pdf", 10,4)
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=N_species))+geom_text(aes(label=hgnc_symbol),size=0.7)+
  scale_color_gradient(high='red',low='blue')+coord_fixed()
ggplot(umap_dt,aes(x=umap_1,y=umap_2,color=as.factor(leiden)))+geom_text(aes(label=hgnc_symbol),size=0.7)+
  guides(color = guide_legend(override.aes = list(size = 4)))+coord_fixed()
dev.off()


library("gprofiler2")
gost_enrich=function(fg,sp){
  fg=fg[!is.na(fg)]
  print(sp)
  if(length(fg)<5){ return(NULL)}
  gp=gost(query = fg,evcodes=F,user_threshold=0.05,correction_method='fdr',source=c("GO:CC","GO:MF","GO:BP") )
  return(gp$result)
}

cl_enrich=umap_dt[,as.data.table(gost_enrich(fg=unique(hgnc_symbol),bg=bg,sp=leiden)),by=c('leiden')]
cl_enrich=cl_enrich[order(p_value)]
cl_enrich[,N:=.N,by="term_name"]

cl_enrich_filt=cl_enrich[,.SD[1:3],by="leiden"][!is.na(term_name)]
cl_enrich_filt[,order_cl:=min(leiden),by="term_name"]
cl_enrich_filt[,term_name:=factor(term_name,levels=unique(term_name[order(order_cl)])),]

pdf("umap_genes_enrich.pdf",7,2.5)
ggplot(cl_enrich_filt,aes(x=leiden,y=term_name,fill=-log10(p_value)))+geom_tile()+scale_fill_gradient(low="blue",high="gold",na.value = "grey")
dev.off()

umap_profiles=merge(umap_dt,profiles_prom,by=c("GeneID","hgnc_symbol"))
#for class
umap_profiles_red=umap_profiles[sample%in%row.names(mat_red)&dip>0,.(meth=mean(meth)),by=c("leiden","color_class")]
umap_profiles_red[,meth_z:=(meth-mean(meth))/sd(meth),by=color_class]

umap_profiles_mat=data.frame(dcast(umap_profiles_red,leiden~color_class,value.var = "meth_z"),row.names = 'leiden')

umap_profiles_red[,color_class:=factor(color_class,levels=names(class_colors)),]

pdf("umap_genes_leiden_meth.pdf",4,2.5)
ggplot(umap_profiles_red,aes(x=leiden,y=color_class,fill=meth_z))+geom_tile()+scale_fill_gradient2(high="red",low="blue",mid="white")
dev.off()

#for tissue
umap_profiles_red=umap_profiles[sample%in%row.names(mat_red)&dip>0,.(meth=mean(meth)),by=c("leiden","tissue")]
umap_profiles_red[,meth_z:=(meth-mean(meth))/sd(meth),by=tissue]

umap_profiles_mat=data.frame(dcast(umap_profiles_red,leiden~tissue,value.var = "meth_z"),row.names = 'leiden')
hcl=hclust(dist(umap_profiles_mat))
hcl_cl=hclust(dist(t(umap_profiles_mat)))

umap_profiles_red[,tissue:=factor(tissue,levels=hcl_cl$labels[hcl_cl$order]),]

pdf("umap_genes_leiden_tissues_meth.pdf",4,3)
ggplot(umap_profiles_red,aes(x=leiden,y=tissue,fill=meth_z))+geom_tile()+scale_fill_gradient2(high="red",low="blue",mid="white")
dev.off()


####Predict sample features across species based on the ortholog gene space

##Functions
get_pred=function(mat,tissues=c('H','L'),classes=c('Mammalia'),pred_col='tissue',min_dip=0.6,iter=100,train_spec=100,scramble=FALSE){
  annot_cols=c('sample', 'color_class', 'tissue', 'species', 'mapped_genome', 'dip', 'perc_mapped_passed')  
  mat_num_red=mat[tissue%in%tissues&dip>min_dip&color_class%in%classes,]
  cs=colSums(!is.na(mat_num_red))/nrow(mat_num_red)
  rs=rowSums(!is.na(mat_num_red))/ncol(mat_num_red)
  mat_num_red=mat_num_red[rs>0.12,cs>0.6,with=F]#0.6 #0.08
  mat_num_df=mat_num_red[,-annot_cols[!annot_cols%in%pred_col],with=FALSE]
  mat_annot_df=mat_num_red[,annot_cols[!annot_cols%in%pred_col],with=FALSE]
  print(dim(mat_num_df))
  
  if (scramble==TRUE){
    ##make mock with scrambled methylation
    mat_num_df1=copy(mat_num_df)
    
    cl=mat_num_df[,pred_col,with=FALSE]
    vals=as.numeric(unlist(mat_num_df))
    vals=vals[!is.na(vals)]
    mat_num_df[!is.na(mat_num_df)]=sample(vals,replace = TRUE,sum(!is.na(mat_num_df)))
    ##mat_num_df[is.na(mat_num_df)]=sample(0:100,replace=TRUE,sum(is.na(mat_num_df)))
    mat_num_df[,(pred_col):=cl,]
    stopifnot(all(is.na(mat_num_df)==is.na(mat_num_df1)))
  }
  
  pred_red=data.table()
  imp=data.table()
  for (i in c(1:iter)){
    #train_idx= which(mat_annot_df$color_class=='Mammalia')#choose all the birds as train
    #train_idx=sample(1:nrow(mat_num_df),300) #choose train randomly
    set.seed(i)
    train_idx= which(mat_annot_df$species%in%sample(unique(mat_annot_df$species),train_spec))#choose  some bird species as train

    train=mat_num_df[train_idx,]
    test=mat_num_df[-train_idx,]

    fit <- rpart(as.formula(paste0(pred_col,'~.')), data = train)

    pred_tmp=data.table(pred=predict(fit,test, type = "prob")[, 2], true=unlist(test[,pred_col,with=FALSE]))
    pred_tmp[,i:=i]
    pred_red=rbindlist(list(pred_red,pred_tmp))

    imp_tmp=as.data.table(fit$variable.importance,keep.rownames = TRUE)
    imp_tmp[,i:=i]
    imp=rbindlist(list(imp,imp_tmp))

  }

  pred <- prediction(pred_red$pred, pred_red$true)
  auc_ROCR <- performance(pred, measure = "auc")
  perf=performance(pred, "tpr", "fpr")

  return(list(roc=data.table(fpr=unlist(perf@x.values),tpr=unlist(perf@y.values),auc=signif(auc_ROCR@y.values[[1]],3)),imp=imp,bg=sub("_.*","",colnames(mat_num_df))))
}


make_summary=function(preds,cols,annots){
  roc_list=list()
  imp_list=list()
  bg=c()
  for (i in 1:length(preds)){
    
    annot_text=paste(cols,':=',annots[[i]])
    roc_list[[i]]=preds[[i]]$roc[,eval(parse(text=annot_text)),]
    imp_list[[i]]=preds[[i]]$imp[,eval(parse(text=annot_text)),]
    bg=unique(c(bg,preds[[i]]$bg))
  
  }
  
  roc=rbindlist(roc_list)
  imp=rbindlist(imp_list)
  return(list(roc=roc,imp=imp,bg=bg))
}

assess_imp=function(all_imp,mat,N_thres=40,mean_thres=5,rel_thres=0.01,test_col='color_class',comp_col='tissue',g1='H',g2='L',topN=4){
  imp_red=all_imp[,.(mean=mean(V2),sd=sd(V2),N=.N),by=c('V1','group')]
  imp_red[,mean_rel:=mean/sum(mean),by='group']
  
  imp_sub=imp_red[N>N_thres&mean>mean_thres&mean_rel>rel_thres]
  imp_sub[,order:=rank(-mean),by='group']
  imp_sub=imp_sub[order<=topN]
  imp_sub[,V1:=factor(V1,levels=unique(V1[order(group,N)])),]

  
  pl=ggplot(imp_sub,aes(y=V1,x=group,fill=N))+geom_tile()+
    scale_fill_gradient(low='lightblue',high='blue')+
    geom_text(aes(label=signif(mean,3)))+
  xlab('')+ylab('')
  
  pred_genes=mat[color_class%in%c('Mammalia','Aves','Reptilia')&tissue%in%c('H','L'),]
  
  pred_genes_long=melt(pred_genes,id.vars = c('color_class','tissue','species','sample'),measure.vars = as.character(imp_sub$V1),variable.name ='gene' ,value.name = 'meth')
  pred_genes_long=merge(pred_genes_long,imp_sub,by.x='gene',by.y='V1',allow.cartesian=TRUE)
  pred_genes_long[,gene:=sub('_.*','',gene),]
  pred_genes_long[,facet:=paste0(group,': ',gene,' ',signif(mean,3)),]
  pred_genes_long[,facet:=factor(facet,levels=unique(facet[order(order)])),]
  pred_genes_long[,color_class:=factor(color_class,levels=names(class_colors))]
  
  
  txt1=paste0(comp_col,"==",g1)
  txt2=paste0(comp_col,"==",g2)
  print(txt1)
  print(txt2)
  
  
  get_pval=function(x,y){
    if(sum(!is.na(x))>2&sum(!is.na(y))>2){
      p=wilcox.test(x = x,y=y)$p.value
    }else(return(1))
    
    return(signif(p,3))
  }
  
  pvals=pred_genes_long[,.(pval=get_pval(
    x = meth[eval(parse(text=paste0(comp_col,"==",g1)))],
    y=meth[eval(parse(text=paste0(comp_col,"==",g2)))])),
    by=c('gene',test_col,'facet')]
  
  return(list(genes=pred_genes_long,pvals=pvals,pl=pl))

}

### stat prediction analysis
##predict class
#pre-filter mat to only species with the same number of heart and liver samples
HL_count=mat[,.(N_L=sum(tissue=='L'),N_H=sum(tissue=='H')),by=c('species')]

H_MA=get_pred(mat[species%in%HL_count[N_L==N_H]$species],tissues=c('H'),classes=c('Mammalia','Aves'),pred_col='color_class',min_dip=0.6,iter=100,train_spec=90)
L_MA=get_pred(mat[species%in%HL_count[N_L==N_H]$species],tissues=c('L'),classes=c('Mammalia','Aves'),pred_col='color_class',min_dip=0.6,iter=100,train_spec=90)
H_MA_s=get_pred(mat[species%in%HL_count[N_L==N_H]$species],tissues=c('H'),classes=c('Mammalia','Aves'),pred_col='color_class',min_dip=0.6,iter=100,train_spec=90,scramble = T)
L_MA_s=get_pred(mat[species%in%HL_count[N_L==N_H]$species],tissues=c('L'),classes=c('Mammalia','Aves'),pred_col='color_class',min_dip=0.6,iter=100,train_spec=90, scramble = T)

MA=make_summary(list(H_MA,L_MA,H_MA_s,L_MA_s),"c('group','type')",list("list('H','t')","list('L','t')","list('H','s')","list('L','s')"))

ann_text=paste0(unique(MA$roc[,c('group','type','auc')])[,.(lab=paste0(group,': ',auc[type=='t'],'/',auc[type=='s'])),by='group']$lab,collapse = '\n')

pdf('pred_roc_MA.pdf',6,4)
ggplot(MA$roc,aes(x=fpr,y=tpr,color=group,lty=type,group=paste0(group,type)))+
  geom_line(lwd=1)+annotate(geom='text',x=1,y=0,label=ann_text, col='black',vjust=-1,hjust=1)+
  scale_linetype_manual(values = c('s'=21,'t'=1))+
  scale_color_manual(values = tissue_colors[names(tissue_colors)%in%c('H','L')])+
  xlab('False positive rate')+ylab('True positive rate')+coord_fixed()
dev.off()

MA_imp=assess_imp(MA$imp[type!='s'],mat[species%in%HL_count[N_L==N_H]$species],comp_col='color_class',test_col = 'tissue',g1="'Mammalia'",g2="'Aves'")

MA_imp$pl

MA_imp_genes_red=MA_imp$genes[,.(facet=paste0(unique(facet),collapse = '\n')),by=c('gene','tissue','color_class','meth','species','sample')]
MA_imp_pvals_red=MA_imp$pvals[,.(facet=paste0(unique(facet),collapse = '\n')),by=c('gene','tissue','pval')]

pdf('pred_genes_MA.pdf',9,6)
ggplot(MA_imp_genes_red,
       aes(x=tissue,y=meth,color=color_class))+
  #  geom_point(aes(group=get(sel_col)),position = position_jitterdodge())+
  geom_text(aes(label=species),position = position_jitterdodge(jitter.width = 0.3),alpha=0.8,size=1)+
  geom_boxplot(outlier.colour = NA,fill='transparent')+
  geom_text(data = MA_imp_pvals_red,size=3,aes(x=tissue,y=105,label=paste0('p=',formatC(pval, format = "e", digits = 2))),color='black')+
  ylab('CpG methylation (%)')+facet_wrap(~facet,ncol=3)+
  scale_color_manual(values = class_colors[names(class_colors)%in%c('Mammalia','Aves','Reptilia')])+
  stat_summary(fun.data = give.n,fun.args = c(y=-6), geom = "text",size=2.5,position = position_dodge(width = 1))+
  theme(legend.position="bottom")
dev.off()

##predict tissue
HL_M=get_pred(mat,tissues=c('H','L'),classes=c('Mammalia'),pred_col='tissue',min_dip=0.6,iter=100,train_spec=80)
HL_M_s=get_pred(mat,tissues=c('H','L'),classes=c('Mammalia'),pred_col='tissue',min_dip=0.6,iter=100,train_spec=80,scramble = T)
HL_A=get_pred(mat,tissues=c('H','L'),classes=c('Aves'),pred_col='tissue',min_dip=0.6,iter=100,train_spec=80)
HL_A_s=get_pred(mat,tissues=c('H','L'),classes=c('Aves'),pred_col='tissue',min_dip=0.6,iter=100,train_spec=80,scramble = T)

HL=make_summary(list(HL_M,HL_A,HL_M_s,HL_A_s),"c('group','type')",list("list('Mammalia','t')","list('Aves','t')","list('Mammalia','s')","list('Aves','s')"))

ann_text=paste0(unique(HL$roc[,c('group','type','auc')])[,.(lab=paste0(group,': ',auc[type=='t'],'/',auc[type=='s'])),by='group']$lab,collapse = '\n')

pdf('pred_roc_HL.pdf',6,4)
ggplot(HL$roc,aes(x=fpr,y=tpr,color=group,lty=type,group=paste0(group,type)))+
  geom_line(lwd=1)+annotate(geom='text',x=1,y=0,label=ann_text, col='black',vjust=-1,hjust=1)+
  scale_linetype_manual(values = c('s'=21,'t'=1))+
  scale_color_manual(values = class_colors[names(class_colors)%in%c('Mammalia','Aves')])+
  xlab('False positive rate')+ylab('True positive rate')+coord_fixed()
dev.off()

HL_imp=assess_imp(HL$imp[type!='s'],mat,comp_col='tissue',test_col = 'color_class',g1="'H'",g2="'L'")

HL_imp$pl

pdf('pred_genes_HL.pdf',9,6)
ggplot(HL_imp$genes,
       aes(x=color_class ,y=meth,color=tissue))+
  #  geom_point(aes(group=get(sel_col)),position = position_jitterdodge())+
  geom_text(aes(label=species),position = position_jitterdodge(jitter.width = 0.1),alpha=0.8,size=1)+
  geom_boxplot(outlier.colour = NA,fill='transparent')+
  geom_text(data = HL_imp$pvals,size=2.5,aes(x=color_class,y=105,label=paste0('p=',formatC(pval, format = "e", digits = 2))),color='black')+
  ylab('CpG methylation (%)')+facet_wrap(~facet,ncol=4)+
  scale_color_manual(values = tissue_colors[names(tissue_colors)%in%c('H','L')])+
  stat_summary(fun.data = give.n,fun.args = c(y=-6), geom = "text",size=2.5,position = position_dodge(width = 1))+
  theme(legend.position="bottom")
dev.off()


####plot methylation levels at specific genes (genes that are captured in all classes)
profiles_prom[,detailed_classes:=ifelse(species%in%c("JL"),species,color_class),]
med_meth=profiles_prom[,.(median_meth=median(meth),N=.N),
                       by=c('detailed_classes',"color_class",'hgnc_symbol','GeneID')][,c('N_classes','q_25'):=list(length(unique(detailed_classes)),quantile(N,0.25)),by=c('hgnc_symbol','GeneID')]
med_meth[,detailed_classes:=factor(detailed_classes,levels=c("Invertebrata","FLA","JL",names(class_colors)[-1])),]

med_meth_sub=med_meth[N_classes>=8&q_25>3]

med_meth_sub[,overall_mean:=median(median_meth),by="hgnc_symbol"]
med_meth_sub[,N_high:=sum(median_meth>0.8),by="hgnc_symbol"]

med_meth_sub[,hgnc_symbol:=factor(hgnc_symbol,levels=unique(hgnc_symbol[order(overall_mean)])),]

pdf("common_genes_meth.pdf",6.5,7.5)
ggplot(med_meth_sub,aes(y=hgnc_symbol,x=detailed_classes,fill=median_meth))+
  geom_tile()+geom_text(aes(label=N,size=N))+
  scale_fill_gradient(low='dodgerblue',high='gold')+scale_size_continuous(range=c(2.5,4.5))+
  scale_x_discrete(labels=c("JL"="Lamprey",class_short))+
  xlab("")+ylab("")
dev.off()
