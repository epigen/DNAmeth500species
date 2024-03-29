source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

## reading in which species we are working with
options <- commandArgs(trailingOnly = TRUE)  
species = options[1]

wd=file.path(analysis_dir, "01_basicStats", "01.9_coverage_analysis", species)
dir.create(wd)
setwd(wd)

meth_file=file.path(processed_dir, species, "toSelf_filtered_0.08mm_final_concat","diffMeth_cpg", paste0(species, "_mean_meth.tsv"))
meth=fread(meth_file)


meth_long=melt(meth,id.vars = c("meta", "score", "consN", "mean_diffNts", "max_diffNts", "min_diffNts"),measure.vars = patterns(meth="*.meth",cov="*.cov"),variable.name = "sample",na.rm = TRUE)


meth_long$sample <- unique(sub(".meth|.cov","",grep(".meth|.cov",names(meth),value = TRUE)))[meth_long$sample]

meth_long[,individual:=unlist(strsplit(sample,"_"))[2],by="sample"]
meth_long[,tissue:=unlist(strsplit(sample,"_"))[3],by="sample"]

#use average/mean because it represents the coverage each fragment 
#would have if it was evenly distributed.
#This is calculated per sample and the classification of each fragment is done 
#relative to this value.
meth_long[,norm_cov:=cov/mean(cov),by="sample"] 
meth_long[,avg_cov:=sum(cov)/.N,by="sample"] 


#check how the distribution looks (no easy natural thresholds)
pdf("cov_distibution.pdf", height=6,width=10)
ggplot(meth_long,aes(x=log2(norm_cov+10^-2),fill=sample))+geom_histogram(bins = 100)+facet_wrap(~sample,scale="free")+geom_vline(xintercept = c(-4:5))
dev.off()


#classify fragments according to their coverage relative to the average coverage in each sample
cov_avg_factor_low=0.5
cov_avg_factor_hi=4
meth_long[,Min:=cov>=avg_cov*cov_avg_factor_low,]
meth_long[,Rep:=cov>=avg_cov*cov_avg_factor_hi,]
#calculate in how many samples a fragment is classified as covered (min) and/or repetitive (rep) 
meth_long[,nMin:=sum(Min),by="meta"]
meth_long[,nRep:=sum(Rep),by="meta"]

#Mark fragments that are covered in >80% (pos) of one individual's samples or < 20% (neg)
meth_long[,Min_indi_pos:=sum(Min)>.N*0.8&Min==TRUE,by=c("meta","individual")]
meth_long[,Min_indi_neg:=sum(Min)<=max(.N*0.2,1)&sum(Min)>0&Min==TRUE,by=c("meta","individual")] ## add exclusive overlap among differnet individuals as neg. control

#Mark fragments that are covered only in one individual or across individuals (bg)
indivN=length(unique(meth_long$individual))
meth_long[,Min_indi:=length(unique(individual[Min==TRUE]))<=max(0.2*indivN,1)&Min_indi_pos==TRUE,by=c("meta")]
meth_long[,Min_indi_bg:=length(unique(individual[Min==TRUE]))>max(0.2*indivN,1)&Min_indi_neg==TRUE,by=c("meta")] 

#count number of fragments in each sample that are classified as repeat, single, individual, individual_bg, or amplified
totalN=length(unique(meth_long$sample))
stats=meth_long[,.(N_frags=.N,total_cov=sum(cov),N_rep=sum(Rep==TRUE&nRep>=totalN*0.8),N_single=sum(Min==TRUE&nMin<=max(c(totalN*0.2,1))),N_ind=sum(Min_indi==TRUE),N_ind_bg=sum(Min_indi_bg==TRUE),N_amp=sum(Rep==TRUE&nRep<=max(c(totalN*0.2,1)))),by=c("sample")]

write.table(stats,"coverage_stats.tsv", row.names=FALSE, sep="\t", quote=FALSE)
