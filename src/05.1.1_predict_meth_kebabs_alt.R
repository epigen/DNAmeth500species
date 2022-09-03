#!/usr/bin/env Rscript
#classical initiation
kmer=NA

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/05.0_predict_initialization.R"))
args = commandArgs(trailingOnly=TRUE)
species=args[1]
subdir_name=args[2]
kmer=as.numeric(args[3])

if (is.na(kmer)){kmer=c(1:10)}
sessionInfo()


options(cores=1)

subdir_name="screen_alt"

wd=file.path(processed_dir,species)
setwd(wd)

subdir=file.path(analysis_dir,"05_predict_meth/05.1_within_species",subdir_name,species)
print(subdir)
dir.create(subdir,recursive=TRUE)
dir.create(file.path(subdir, "sequences"), recursive=TRUE)
dir.create(file.path(subdir, "stats"), recursive=TRUE)

#load data
meth_data_mean=fread(list.files(path="toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/",pattern=".*_mean_meth.tsv",full.names=TRUE))
uc_map=fread(list.files(path="toSelf_filtered_0.08mm_final_concat",
                        pattern=".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",recursive=TRUE,full.names=TRUE))
ded_ref=readDNAStringSet("reduced/consensus//toSelf_filtered_0.08mm_final.fa")

ded_ref_cgg=as.character(substring(ded_ref,1,3))=="CGG"
cgg_names=names(ded_ref[ded_ref_cgg])

#select features (methylated/unmethylated)
meth_data_mean_long=melt(meth_data_mean,id.vars=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),
                         measure=patterns(".cov", ".meth"),variable.factor=FALSE,variable.name="sample",
                         value.name=c("cov","meth"),na.rm=TRUE)
meth_data_mean_cond=meth_data_mean_long[,list(Nsamples=.N,mean_cov=mean(cov),min_cov=min(cov),max_cov=max(cov),
                                              mean_meth=mean(meth),min_meth=min(meth),max_meth=max(meth)),
                                        by=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
meth_data_mean_cond_red=meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
##CHANGED: label 1 criteria stays the same, while label -1 is assign if AT LEAST one of the 
##samples has methylation levels lower than 20
meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(min_meth<20,-1,NA)),]

res<-create_datasets(meth_data_mean_cond_red, ded_ref, species, SEED = FALSE)
#minimum features of either group --> maximum number of features for each group
x<-res$x
y<-res$y


#original labels
simpleCache(cacheName="methPred_noRand_uc",instruction={train_test(x=x,y=y,type=species, 
                                                      ifRand='noRand',run=0,k=kmer, subdir=subdir)},
            cacheDir=paste0(subdir,"/RCache"),assignToVariable="res",recreate=FALSE)
roc_res=res$roc_dt




#random labels (several iterations)
#rand_labs=lapply(seq_len(5), function(x) {return(list(name=x,y_rand=sample(y)))})
#CHANGED to one iteration:
y_rand = sample(y)

rand_res={simpleCache(cacheName="methPred_Rand_uc",instruction={train_test(x=x,y=y_rand,type=species, 
                                                              ifRand='rand',run=1,k=kmer, subdir=subdir)},
                      assignToVariable="rand_res",cacheDir=paste0(subdir,"/RCache"),recreate=FALSE)}
rand_res_dt = rand_res$roc_dt

#merge rand and original + plot
roc_res=rbind(roc_res,rand_res_dt)

auc_res <- roc_res[,list(auc=mean(auc)),by=ifRand]
##uploading the results from the previous:
##uploading the initial AUC values(by stronger method)
auc_res_initial = read.csv(file.path(analysis_dir, "05_predict_meth/05.1_within_species/screen", species, 
                                     paste0(species, "_stats.tsv" )), sep = "\t")
auc_res <- rbind(auc_res, t(data.frame(c("ifRand" = "initial","auc" = as.numeric(auc_res_initial$AUC)))))
auc_res$auc <- as.numeric(auc_res$auc)
auc_res[,x:=0.85,]
auc_res[,y:=c(0.01,0.07, 0.13),]

pdf(paste0(subdir,"/ROC.pdf"),height=4,width=5)
ggplot(roc_res, aes(x=fdr,y=tpr,col=ifRand))+
  geom_line(aes(group=run, alpha=ifRand))+
  geom_text(data=auc_res,aes(x=x,y=y,label=paste0("auc=",round(auc,3))))+
  scale_color_manual(values=c("rand"="grey","noRand"="blue","initial"="red" ))+
  scale_alpha_manual(values=c("rand"=0.5,"noRand"=1, "initial"=1))
dev.off()

ids <- which(stats_annot$species == species & stats_annot$conversion_type == "converted")
roc_res$N_conv_samples <- length(ids)

write.csv(roc_res, paste0(subdir,"/roc_res.csv"))


##uploading the initial AUC values(by stronger method)
if(file.exists(file.path(analysis_dir, "05_predict_meth/05.1_within_species/screen1", species, 
                         "roc_res.csv" ))){
  roc_res_initial = read.csv(file.path(analysis_dir, "05_predict_meth/05.1_within_species/screen1", species, 
                                       "roc_res.csv" ), row.names = 1)
  
  roc_res_initial <- roc_res_initial[roc_res_initial$ifRand=="noRand",]
  roc_res_initial$ifRand <- "initial"
  roc_res_initial$N_conv_samples <- 0
  roc_res <- rbind(roc_res, roc_res_initial, fill = T)
  pdf(paste0(subdir,"/ROC_with_initial.pdf"),height=4,width=5)
  print(ggplot(roc_res, aes(x=fdr,y=tpr,col=ifRand))+
    geom_line(aes(group=run, alpha=ifRand))+
    geom_text(data=auc_res,aes(x=x,y=y,label=paste0("auc=",round(auc,3))))+
    scale_color_manual(values=c("rand"="grey","noRand"="blue","initial"="red" ))+
    scale_alpha_manual(values=c("rand"=0.5,"noRand"=1, "initial"=1)))
  dev.off()
}

