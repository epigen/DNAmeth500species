#!/usr/bin/env Rscript
library(pryr)

kmer=NA

if (is.na(kmer)){kmer=c(1:10)}
args = commandArgs(trailingOnly=TRUE)
species=args[1]
if (is.na(kmer)){kmer=c(1:10)}
subdir_name="sample_check"
source(file.path(Sys.getenv("CODEBASE"),"compEpi/src/02.0_predict_initialization.R"))


#species = "WHH"
wd=file.path(processed_dir)
mywd="/scratch/lab_bock/shared/projects/compEpi/results_analysis/02_predict_meth"
setwd(wd)
#previously results were stored in the results_pipeline folder (paste0("/scratch/lab_bock/jklughammer/projects/compEpi_uc/results_pipeline/",species,"/methPred/",subdir_name)), now moving to results_analysis (more appropriate location) 
subdir=file.path(mywd,"02.4_verify_inverted",subdir_name,species)
dir.create(subdir,recursive=TRUE)
dir.create(file.path(subdir, "sequences"), recursive=TRUE)
dir.create(file.path(subdir, "stats"), recursive=TRUE)

##training on each tissue
load_data<-function(SP){
  #load data
  meth_data_mean=fread(list.files(path=paste0(SP,"/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/"),pattern=".*_mean_meth.tsv",full.names=TRUE))
  
  uc_map=fread(list.files(path=paste0(SP,"/toSelf_filtered_0.08mm_final_concat"),pattern=".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",recursive=TRUE,full.names=TRUE))
  ded_ref=readDNAStringSet(paste0(SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  
  ded_ref_cgg=as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names=names(ded_ref[ded_ref_cgg]) 
  
  #load annotation
  sampleAnnotation<-fread(paste0(SP,"/meta/sampleAnnotation.tsv"))
  sampleAnnotation<-sampleAnnotation[sampleAnnotation$Select==1]
  
  sampleAnnotation$isUC<-as.character(lapply(strsplit(sampleAnnotation$Sample_Name, split="_"), "[", 4))
  
  sampleAnnotation[is.na(sampleAnnotation)]<-0
  sampleAnnotation<-sampleAnnotation[sampleAnnotation$isUC==0]
  
  sampleAnnotation$rep <- as.character(lapply(strsplit(sampleAnnotation$Sample_Name, split="_"), "[", 2))
  samples<-unique(sampleAnnotation$rep)
  sample_ids<-lapply(samples, function(x) sampleAnnotation[sampleAnnotation$rep==x]$Sample_Name)
  print(paste0("data read about ", SP))
  return(list(meth_data_mean=meth_data_mean, ded_ref=ded_ref, cgg_names=cgg_names, uc_map=uc_map, samples=samples, sample_ids=sample_ids))
}

get_dataset<-function(ids, meth_data_mean, ded_ref, cgg_names, uc_map, isTest=FALSE, train_ids=NA){
  
  ## subselection of columns for the id
  names<-c(c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),c(paste0(ids, ".cov"), paste0(ids, ".meth")))
  sub_df<-meth_data_mean[, ..names]
  sub_df_long<-melt(sub_df,id.vars=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),measure=patterns(".cov", ".meth"),variable.factor=FALSE,variable.name="sample",value.name=c("cov","meth"),na.rm=TRUE)
  sub_df_cond<-sub_df_long[,list(Nsamples=.N,mean_cov=mean(cov),min_cov=min(cov),max_cov=max(cov),mean_meth=mean(meth),min_meth=min(meth),max_meth=max(meth)),by=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
  sub_df_cond_red<-sub_df_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
  sub_df_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
  
  #making sure that we don have same reference ids in test and train datasets
  if(isTest){
    print("filtering")
    idx<-!unlist(lapply(sub_df_cond_red$meta, is.element, set=train_ids))
    sub_df_cond_red<-sub_df_cond_red[idx,]
    N=1000
  }else{
    N=2000
  }
  
  ds<-create_datasets(sub_df_cond_red, ded_ref, strsplit(ids[[1]], "_")[[1]][3], TO_SAVE=!FALSE, N = N )
  return(ds)
}

##CHANGE FOR PARALEL ON THE FULL DATASET
data <- load_data(species[[1]])

if(isDev){
  train_sample_ids<-data$sample_ids[c(1:2)]
  train_samples<-data$samples[c(1:2)]
}else{
  train_sample_ids<-data$sample_ids
  train_samples<-data$samples
}
train_data <- mclapply(data$sample_ids, get_dataset, meth_data_mean=data$meth_data_mean, cgg_names=data$cgg_names, ded_ref=data$ded_ref, uc_map=data$uc_map)

##getting train_test models for all the tissues:
train_models <- mclapply(seq_along(train_samples), function(i) {simpleCache(cacheName=paste0("methTrain_",train_samples[[i]]),instruction = {train_test(train_data[[i]]$x, train_data[[i]]$y, type=train_sample_ids[[i]][[1]], ifRand=paste0("noRand", "Train"), run=0, k=kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},cacheDir=paste0(subdir,"/RCache"),assignToVariable="res",recreate=FALSE)})

##training on random labels for comparison
train_models_rand<-mclapply(seq_along(train_samples), function(i) {simpleCache(cacheName=paste0("methTrainRand_",train_samples[[i]]),instruction = {train_test(train_data[[i]]$x, sample(train_data[[i]]$y), type=train_sample_ids[[i]][[1]], ifRand=paste0("rand", "Train"), run=0, k=kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},cacheDir=paste0(subdir,"/RCache"),assignToVariable="resRand",recreate=FALSE)}) 

##test on the average methylation values for all the species (or on a subset for development)
test_species<-read.table("/scratch/lab_bock/dromanovskaia/CompEpi/compEpi/meta/species_list.txt")$V1

if(isDev){
  #taking random species
  test_species<-c(as.character(sample(test_species)[c(2:12)]),species)
}


in_other_sample<-function(SP){
  print(paste0("reading data about ", SP))
  meth_data_mean = fread(list.files(path = paste0(SP, "/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/"),pattern = ".*_mean_meth.tsv",
                                    full.names = TRUE))
  uc_map = fread(list.files(path = paste0(SP, "/toSelf_filtered_0.08mm_final_concat"),
                            pattern = ".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",recursive = TRUE,full.names = TRUE))
  ded_ref = readDNAStringSet(paste0(SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  ded_ref_cgg = as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names = names(ded_ref[ded_ref_cgg]) 
  
  ##all we need is an average dataset
  meth_data_mean_long = melt(meth_data_mean,id.vars = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),
                             measure = patterns(".cov", ".meth"),variable.factor = FALSE,variable.name = "sample",
                             value.name = c("cov","meth"),na.rm = TRUE)
  
  meth_data_mean_cond = meth_data_mean_long[,list(Nsamples = .N,mean_cov = mean(cov),min_cov = min(cov),max_cov = max(cov),
                                                  mean_meth = mean(meth),min_meth = min(meth),max_meth = max(meth)),
                                            by = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
  
  meth_data_mean_cond_red = meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
  
  meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
  
  if(NROW(meth_data_mean_cond_red[meth_data_mean_cond_red$category==1, ])==0 | NROW(meth_data_mean_cond_red[meth_data_mean_cond_red$category==-1, ])==0){
    return(0)
  }else{
  test_ds <- create_datasets(meth_data_mean_cond_red = meth_data_mean_cond_red, ded_ref = ded_ref, SP, FALSE, 2000)
  
  preddec_list <- lapply(seq_along(train_models), function(i) predict(test_ds$x, object = train_models[[i]]$model[[1]]))
  pref <- lapply(preddec_list, evaluatePrediction, test_ds$y,allLabels=c(-1,1), print=FALSE)
  f1_list <- lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec <- lapply(seq_along(train_models), function(i) predict(test_ds$x, object=train_models[[i]]$model[[1]], predictionType="decision"))
  rocdata_list <- lapply(preddec_list_dec, computeROCandAUC, test_ds$y, allLabels=c(-1,1))
  roc_dt_list=lapply(seq_along(rocdata_list), function(i) data.table(fdr=unlist(rocdata_list[[i]]@FPR),tpr=unlist(rocdata_list[[i]]@TPR),auc=unlist(rocdata_list[[i]]@AUC),f1=f1_list[[i]], run=0, type=SP,ifRand="noRandTest", k=0, C=0,min_motif="NN",max_motif="NN", sample_N=train_samples[[i]]))
  roc_dt_sample<-rbindlist(roc_dt_list)
  
  ##prediction the random model on random labels
  #building predictions
  y_rand <- sample(test_ds$y)
  preddec_list_rand <- lapply(seq_along(train_models), function(i) predict(test_ds$x, object = train_models_rand[[i]]$model[[1]]))
  pref_rand <- lapply(preddec_list_rand, evaluatePrediction, y_rand, allLabels=c(-1,1), print=FALSE)
  f1_list <- lapply(pref_rand, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec_rand <- lapply(seq_along(train_models_rand), function(i) predict(test_ds$x, object=train_models_rand[[i]]$model[[1]], predictionType="decision"))
  rocdata_list_rand <- lapply(preddec_list_dec_rand, computeROCandAUC, y_rand, allLabels=c(-1,1))
  roc_dt_list_rand=lapply(seq_along(rocdata_list_rand), function(i) data.table(fdr=unlist(rocdata_list_rand[[i]]@FPR),tpr=unlist(rocdata_list[[i]]@TPR),auc=unlist(rocdata_list_rand[[i]]@AUC),f1=f1_list[[i]], run=0, type=SP,ifRand="RandTest", k=0, C=0,min_motif="NN",max_motif="NN", sample_N=train_samples[[i]]))
  roc_dt_sample_rand<-rbindlist(roc_dt_list_rand)
  
  rm(test_ds)
  return(list(pred=roc_dt_sample, pred_rand=roc_dt_sample_rand))
  }
}


res_list <- sapply(test_species, in_other_sample, simplify = F)
print(sapply(res_list, length))
res_list <- res_list[sapply(res_list, length) > 1]
res <- rbindlist(lapply(res_list, function(x) x$pred))
res_rand <- rbindlist(lapply(res_list, function(x) x$pred_rand))
print(mem_used())

res_full <- rbind(res, res_rand)
write.table(res_full, file.path(subdir, "all_data.csv"), row.names = F, quote = F)

res_sum <- res_full %>%
  dplyr::group_by(type, ifRand, sample_N) %>% 
  dplyr::summarise(AUC = first(auc), F1 = first(f1))
self_auc <- res_sum[res_sum$type==species & res_sum$ifRand=="noRandTest", ]
res_sum$color_class <- unlist(sapply(res_sum$type, function(x) sp_df[x, ]$color_class))
res_sum$class_short <- unlist(sapply(res_sum$color_class, function(x) class_short[[x]]))
res_sum$class_short <- factor(res_sum$class_short, levels = class_short)

p <- ggplot(res_sum, aes(x = class_short, y = AUC, fill = ifRand)) + geom_boxplot(outlier.shape = 21) +
  facet_wrap(~sample_N, ncol = 1) + xlab("taxonomy group") + 
  scale_fill_manual(values = c("noRandTest"="#4682B4", "RandTest" = "grey")) + 
  geom_hline(data =  self_auc, mapping = aes(yintercept = AUC), linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey", alpha = 0.5)

ggsave(file.path(subdir,  "AUC_in_other_sp.pdf"), p, width = 6, height = 2*length(unique(res_sum$sample_N)), device = "pdf")
