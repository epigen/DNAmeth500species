#!/usr/bin/env Rscript

kmer=NA
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/02.0_predict_initialization.R"))
if (is.na(kmer)){kmer=c(1:10)}
args = commandArgs(trailingOnly=TRUE)
species=args[1]
subdir_name="screen"

library(pryr)


options(cores=8)


wd=file.path(processed_dir)
setwd(wd)

mywd="/scratch/lab_bock/shared/projects/compEpi/results_analysis/02_predict_meth"
subdir=file.path(mywd,"02.3_tissue_vs_sample",subdir_name,species)
dir.create(subdir,recursive=TRUE)

dir.create(file.path(subdir, "sequences"), recursive=TRUE)
dir.create(file.path(subdir, "stats"), recursive=TRUE)

####initializing a color libraary
### creating a tissue color list
tissues<-unique(sampleAnnot$Tissue)
tissues<-c("trained on", "rand", tissues[1:(length(tissues)-2)])

#library(RColorBrewer)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
col_vector<-c("#000000", "#D3D3D3", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
              "#666666", "#1B9E77",  "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", 
              "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 
              "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", 
              "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
              "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854")

color_tissues<-setNames(col_vector, tissues)

#####
load_data<-function(SP){
  #load data
  meth_data_mean=fread(list.files(path=paste0(SP,"/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/"),
                                  pattern=".*_mean_meth.tsv",full.names=TRUE))
  uc_map=fread(list.files(path=paste0(SP,"/toSelf_filtered_0.08mm_final_concat"),
                          pattern=".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",recursive=TRUE,full.names=TRUE))
  ded_ref=readDNAStringSet(paste0(SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  
  ded_ref_cgg=as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names=names(ded_ref[ded_ref_cgg]) 
  
  #load annotation
  sampleAnnotation<-fread(paste0(SP,"/meta/sampleAnnotation.tsv"))
  sampleAnnotation<-sampleAnnotation[sampleAnnotation$Select==1]
  
  sampleAnnotation$isUC<-as.character(lapply(strsplit(sampleAnnotation$Sample_Name, split="_"), "[", 4))
  
  sampleAnnotation[is.na(sampleAnnotation)]<-0
  sampleAnnotation<-sampleAnnotation[sampleAnnotation$isUC==0]
  
  samples<-unique(sampleAnnotation$Tissue)
  sample_ids<-lapply(samples, function(x) sampleAnnotation[sampleAnnotation$Tissue==x]$Sample_Name)
  print(paste0("data read about ", SP))
  return(list(meth_data_mean=meth_data_mean, ded_ref=ded_ref, cgg_names=cgg_names, uc_map=uc_map, samples=samples, 
              sample_ids=sample_ids))
}

get_dataset<-function(ids, meth_data_mean, ded_ref, cgg_names, uc_map, isTest=FALSE, train_ids=NA){
  
  ## subselection of columns for the id
  names<-c(c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),c(paste0(ids, ".cov"), paste0(ids, ".meth")))
  sub_df<-meth_data_mean[, ..names]
  sub_df_long<-melt(sub_df,id.vars=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),
                    measure=patterns(".cov", ".meth"),variable.factor=FALSE,variable.name="sample",
                    value.name=c("cov","meth"),na.rm=TRUE)
  sub_df_cond<-sub_df_long[,list(Nsamples=.N,mean_cov=mean(cov),min_cov=min(cov),max_cov=max(cov),mean_meth=mean(meth),
                                 min_meth=min(meth),max_meth=max(meth)),
                           by=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
  sub_df_cond_red<-sub_df_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
  sub_df_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
  
  #making sure that we don have same reference ids in test and train datasets
  if(isTest){
    print("filtering")
    idx<-!unlist(lapply(sub_df_cond_red$meta, is.element, set = train_ids))
    sub_df_cond_red<-sub_df_cond_red[idx,]
    N=1000
  }else{
    N=2000
  }
  
  ds <- create_datasets(sub_df_cond_red, ded_ref, strsplit(ids[[1]], "_")[[1]][3], TO_SAVE=!FALSE, N = N )
  return(ds)
  
}

data<-load_data(species[[1]])
train_sample_ids<-data$sample_ids
train_samples<-data$samples
train_data<-mclapply(data$sample_ids, get_dataset, meth_data_mean=data$meth_data_mean, 
                     cgg_names=data$cgg_names, ded_ref=data$ded_ref, uc_map=data$uc_map)
##getting train_test models for all the tissues:
train_models<-mclapply(seq_along(train_samples), function(i) {simpleCache(cacheName=paste0("methTrain_",train_samples[[i]]),
               instruction = {train_test(train_data[[i]]$x, train_data[[i]]$y, type=train_sample_ids[[i]][[1]], 
                ifRand=paste0("noRand", "Train"), run=0, k=kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
               cacheDir=paste0(subdir,"/RCache"),assignToVariable="res",recreate=FALSE)})

##training on random labels for comparison
train_models_rand<-mclapply(seq_along(train_samples), function(i) {simpleCache(cacheName=paste0("methTrainRand_",train_samples[[i]]),
                                                                  instruction = {train_test(train_data[[i]]$x, 
                                                                                sample(train_data[[i]]$y), 
                                                                                type=train_sample_ids[[i]][[1]], 
                                                                                ifRand=paste0("rand", "Train"), run=0, k=kmer, 
                                                                                subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
                                                                  cacheDir=paste0(subdir,"/RCache"),
                                                                  assignToVariable="resRand",recreate=FALSE)}) 

##getting the predictions in the same tissues
pred_in_other_tissues<-function(model, data, sample_ids, train_samples, I, ifRand ){
  #creating test datasets from the same species
  test_data_tisssues<-mclapply(sample_ids, get_dataset, meth_data_mean=data$meth_data_mean, ded_ref=data$ded_ref,
                               cgg_names=data$cgg_names,uc_map=data$uc_map, isTest=TRUE, train_ids=model$train_ids)
  test_data_tisssues_x<-lapply(test_data_tisssues, function(x) x$x)
  if(ifRand!="rand"){
  test_data_tisssues_y<-lapply(test_data_tisssues, function(x) x$y)
  }else{
    test_data_tisssues_y<-lapply(test_data_tisssues, function(x) sample(x$y))
  }
  rm(test_data_tisssues)
  
  if(NROW(test_data_tisssues_x)>0){
  #building predictions
  preddec_list<-mclapply(test_data_tisssues_x, predict, object=model$model[[1]])
  pref<-mapply(evaluatePrediction, preddec_list, test_data_tisssues_y, SIMPLIFY=FALSE,  MoreArgs=(list(allLabels=c(-1,1), 
                                                                                                       print=FALSE)))
  f1_list<-lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec<-mclapply(test_data_tisssues_x, predict, object=model$model[[1]], predictionType="decision")
  rocdata_list<-mapply(computeROCandAUC, preddec_list_dec, test_data_tisssues_y, MoreArgs=(list(allLabels=c(-1, 1))))
  
  roc_dt_list=lapply(seq_along(rocdata_list), function(i) data.table(fdr=unlist(rocdata_list[[i]]@FPR),
                        tpr=unlist(rocdata_list[[i]]@TPR),auc=unlist(rocdata_list[[i]]@AUC),f1=f1_list[[i]], 
                        run=0, type=sample_ids[[i]][[1]],ifRand=paste0(ifRand, "Test"), k=0, C=0,min_motif="NN",
                        max_motif="NN", test_tissue=train_samples[-I][[i]]))
  
  
  roc_dt_tissues<-rbindlist(roc_dt_list)
  roc_dt_tissues[,train_tissue:=train_samples[[I]]]
  }
  else{
    roc_dt_tissues<-data.frame()
  }
  roc_dt_train<-model$roc_dt
  roc_dt_train[,test_tissue:="trained on"]
  roc_dt_train[,train_tissue:=train_samples[[I]]]
  print("got predictions for all tissues within one sample")
  roc<-rbind(roc_dt_train, roc_dt_tissues)
  return(roc)
}
if(length(train_samples)>1){
roc_tissues_list<-mclapply(seq_along(train_samples), function(i) pred_in_other_tissues(train_models[[i]], data, 
                                                                train_sample_ids[-i],train_samples,i, "noRand"))

roc_tissues<-rbindlist(roc_tissues_list)


##prediction using random labels
roc_tissues_list_rand<-mclapply(seq_along(train_samples), function(i) pred_in_other_tissues(train_models_rand[[i]], 
                                                                        data, train_sample_ids[-i],train_samples,i, "rand"))
roc_tissues_rand<-rbindlist(roc_tissues_list_rand)
}else{
  roc_tissues<-train_models[[1]]$roc_dt
  roc_tissues[,test_tissue:="trained on"]
  roc_tissues[,train_tissue:=train_samples[[1]]]
  roc_tissues_rand<-train_models_rand[[1]]$roc_dt
  roc_tissues_rand[,test_tissue:="trained on"]
  roc_tissues_rand[,train_tissue:=train_samples[[1]]]
}
rm(data)

print(mem_used())
###uploading the other test species (the short-list)
test_species<-read.table("/scratch/lab_bock/dromanovskaia/CompEpi/compEpi/meta/selected_abbr.txt")$V1
test_species<-test_species[test_species!=species[[1]]]
if(isDev){
  test_species<-test_species[c(2:12)]
}


in_other_sample<-function(SP, train_samples, x){
  test_sample<-load_data(SP)
  
  ##getting lists of indexes to analyse
  test_ids<-test_sample$sample_ids[is.element(test_sample$samples, train_samples)]
  sel_model<-train_models[is.element(train_samples,test_sample$samples)]
  sel_model_rand<-train_models_rand[is.element(train_samples,test_sample$samples)]
  samples_test<-test_sample$samples[is.element(test_sample$samples, train_samples)]
  #create ds
  test_ds<-lapply(seq_along(test_ids),function(i) get_dataset(test_ids[[i]], meth_data_mean = test_sample$meth_data_mean,
                                                            cgg_names = test_sample$cgg_names, ded_ref = test_sample$ded_ref, 
                                                            uc_map=test_sample$uc_map, isTest=TRUE, 
                                                            train_ids=train_models[[i]]$train_ids))
  test_data_sample_x<-lapply(test_ds, function(x) x$x)
  test_data_sample_y<-lapply(test_ds, function(x) x$y)
  test_data_sample_y_rand<-lapply(test_ds, function(x) sample(x$y))
  rm(test_ds)
  
  #building predictions
  preddec_list<-lapply(seq_along(test_data_sample_x), function(x) predict(test_data_sample_x[[x]], 
                                                                          object=train_models[[x]]$model[[1]]))
  pref<-mapply(evaluatePrediction, preddec_list, test_data_sample_y, SIMPLIFY=FALSE,  MoreArgs=(list(allLabels=c(-1,1), print=FALSE)))
  f1_list<-lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  preddec_list_dec<-lapply(seq_along(test_data_sample_x), function(x) predict(test_data_sample_x[[x]], 
                                                                    object=train_models[[x]]$model[[1]], predictionType="decision"))
  rocdata_list<-mapply(computeROCandAUC, preddec_list_dec, test_data_sample_y, MoreArgs=(list(allLabels=c(-1,1))))
  roc_dt_list=lapply(seq_along(rocdata_list), function(i) data.table(fdr=unlist(rocdata_list[[i]]@FPR),
                                                                    tpr=unlist(rocdata_list[[i]]@TPR),
                                                                    auc=unlist(rocdata_list[[i]]@AUC),f1=f1_list[[i]], 
                                                                    run=x, type=test_ids[[i]][[1]],ifRand="noRandTest", k=0, C=0,
                                                                    min_motif="NN",max_motif="NN", train_tissue=samples_test[[i]],
                                                                    test_tissue=samples_test[[i]]))
  roc_dt_sample<-rbindlist(roc_dt_list)
  
  ##prediction the random model on random labels
  #building predictions
  pref_rand<-mapply(evaluatePrediction, preddec_list, test_data_sample_y_rand, SIMPLIFY=FALSE,  
                    MoreArgs=(list(allLabels=c(-1,1), print=FALSE)))
  f1_list<-lapply(pref_rand, function(x) {get_f1(x$SENS, x$PREC)})
  
  rocdata_list_rand<-mapply(computeROCandAUC, preddec_list_dec, test_data_sample_y_rand, MoreArgs=(list(allLabels=c(-1,1))))
  roc_dt_list_rand=lapply(seq_along(rocdata_list_rand), function(i) data.table(fdr=unlist(rocdata_list_rand[[i]]@FPR),
                                                                      tpr=unlist(rocdata_list_rand[[i]]@TPR),
                                                                      auc=unlist(rocdata_list_rand[[i]]@AUC),f1=f1_list[[i]], 
                                                                      run=x, type=test_ids[[i]][[1]],ifRand="randTest", k=0, C=0,
                                                                      min_motif="NN",max_motif="NN", train_tissue=samples_test[[i]], 
                                                                      test_tissue=samples_test[[i]]))
  roc_dt_sample_rand<-rbindlist(roc_dt_list_rand)
  rm(test_sample)
  return(list(pred=roc_dt_sample, pred_rand=roc_dt_sample_rand))
}


res_list<-mclapply(seq_along(test_species), function(x) in_other_sample(test_species[[x]], train_samples, x))

res<-rbindlist(lapply(res_list, function(x) x$pred))
res_rand<-rbindlist(lapply(res_list, function(x) x$pred_rand))
print(mem_used())
write.table(rbind(res, res_rand, roc_tissues, roc_tissues_rand), paste0(subdir, "/all_points.csv"),
            quote=FALSE,sep="\t",row.names=FALSE)

full<-rbind(res, res_rand, roc_tissues, roc_tissues_rand)


train<-full[full$test_tissue=="trained on"]

test<-full[full$test_tissue!="trained on"]
if(NROW(test)>0){
test[,species:=tstrsplit(type, "_", fixed=TRUE)[[1]]]
}
n_pages<-ceiling(length(unique(test$species))/8)
ncol_my=length(train_samples)
library(ggforce)

pdf(paste0(subdir,"/ROC.pdf"),height=10,width=8)

for (i in c(1:n_pages)){
  print(ggplot(test, aes(x=fdr, y=tpr, col=test_tissue, alpha=ifRand)) + geom_line(aes(group=ifRand)) + 
          facet_grid_paginate(species~train_tissue, nrow=8,ncol=ncol_my,page=i)+geom_line(data=train, aes(group=ifRand)) +
          scale_color_manual(values = color_tissues)+
          scale_alpha_manual(values = c("noRandTest"=1, "randTest"=0.3, "noRandTrain"=1, "randTrain" =0.5)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)))
}
dev.off()

library(dplyr)
auc_res<-test %>%
  group_by(.dots=c("type", "train_tissue", "species", "ifRand")) %>%
  summarize(auc=mean(auc), f1=mean(f1), test_tissue=first(test_tissue))

auc_res_train<-train %>%
  group_by(.dots=c("type", "train_tissue", "ifRand")) %>%
  summarize(auc=mean(auc), f1=mean(f1), test_tissue=first(test_tissue), species=species[[1]])

write.table(rbind(auc_res_train, auc_res), paste0(subdir, "/all_aucs.csv"),quote=FALSE,sep="\t",row.names=FALSE)

print("saved stuff!")
