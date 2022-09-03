#!/usr/bin/env Rscript
library(pryr)
kmer=NA

args = commandArgs(trailingOnly=TRUE)
species=args[1]
if (is.na(kmer)){kmer=c(1:10)}

subdir_name="bootstrap"

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/05.0_predict_initialization.R"))

wd=file.path(processed_dir)
setwd(wd)

subdir=file.path(analysis_dir,"06_inv", "06.2_verify_inverted",subdir_name,species)

dir.create(subdir,recursive=TRUE)
dir.create(file.path(subdir, "sequences"), recursive=TRUE)
dir.create(file.path(subdir, "stats"), recursive=TRUE)

read_data  <- function(SP, isTest = FALSE, train_ids = NA, SEED = TRUE){
  print(paste0("reading data about ", SP))
  meth_data_mean = fread(list.files(path = paste0(SP, "/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/"),
                                    pattern = ".*_mean_meth.tsv",
                                    full.names = TRUE))
  
  uc_map = fread(list.files(path = paste0(SP, "/toSelf_filtered_0.08mm_final_concat"),
                            pattern = ".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",
                            recursive = TRUE,full.names = TRUE))
  
  ded_ref = readDNAStringSet(paste0(SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  ded_ref_cgg = as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names = names(ded_ref[ded_ref_cgg]) 
  
  #select features (methylated/unmethylated)
  meth_data_mean_long = melt(meth_data_mean,id.vars = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),
                             measure = patterns(".cov", ".meth"),variable.factor = FALSE,variable.name = "sample",
                             value.name = c("cov","meth"),na.rm = TRUE)
  meth_data_mean_cond = meth_data_mean_long[,list(Nsamples = .N,mean_cov = mean(cov),
                                                  min_cov = min(cov),max_cov = max(cov),
                                                  mean_meth = mean(meth),min_meth = min(meth),
                                                  max_meth = max(meth)),
                                            by = c("meta","score","consN","mean_diffNts",
                                                   "max_diffNts","min_diffNts")]
  meth_data_mean_cond_red = meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
  meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
  print(NROW(meth_data_mean_cond_red[meth_data_mean_cond_red$category ==1, ]))
  N = 2000
  
  if(isTest){
    print("filtering")
    idx <- !unlist(lapply(meth_data_mean_cond_red$meta, is.element, set = train_ids))
    meth_data_cond_red <- meth_data_mean_cond_red[idx,]
    N = 1000
  }
  
  print("features mapped")
  #minimum features of either group --> maximum number of features for each group
  #  if (length(meth_data_mean_cond_red$meta) > 100){
  if (length(meth_data_mean_cond_red[meth_data_mean_cond_red$category == 1]$meta) > 0 & 
      length(meth_data_mean_cond_red[meth_data_mean_cond_red$category == -1]$meta) > 0){
    a <- create_datasets(meth_data_mean_cond_red, ded_ref, SP, FALSE, N, SEED)
    print("dataset formed")
    return(a)
    
  } else{
    print("with current thresholds can not create a test dataset")
    return (0)
  }
}

train_data <- sapply(c(1:10), function(x) read_data(species, SEED = F), simplify = F)



train_models <- mclapply(seq_along(train_data), function(i) {simpleCache(cacheName=paste0("methTrain_",i),
                instruction = {train_test(train_data[[i]]$x, train_data[[i]]$y, type=i,
                ifRand=paste0("noRand", "Train"), run=0, k=kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
                cacheDir=paste0(subdir,"/RCache"), assignToVariable="res", recreate=FALSE)})


train_models_df <- rbindlist(lapply(train_models, function(x) x$roc_dt))
train_models_df <- unique(train_models_df[, c("type", "auc", "run", "ifRand")])
                                    
train_models_rand <- mclapply(seq_along(train_data), function(i) {simpleCache(cacheName=paste0("methTrainRand_",i),
                     instruction = {train_test(train_data[[i]]$x, sample(train_data[[i]]$y), type=i,
                    ifRand=paste0("rand", "Train"), run=0, k=kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
                    cacheDir=paste0(subdir,"/RCache"), assignToVariable="resRand", recreate=FALSE)}) 

train_models_df_rand <- rbindlist(lapply(train_models_rand, function(x) x$roc_dt))
train_models_df_rand <- unique(train_models_df_rand[, c("type", "auc", "run", "ifRand")])

df <- rbind(train_models_df, train_models_df_rand)                                         
my_wt(df, file.path(subdir, "bootstrap_AUC.tsv"))
