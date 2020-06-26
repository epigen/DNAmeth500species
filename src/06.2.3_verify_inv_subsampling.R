#!/usr/bin/env Rscript
library(pryr)

args = commandArgs(trailingOnly=TRUE)
species=args[1]
if (is.na(kmer)){kmer=c(1:10)}

subdir_name="bootstrap"

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/02.0_predict_initialization.R"))


wd=file.path(processed_dir)
setwd(wd)

mywd=file.path(analysis_dir, "02_predict_meth")

subdir=file.path(mywd,"02.4_verify_inverted",subdir_name,species)

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




train_models_rand <- mclapply(seq_along(train_data), function(i) {simpleCache(cacheName=paste0("methTrainRand_",i),
                     instruction = {train_test(train_data[[i]]$x, sample(train_data[[i]]$y), type=i,
                    ifRand=paste0("rand", "Train"), run=0, k=kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
                    cacheDir=paste0(subdir,"/RCache"), assignToVariable="resRand", recreate=FALSE)}) 



test_species <- read.table(file.path(Sys.getenv("CODEBASE"),
                                     "DNAmeth500species/meta/species_list.txt")$V1)


train_ids <- unique(sapply(train_models, function(x) x$train_ids))

test_data <- lapply(test_species, read_data, isTest = TRUE, train_ids = train_ids)
test_data <- test_data[sapply(test_data, length) > 1]

test_data_x <- lapply(test_data, function(x) x$x)
test_data_y <- lapply(test_data, function(x) x$y)
rm(test_data)

test_on_other <- function(fit,test_data_x, test_data_y, run, ifRand, kmer, subdir){
  print(run)
  preddec_list <- mclapply(test_data_x, predict, object = fit$model[[1]])
  
  pref <- mapply(evaluatePrediction, preddec_list, test_data_y, SIMPLIFY = FALSE,  
                 MoreArgs = (list(allLabels = unique(test_data_y[[1]]), print = FALSE)))
  f1_list <- lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec <- mclapply(test_data_x, predict, object = fit$model[[1]], predictionType = "decision")
  rocdata_list <- mapply(computeROCandAUC, preddec_list_dec, test_data_y, 
                         MoreArgs = (list(allLabels = unique(test_data_y[[1]]))))
  
  roc_dt_list <- lapply(seq_along(rocdata_list), function(i) data.table(fdr = unlist(rocdata_list[[i]]@FPR),
                  tpr = unlist(rocdata_list[[i]]@TPR),auc = unlist(rocdata_list[[i]]@AUC),f1 = f1_list[[i]], 
                  run = run, type = test_species[[i]],ifRand = paste0(ifRand, "Test"), k = 0, C = 0,
                                                                        min_motif = "NN",max_motif = "NN"))
  
  roc_dt_full <- rbindlist(roc_dt_list)
  
  return(list(roc_dt = roc_dt_full))
}

res_list <- sapply(seq_along(train_models), function(i) test_on_other(train_models[[i]],
                                                                      test_data_x = test_data_x, 
                    test_data_y = test_data_y, ifRand = "noRand",run = i, subdir = subdir,k = kmer), 
                          simplify = F)

res <- rbindlist(lapply(res_list, function(x) x$roc_dt ))

test_data_y_rand <- lapply(test_data_y, function(y) y_rand = sample(y))
res_rand_list <- sapply(seq_along(train_models_rand), function(i) test_on_other(train_models_rand[[i]],
                      test_data_x = test_data_x, test_data_y = test_data_y_rand, ifRand = "Rand",
                      run = i, subdir = subdir, k = kmer), simplify = F)
res_rand <- rbindlist(lapply(res_rand_list, function(x) x$roc_dt ))
res_full <- rbind(res, res_rand)
write.table(res_full, file.path(subdir, "all_data.csv"), row.names = F, quote = F)

res_sum <- res_full %>%
  dplyr::group_by(type, ifRand, run) %>% 
  dplyr::summarise(AUC = first(auc), F1 = first(f1))

self_auc <- res_sum[res_sum$type==species & res_sum$ifRand=="noRandTest", ]

res_sum$color_class <- unlist(sapply(res_sum$type, function(x) sp_df[x, ]$color_class))

ggplot(res_sum, aes(x = type, y = AUC, fill = ifRand)) + geom_boxplot() +
  facet_wrap(~color_class, ncol = 1) + xlab("") + 
  scale_fill_manual(values = c("noRandTest"="#4682B4", "RandTest" = "grey")) + 
  geom_hline(data =  self_auc, mapping = aes(yintercept = AUC), linetype = "dashed", color = "red")

ggsave(file.path(subdir,  "AUC_in_other_sp.pdf"), width = 20, height = 8)
