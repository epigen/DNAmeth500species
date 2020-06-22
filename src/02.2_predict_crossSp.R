#!/usr/bin/env Rscript
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/02.0_predict_initialization.R"))


kmer = NA
if (is.na(kmer)){kmer = c(1:10)}
args = commandArgs(trailingOnly = TRUE)
species = args[1]


sessionInfo()


#subdir_name = args[2]
#kmer = as.numeric(args[3])


test_species <- unique(stats_annot$species)
test_species <- test_species[! test_species%in% c(species)]

if(isDev){
  test_species <- test_species[c(200,10)]
}

subdir_name = "screen"

wd = file.path(processed_dir)
setwd(wd)

mywd = "/scratch/lab_bock/shared/projects/compEpi/results_analysis/02_predict_meth"
subdir = file.path(mywd,"02.2_test_on_other_species",subdir_name,species)

dir.create(subdir, recursive = TRUE)
dir.create(file.path(subdir, "sequences"), recursive = TRUE)
dir.create(file.path(subdir, "stats"), recursive = TRUE)
library(ggforce)

#load data function
read_data  <- function(SP, isTest = FALSE, train_ids = NA){
  print(paste0("reading data about ", SP))
  meth_data_mean = fread(list.files(path = paste0(SP, "/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/"),
                                    pattern = ".*_mean_meth.tsv", full.names = TRUE))
  uc_map = fread(list.files(path = paste0(SP, "/toSelf_filtered_0.08mm_final_concat"),
                            pattern = ".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",recursive = TRUE,full.names = TRUE))
  ded_ref = readDNAStringSet(paste0(SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  ded_ref_cgg = as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names = names(ded_ref[ded_ref_cgg]) 
  
  #select features (methylated/unmethylated)
  meth_data_mean_long = melt(meth_data_mean,id.vars = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),
                             measure = patterns(".cov", ".meth"),variable.factor = FALSE,variable.name = "sample",
                             value.name = c("cov","meth"),na.rm = TRUE)
  meth_data_mean_cond = meth_data_mean_long[,list(Nsamples = .N,mean_cov = mean(cov),min_cov = min(cov),max_cov = max(cov),
                                                  mean_meth = mean(meth),min_meth = min(meth),max_meth = max(meth)),
                                            by = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
  meth_data_mean_cond_red = meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
  meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
  
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
  if (length(meth_data_mean_cond_red[meth_data_mean_cond_red$category == 1]$meta) > 0 & length(meth_data_mean_cond_red[meth_data_mean_cond_red$category == -1]$meta) > 0){
    a <- create_datasets(meth_data_mean_cond_red, ded_ref, SP, FALSE, N)
    print("dataset formed")
    return(a)
    
  } else{
    print("with current thresholds can not create a test dataset")
    return (list(x = 0, y = 0))
  }
}

# data from the main species, we build predictive model on
train_data <- read_data(species)

#building a standart prediction:
simpleCache(cacheName = "methTrain_noRand",instruction = {train_test(x = train_data$x, y = train_data$y,
              type = species[[1]], ifRand = "noRandTrain", run = 0, k = kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
            cacheDir = paste0(subdir,"/RCache"), assignToVariable = "res", recreate = FALSE)

roc_res_train = res$roc_dt
model = res$model[[1]]
train_ids = res$train_ids

#loading all other datasets
# data from the species, we will test on
test_data <- lapply(test_species, read_data, isTest = TRUE, train_ids = train_ids)

## delete the empty datasets
a <- lapply(test_data, function(x) (length(x$y)==1))

if (length(which(unlist(a)))>0){
  test_data <- test_data[-which(unlist(a))]
  test_species <- test_species[-which(unlist(a))]
}

test_data_x <- lapply(test_data, function(x) x$x)
test_data_y <- lapply(test_data, function(x) x$y)
rm(test_data)


test_on_other <- function(fit,test_data_x, test_data_y, run, ifRand, kmer, subdir){
  preddec_list <- mclapply(test_data_x, predict, object = fit)
  
  pref <- mapply(evaluatePrediction, preddec_list, test_data_y, SIMPLIFY = FALSE,  
                 MoreArgs = (list(allLabels = unique(test_data_y[[1]]), print = FALSE)))
  f1_list <- lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec <- mclapply(test_data_x, predict, object = fit, predictionType = "decision")
  rocdata_list <- mapply(computeROCandAUC, preddec_list_dec, test_data_y, 
                         MoreArgs = (list(allLabels = unique(test_data_y[[1]]))))
  
  roc_dt_list <- lapply(seq_along(rocdata_list), function(i) data.table(fdr = unlist(rocdata_list[[i]]@FPR),
                      tpr = unlist(rocdata_list[[i]]@TPR),auc = unlist(rocdata_list[[i]]@AUC),f1 = f1_list[[i]], 
                      run = run, type = test_species[[i]],ifRand = paste0(ifRand, "Test"), k = 0, C = 0,
                      min_motif = "NN",max_motif = "NN"))
  
  roc_dt_full <- rbindlist(roc_dt_list)
  
  return(list(roc_dt = roc_dt_full))
}

#predicting on other species
simpleCache(cacheName = "methPred_noRand_uc",instruction = {test_on_other(fit = model, test_data_x = test_data_x, 
            test_data_y = test_data_y, ifRand = "noRand",run = 0,subdir = subdir,k = kmer)},
            cacheDir = paste0(subdir,"/RCache"),assignToVariable = "resTest",recreate = FALSE)



#random labels 1 times for the train dataset
rand_labs <- lapply(seq_len(1), function(x) {return(list(name = x,y_rand = sample(train_data$y)))})
rand_res_train <- mclapply(seq_along(rand_labs), function(x) 
  {simpleCache(cacheName = paste0("methTrain_Rand",rand_labs[[x]]$name), 
               instruction = {train_test(x = train_data$x, y = rand_labs[[x]]$y_rand, type = species[[1]], 
                          ifRand = "RandTrain", run = x, k = kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},
               cacheDir = paste0(subdir,"/RCache"),assignToVariable = "resRand",recreate = FALSE)})

rand_labels_test <- lapply(seq_len(1), function(y){lapply(test_data_y, function(y) y_rand = sample(y))})



rand_res_test <- mclapply(seq_len(1), function(x) {simpleCache(cacheName = paste0("methPred_Rand", x),
                      instruction = {test_on_other(fit = rand_res_train[[x]]$model[[1]], 
                      test_data_x = test_data_x, test_data_y = rand_labels_test[[x]], ifRand = "Rand", run = x, 
                      subdir = subdir, k = kmer)},cacheDir = paste0(subdir,"/RCache"),
                      assignToVariable = "resRandTest",recreate = FALSE)})

rand_res_test_list <- lapply(rand_res_test, function(x){return (x$roc_dt)})
rand_res_test <- rbindlist(rand_res_test_list)
res_test <- rbind(resTest$roc_dt, rand_res_test)

rand_res_train_list <- lapply(rand_res_train, function(x){return (x$roc_dt)})
rand_res_train_long <- rbindlist(rand_res_train_list)
res_train <- rbind(res$roc_dt, rand_res_train_long)


auc_res <- res_test%>%
  group_by(.dots = c("ifRand", "type")) %>%
  summarize(auc = mean(auc), f1 = mean(f1))

auc_res_train <- res_train %>%
  group_by(.dots = c("ifRand", "type")) %>%
  summarize(auc = mean(auc), f1 = mean(f1))

tit <- paste0(species[[1]], ": ",paste(round(auc_res_train[auc_res_train$ifRand=="noRandTrain", ]$auc,3), 
                                       round(auc_res_train[auc_res_train$ifRand=="RandTrain", ]$auc,3), sep = "("), ")")

write.table(union(auc_res, auc_res_train), paste0(subdir, "/all_aucs.csv"),quote = FALSE,sep = "\t",row.names = FALSE)


auc_res <- as.data.table(auc_res)
auc_res[,x:=0.75,]
auc_res[,y:=ifelse(ifRand=="RandTest",0.05,0.2),]


res_train <- select(res_train, -c(type))
n_pages <- ceiling(length(unique(res_test$type))/24)

pdf(paste0(subdir,"/ROC.pdf"),height = 10,width = 8)
for (i in c(1:n_pages)){
  print(ggplot(res_test, aes(x = fdr,y = tpr, col = ifRand))+geom_line(aes(group = run))+
          facet_wrap_paginate(~type, ncol = 4, nrow = 6, page = i)+
          geom_line(data = res_train, aes(group = run, col = ifRand))+ggtitle(tit)+
          scale_color_manual(values = c("RandTest" = "lightgrey","noRandTest" = "blue", "noRandTrain" = "red", 
                                        "RandTrain" = "darkgrey"))+
          geom_text(data = auc_res,aes(x = x,y = y,label = paste0("auc=",round(auc,3)))))+ 
          theme(panel.spacing.x = unit(0.5, "lines"),panel.spacing.y = unit(1, "lines"))
}
dev.off()

