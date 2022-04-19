#!/usr/bin/env Rscript
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/05.0_predict_initialization.R"))

library(ggforce)

isDev = FALSE

args = commandArgs(trailingOnly = TRUE)
species = args[1]


#subdir_name = args[2]
#kmer = as.numeric(args[3])

kmer = NA
if (is.na(kmer)){kmer = c(1:10)}

test_species <- unique(stats_annot$species)
test_species <- test_species[! test_species%in% c(species)]

if(isDev){
  test_species <- test_species[c(200,10)]
}

subdir_name = "screen"

wd = file.path(processed_dir)
setwd(wd)

mywd = file.path(analysis_dir,"05_predict_meth")
subdir = file.path(mywd,"05.2_test_on_other_species",subdir_name,species)
model_dir = file.path(mywd,"05.1_within_species",subdir_name,species)

dir.create(subdir, recursive = TRUE)


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

#loading the standart prediction model:
simpleCache(cacheName = "methPred_noRand_uc",instruction = {train_test(x = train_data$x, y = train_data$y, type = species, ifRand = "noRandTrain", run = 0, k = kmer, subdir = subdir, SAVE_TRAIN_IDS = TRUE)},cacheDir = paste0(model_dir,"/RCache"), assignToVariable = "res", recreate = FALSE)

roc_res_train = res$roc_dt
model = res$model[[1]]
#train_ids = res$train_ids

#loading all other datasets
# data from the species, we will test on
test_data <- lapply(test_species, read_data, isTest = TRUE)

## delete the empty datasets
a <- lapply(test_data, function(x) (length(x$y)==1))

if (length(which(unlist(a)))>0){
  test_data <- test_data[-which(unlist(a))]
  test_species <- test_species[-which(unlist(a))]
}

test_data_x <- lapply(test_data, function(x) x$x)
test_data_y <- lapply(test_data, function(x) x$y)
rm(test_data)


test_on_other <- function(fit,test_data_x, test_data_y, run, ifRand, kmer, subdir, test_sp = test_species){
  preddec_list <- mclapply(test_data_x, predict, object = fit)
  
  pref <- mapply(evaluatePrediction, preddec_list, test_data_y, SIMPLIFY = FALSE,  
                 MoreArgs = (list(allLabels = unique(test_data_y[[1]]), print = FALSE)))
  f1_list <- lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec <- mclapply(test_data_x, predict, object = fit, predictionType = "decision")
  rocdata_list <- mapply(computeROCandAUC, preddec_list_dec, test_data_y, 
                         MoreArgs = (list(allLabels = unique(test_data_y[[1]]))))
  
  roc_dt_full= data.table()
    
   for(i in 1:length(test_sp)){
       print(i)
       df = data.table(fdr = rocdata_list[[i]]@FPR, tpr = rocdata_list[[i]]@TPR,
                      auc=rocdata_list[[i]]@AUC, run = run, 
                        type = test_sp[[i]], ifRand = ifRand)
       roc_dt_full = rbind(roc_dt_full,df)
   }
#  roc_dt_list <- lapply(seq_along(rocdata_list), function(i) data.table(fdr = unlist(rocdata_list[[i]]@FPR),
#                      tpr = unlist(rocdata_list[[i]]@TPR),
#                      auc = unlist(rocdata_list[[i]]@AUC),f1 = f1_list[[i]], 
#                      run = run, type = test_species[[i]],ifRand = paste0(ifRand, "Test"), k = #0, C = 0,
 #                     min_motif = "NN",max_motif = "NN"))
#  
 # roc_dt_full <- rbindlist(roc_dt_list)
  
  return(list(roc_dt = roc_dt_full))
}
                      
                      
print("testing on the model")
#predicting on other species
simpleCache(cacheName = "methPred_noRand",instruction = {test_on_other(fit = model, 
              test_data_x = test_data_x, 
            test_data_y = test_data_y, 
            ifRand = "noRandTest",run = 0,subdir = subdir,k = kmer)},
            cacheDir = paste0(subdir,"/RCache"),assignToVariable = "resTest",recreate = FALSE)



#random labels 1 times for the train dataset
rand_labs <- sample(train_data$y)

##loading the random model: 
simpleCache(cacheName = "methPred_Rand_uc_1", 
               instruction = {train_test(x = train_data$x, y = rand_labs[[x]]$y_rand, 
                                         type = species[[1]], 
                                          ifRand = "RandTrain", 
                                           run = 1, k = kmer, 
                                         subdir = subdir)},
               cacheDir = paste0(model_dir,"/RCache"),
            assignToVariable = "resRand",recreate = FALSE)
                        
                        
                        
rand_model =  resRand$model[[1]]

 simpleCache(cacheName = "methPred_Rand_1",
                      instruction = {test_on_other(fit = rand_model, 
                      test_data_x = test_data_x, 
                       test_data_y = lapply(test_data_y, function(x) sample(x)), 
                        ifRand = "RandTest", run = 1)},cacheDir = paste0(subdir,"/RCache"),
                      assignToVariable = "resRandTest",recreate = FALSE)

#rand_res_test_list <- lapply(rand_res_test, function(x){return (x$roc_dt)})
#rand_res_test <- rbindlist(rand_res_test_list)

res_test <- rbind(resTest$roc_dt, resRandTest$roc_dt)

#rand_res_train_list <- lapply(rand_res_train, function(x){return (x$roc_dt)})
#rand_res_train_long <- rbindlist(rand_res_train_list)
res_train <- rbind(res$roc_dt, resRand$roc_dt)

auc_res <-unique(res_test[, c("ifRand", "type", "auc")])
                                    
auc_res_train <- unique(res_train[, c("ifRand", "type", "auc")])
                                     

tit <- paste0(species, ": ",paste(round(auc_res_train[ifRand=="noRand", ]$auc,3), round(auc_res_train[auc_res_train$ifRand=="rand", ]$auc,3), sep = " ("), ")")

write.table(union(auc_res, auc_res_train), paste0(subdir, "/all_aucs.csv"),quote = FALSE,sep = "\t",row.names = FALSE)


auc_res <- setDT(auc_res)
auc_res[,x:=0.75,]
auc_res[,y:=ifelse(ifRand=="RandTest",0.05,0.1),]

res_train <- select(res_train, -c(type))
                                     
n_pages <- ceiling(length(unique(res_test$type))/24)

pdf(paste0(subdir,"/ROC.pdf"),height = 10,width = 8)
for (i in c(1:n_pages)){
  print(ggplot(res_test, aes(x = fdr,y = tpr, col = ifRand))+geom_line(aes(group = run))+
          facet_wrap_paginate(~type, ncol = 4, nrow = 6, page = i) +
          geom_line(data = res_train, aes(group = run, col = ifRand)) + ggtitle(tit)+
          scale_color_manual(values = c("RandTest" = "lightgrey","noRandTest" = "blue", "noRand" = "red",  "rand" = "darkgrey"))+
          geom_text(data = auc_res,aes(x = x,y = y,label = paste0("auc=",round(auc,3)))))+ 
          theme(panel.spacing.x = unit(0.5, "lines"),panel.spacing.y = unit(1, "lines"))
}
dev.off()
