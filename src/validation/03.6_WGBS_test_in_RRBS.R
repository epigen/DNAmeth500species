source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/05.0_predict_initialization.R"))

library(kebabs)
library(caret)
library(ROCR)

args = commandArgs(trailingOnly=TRUE)

species=args[[1]]

path_to_folder <- file.path(analysis_dir, "validation", "03_WGBS", "03.6_WGBS_test_in_RRBS", species)
dir.create(path_to_folder, recursive = T)
setwd(path_to_folder)
                        
#load data function
read_data  <- function(SP){
  print(paste0("reading data about ", SP))

    
  meth_data_mean = fread(list.files(path = paste0(processed_dir,"/", SP, "/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/"),pattern = ".*_mean_meth.tsv", full.names = TRUE))
    
  uc_map = fread(list.files(path = paste0(processed_dir, "/",SP, "/toSelf_filtered_0.08mm_final_concat"),
                            pattern = ".*_uc.all.aligned.bsmap.0.08.mism.1.r.cov$",recursive = TRUE,full.names = TRUE))
  ded_ref = readDNAStringSet(paste0(processed_dir,"/", SP,"/reduced/consensus//toSelf_filtered_0.08mm_final.fa"))
  ded_ref_cgg = as.character(substring(ded_ref,1,3))=="CGG"
  cgg_names = names(ded_ref[ded_ref_cgg]) 
  
  #select features (methylated/unmethylated)
  meth_data_mean_long = melt(meth_data_mean,id.vars = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"), measure = patterns(".cov", ".meth"),variable.factor = FALSE,variable.name = "sample", value.name = c("cov","meth"),na.rm = TRUE)
  
    meth_data_mean_cond = meth_data_mean_long[,list(Nsamples = .N,mean_cov = mean(cov),min_cov = min(cov),max_cov = max(cov), mean_meth = mean(meth),min_meth = min(meth),max_meth = max(meth)),by = c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
    
    meth_data_mean_cond_red = meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000&(min_meth>80|max_meth<20)&meta%in%uc_map[V13==1]$V4&meta%in%cgg_names]
    
    meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
  
  N = 1000
    
  print("features mapped")
  #minimum features of either group --> maximum number of features for each group
#  if (length(meth_data_mean_cond_red$meta) > 100){
  if (min(table(meth_data_mean_cond_red$category)) > 10){
    a <- create_datasets(meth_data_mean_cond_red, ded_ref, SP, FALSE, N)
    print("dataset formed")
    return(a) 
  }else{
    print("with current thresholds can not create a test dataset")
    return (list(x = 0, y = 0))
  }
}

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
                        
#path_to_data <- file.path(analysis_dir, "validation", "03_WGBS", "03.3_fragments")

simpleCache(cacheName=paste0("methPred_noRand_uc_1234"), instruction={ train_test(x_train=split_ds$x_train,
            x_test = split_ds$x_test, y_train = split_ds$y_train, y_test = split_ds$y_test,
                                                        ifRand='noRand', k=kmer, runid = 0)},
            cacheDir=paste0("../../03.4_prediction/Branchiostoma_lanceolatum/kebabs_model/RCache"), assignToVariable="res", recreate=FALSE)

t <- read_data("AB")
## load data to train on:
test_data <- lapply(sp_df$species, read_data)

## delete the empty datasets
a <- lapply(test_data, function(x) (length(x$y)==1))

if (length(which(unlist(a)))>0){
  test_data <- test_data[-which(unlist(a))]
  test_species <- test_species[-which(unlist(a))]
}

test_data_x <- lapply(test_data, function(x) x$x)
test_data_y <- lapply(test_data, function(x) x$y)
rm(test_data)
                      
## predicting
                      
#predicting on other species
simpleCache(cacheName = "methPred_noRand_uc",instruction = {test_on_other(fit = model, test_data_x = test_data_x, 
            test_data_y = test_data_y, ifRand = "noRand",run = 0,subdir = subdir,k = kmer)},
            cacheDir = paste0(subdir,"/RCache"),assignToVariable = "resTest",recreate = FALSE)
                      
my_wt(resTest$roc_dt, "roc_dt.tsv")
                      
unique_auc <- unique(resTest$roc_dt[, c("type", "auc", "ifRand")])
unique_auc$AUC_train <- res$roc_dt$auc[[1]]
unique_auc$species_train <- species                     
my_wt(unique_auc, "auc_test.tsv")