source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/05.0_predict_initialization.R"))

library(kebabs)
library(caret)
library(ROCR)

args = commandArgs(trailingOnly=TRUE)

species=args[[1]]

read_data <- function(species_name, genomeid, cov_min, N = 1000){
    print(species_name)
    ##loading data
    meth_data_mean = fread(file.path(path_to_data,species_name, "mean_meth_per_fragment.tsv"))
  
    ##formatting
    meth_data_mean_long <- melt(meth_data_mean, measure=patterns("cov", "mean_meth"),
                            variable.factor=TRUE,
                         variable.name="sample",
                            value.name=c("cov","meth"),na.rm=TRUE)
                            
    meth_data_mean_cond <- meth_data_mean_long[,list(Nsamples=.N, mean_cov=mean(cov), min_cov=min(cov), max_cov=max(cov), mean_meth=mean(meth), min_meth=min(meth), max_meth=max(meth)),
                                        by=c("name")]
    
    meth_data_mean_cond_red=meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>cov_min & mean_cov<1000 &(min_meth>80|max_meth<20)]
    
    print(NROW(meth_data_mean_cond_red))

    meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]

    ##loading sequence
    ded_ref=readDNAStringSet(file.path(path_to_data, species_name,paste0(genomeid, "_fragments.fa")))

    N = 1000
    colnames(meth_data_mean_cond_red)[1] <- "meta"  
    if(min(table(meth_data_mean_cond_red$category)) > 10){
        a <- create_datasets(meth_data_mean_cond_red, ded_ref, species_name, FALSE, N)
        return(a)
    }else{
        print("not enough data")
        return(list(x = 0, y = 0))
    }
}

test_on_other <- function(fit,test_data_x, test_data_y, ifRand, test_species){
  preddec_list <- mclapply(test_data_x, predict, object = fit)
  
  pref <- mapply(evaluatePrediction, preddec_list, test_data_y, SIMPLIFY = FALSE,  
                 MoreArgs = (list(allLabels = unique(test_data_y[[1]]), print = FALSE)))
    
  #f1_list <- lapply(pref, function(x) {get_f1(x$SENS, x$PREC)})
  
  preddec_list_dec <- mclapply(test_data_x, predict, object = fit, predictionType = "decision")
  rocdata_list <- mapply(computeROCandAUC, preddec_list_dec, test_data_y, 
                         MoreArgs = (list(allLabels = unique(test_data_y[[1]]))))
  
  roc_dt_list <- lapply(seq_along(rocdata_list), function(i) data.table(fdr = unlist(rocdata_list[[i]]@FPR),
                      tpr = unlist(rocdata_list[[i]]@TPR),
                      auc = unlist(rocdata_list[[i]]@AUC),
                      type = test_species[[i]], ifRand = paste0(ifRand, "Test")))
  
  roc_dt_full <- rbindlist(roc_dt_list)
  
  return(list(roc_dt = roc_dt_full))
}

path_to_folder <- file.path(analysis_dir, "validation", "03_WGBS", "03.5_RRBS_test_in_WGBS", species)
dir.create(path_to_folder, recursive = T)
setwd(path_to_folder)
                        

                        
path_to_data <- file.path(analysis_dir, "validation", "03_WGBS", "03.3_fragments")
                        
                        
## loading the train model
simpleCache(cacheName="methPred_noRand_uc",instruction={train_test(x=x,y=y,type=species, 
                                                        ifRand='noRand',run=0,k=kmer, subdir=subdir)},
            cacheDir=paste0(file.path(analysis_dir, "05_predict_meth","05.1_within_species","screen",species,"RCache")),
            assignToVariable="res",recreate=FALSE)



##loading test data
annot_WGBS <- fread( file.path(analysis_dir, "validation", "03_WGBS", "WGBS_prediction_selection.tsv"))

test_data <- lapply(1:NROW(annot_WGBS), function(x) read_data(annot_WGBS$Species[x], annot_WGBS$genomeId[x], annot_WGBS$thr[x]))
#saveRDS(object = test_data, file = "test_data.RDS")
#print(test_data)
test_data_x <- lapply(test_data, function(x) x$x)
test_data_y <- lapply(test_data, function(x) x$y)
rm(test_data)
         
#predicting on other species
simpleCache(cacheName = "methTest_noRand_uc",instruction = {test_on_other(fit = res$model[[1]], test_data_x = test_data_x, test_data_y = test_data_y, ifRand = "noRand",test_species = annot_WGBS$Species)},
            cacheDir = paste0("RCache"),assignToVariable = "resTest",recreate = FALSE)
                      
my_wt(resTest$roc_dt, "roc_dt.tsv")
unique_auc <- unique(resTest$roc_dt[, c("type", "auc", "ifRand")])
unique_auc$AUC_train <- res$roc_dt$auc[[1]]
unique_auc$species_train <- species                     
my_wt(unique_auc, "auc_test.tsv")

                    
