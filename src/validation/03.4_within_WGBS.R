#from WGBS to WGBS
#use another environment
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

library(kebabs)
library(caret)
library(ROCR)

args = commandArgs(trailingOnly=TRUE)

species <- args[[1]]
genome_id <- args[[2]]
cov_min <- as.integer(args[[3]])
matched_species <- args[[4]]
#assembly <- args[[2]]
#abbreviation_sp <- args[[3]]


#kmer <- args[[4]]
kmer=c(1:10)
#kmer=c(1:3)
print(species)
print(cov_min)

get_train_test_split <- function(meth_data_mean_cond_red, ded_ref, N = 2000, SEED = "1234"){

  ## chromosome for train/split
  meth_data_mean_cond_red$chr <- as.character(sapply(meth_data_mean_cond_red$name, 
                                                function(x) {a <- strsplit(x, "_")[[1]]; 
                                                y <- ifelse(length(a) == 3, a[1], paste0(a[1], "_", a[2]));
                                                return(y)}))
  print(table(meth_data_mean_cond_red$category))
  
    ##adding the chromosome test/train split
  chr_list <- unique(meth_data_mean_cond_red$chr)
  chr_to_train <- sample(chr_list, length(chr_list)/2)
  print(chr_to_train)

  chr_to_test <- setdiff(chr_list, chr_to_train)
  print(chr_to_test)

  #minimum features of either group --> maximum number of features for each group
  Nsel=min(c(min(table(meth_data_mean_cond_red[chr %in% chr_to_train, ]$category)), 
             min(table(meth_data_mean_cond_red[chr %in% chr_to_test, ]$category))))
  
    if (Nsel > N/2) {
    Nsel = N/2
  }
    
  print(paste0("sample size ", Nsel))
  #set up features (sequences) and labels
  set.seed(SEED)

  x_train=c(ded_ref[sample(meth_data_mean_cond_red[category==1 & chr %in% chr_to_train]$name, size=Nsel)],
      ded_ref[sample(meth_data_mean_cond_red[category==-1 & chr %in% chr_to_train]$name,size=Nsel)])
  print(head(x_train))

  x_test=c(ded_ref[sample(meth_data_mean_cond_red[category==1 & chr %in% chr_to_test]$name,size=Nsel)],
      ded_ref[sample(meth_data_mean_cond_red[category==-1 & chr %in% chr_to_test]$name,size=Nsel)])
  print(head(x_test)) 

  y_train=c(sample(meth_data_mean_cond_red[category==1]$category,size=Nsel),
      sample(meth_data_mean_cond_red[category==-1]$category,size=Nsel))
  print(head(y_train))

  y_test=c(sample(meth_data_mean_cond_red[category==1]$category,size=Nsel),
      sample(meth_data_mean_cond_red[category==-1]$category,size=Nsel))
  print(head(y_test))

  return(list(x_train = x_train, x_test = x_test, y_train = y_train, y_test = y_test))
}

train_test <- function(x_train, y_train, x_test, y_test, ifRand, k= kmer, runid = 0){
 print("training the model")
 chr_group <- as.factor(as.character(sapply(names(x_train), 
                                                function(x) {a <- strsplit(x, "_")[[1]]; 
                                                y <- ifelse(length(a) == 3, a[1], paste0(a[1], "_", a[2]));
                                                return(y)})))
 ## defining the number of cross-validations
 N_chr = length(unique(chr_group))
 print(N_chr)
 if(N_chr <= 10){ n_cross = as.integer(0.7*N_chr)
                }else n_cross = 10

 specK = spectrumKernel(k=k)
  ms = kbsvm(x=x_train, y=y_train, kernel=specK,pkg="e1071", svm="C-svc", C=c(0.01,0.1,1,10), 
 # ms = kbsvm(x=x_train, y=y_train, kernel=specK, pkg="e1071", svm="C-svc", C=c(0.1, 1),
             explicit="yes",showProgress=TRUE, cross=n_cross, noCross=1, nestedCross=n_cross, groupBy=chr_group)
  print("model trained...")
  ms_res = modelSelResult(ms)
   ##optimal values
  best_ks = unlist(lapply(ms_res@selGridRow,function(x){return(x@k)}))
  best_Cs = ms_res@selGridCol$C
  sel_k=as.numeric(names(sort(table(best_ks),decreasing=TRUE)[1]))
  sel_C=as.numeric(names(sort(table(best_Cs[best_ks==sel_k]),decreasing=TRUE)[1]))

  ##optimal idea
  specK_sel <- spectrumKernel(k=sel_k)
  fit <- kbsvm(x=x_train, y=y_train, kernel=specK_sel, pkg="e1071", svm="C-svc", C=sel_C, explicit="yes",
               featureWeights="yes")
  ## feature weights
  options( warn = -1 )
  feature_weights=melt(sort(as.data.table(featureWeights(fit))))
  options( warn = 0 )
  
  ## fitting the test
  #preds <- predict( fit, x_test)
  #pref<-evaluatePrediction(preds, y_test,allLabels=unique(y_train), decValues=preds, print=FALSE)
  #print(pref$AUC)

  ##rocauc curve
  preddec <- predict(fit, x_test, predictionType="decision")
  rocdata <- computeROCandAUC(preddec, y_test, allLabels=unique(y_train))

  roc_dt=data.table(fdr=unlist(rocdata@FPR),
    tpr=unlist(rocdata@TPR),
    auc=unlist(rocdata@AUC), 
    ifRand=ifRand, k=sel_k, C=sel_C, run = runid)
  
  return(list(roc_dt=roc_dt, roc=rocdata, feature_weights=feature_weights,
                param=list(k=sel_k, C=sel_C), model=list(fit, ms_res)))
}
 

##setup
path_to_folder <- file.path(analysis_dir, "validation", "03_WGBS", "03.4_prediction", species)
path_to_data <- file.path(analysis_dir, "validation", "03_WGBS", "03.3_fragments", species)
### now, if we alreaady have run the pipeline and want to run it with another seed - we need to redirect it to other folder
##outdir
outdir <- file.path(path_to_folder, "kebabs_model")

if (!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
  SEED = "1234" ## first time, reproducible seed
}else{
  SEED <- as.character(sample(1:1000, 1))
}

##loading data
meth_data_mean = fread(file.path(path_to_data, "mean_meth_per_fragment.tsv"))

##formatting
meth_data_mean_long <- melt(meth_data_mean, measure=patterns("cov", "mean_meth"),
                            variable.factor=TRUE,
                         variable.name="sample",
                            value.name=c("cov","meth"),na.rm=TRUE)
                            
meth_data_mean_cond <- meth_data_mean_long[,list(Nsamples=.N, mean_cov=mean(cov), min_cov=min(cov),
                                              max_cov=max(cov), mean_meth=mean(meth), min_meth=min(meth),
                                              max_meth=max(meth)),
                                        by=c("name")]

print(NROW(meth_data_mean_cond))

meth_data_mean_cond_red=meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>cov_min & mean_cov<1000 &(min_meth>80|max_meth<20)]
print(NROW(meth_data_mean_cond_red))

meth_data_mean_cond_red[,category:=ifelse(min_meth>80,1,ifelse(max_meth<20,-1,NA)),]
print(NROW(meth_data_mean_cond_red[category==1]))
print(NROW(meth_data_mean_cond_red[category==-1]))

##loading sequence
ded_ref=readDNAStringSet(file.path(path_to_data, paste0(genome_id, "_fragments.fa")))


split_ds <- get_train_test_split(meth_data_mean_cond_red, ded_ref, SEED = SEED)


#original labels
simpleCache(cacheName=paste0("methPred_noRand_uc_", SEED), instruction={ train_test(x_train=split_ds$x_train,
            x_test = split_ds$x_test, y_train = split_ds$y_train, y_test = split_ds$y_test,
                                                        ifRand='noRand', k=kmer, runid = 0)},
            cacheDir=paste0(outdir,"/RCache"), assignToVariable="res", recreate=FALSE)

roc_res=res$roc_dt


## randomise labels
#random labels (several iterations)
rand_labs = lapply(seq_len(5), 
  function(x) { return(list(name=x, y_rand_train = sample(split_ds$y_train), 
                                              y_rand_test = sample(split_ds$y_test)))})

rand_res_parallel = mclapply(rand_labs, function(labels_rand) {simpleCache(cacheName=paste0("methPred_Rand_uc_",
        labels_rand$name, "_", SEED),
        instruction = {train_test(x_train = split_ds$x_train, 
                                  x_test = split_ds$x_test,
                                   y_train=labels_rand$y_rand_train, 
                                  y_test=labels_rand$y_rand_test, 
                                  ifRand='rand', runid=labels_rand$name,
                                  k=kmer)}, assignToVariable="rand_res",
                                        cacheDir=paste0(outdir,"/RCache"),
                                        buildEnvir=c(labels_rand),recreate=FALSE); return(rand_res)},
        mc.cores=5, mc.preschedule=FALSE)

rand_roc_res_parallel_list=lapply(rand_res_parallel,function(x){return(x$roc_dt)})
rand_roc_res_parallel=rbindlist(rand_roc_res_parallel_list,idcol="run")
rand_roc_res_parallel[,run:=run,]
rand_roc_res_parallel[,run:=NULL,]

#merge rand and original + plot
roc_res=rbindlist(list(roc_res,rand_roc_res_parallel))


auc_res=roc_res[,list(auc=mean(auc)),by=ifRand]
auc_res[,x:=0.9,]
auc_res[,y:=ifelse(ifRand=="rand", 0.09, 0.13),]


auc_table <- fread(file.path(analysis_dir, "05_predict_meth","05.1_within_species", "summary", "all_aucs.csv"))

if(matched_species == "ALL FROGS"){
    frogs <- unique(stats_annot[grep("frog", stats_annot$English), c("species", "English", "ncbi_name")])
    
    auc_res <- rbind(auc_res, list(ifRand="RRBS", auc = mean(auc_table[species %in% frogs$species, ]$AUC), x = 0.9, y = 0.17))
}else{
auc_res <- rbind(auc_res, list(ifRand="RRBS", auc = auc_table[species == matched_species, ]$AUC, x = 0.9, y = 0.17))}

pdf(paste0(outdir,"/", genome_id, "ROC.pdf"), height=4, width=5)
ggplot(roc_res, aes(x=fdr,y=tpr,col=ifRand)) + geom_line(aes(group=run, alpha=ifRand)) + 
    geom_text(data=auc_res,aes(x=x,y=y,label=paste0("auc=",round(auc,3)))) +
    scale_color_manual(values=c("rand"="grey","noRand"="blue", "RRBS" = "red"))+
    scale_alpha_manual(values=c("rand"=0.5,"noRand"=1)) + 
    theme(text = element_text(size = 15))
dev.off()

write.csv(roc_res, paste0(outdir,"/", species, "roc_res.csv"))

