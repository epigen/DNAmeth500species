source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
sessionInfo()

library(kebabs)
library(caret)
library(ROCR)
library(seqLogo)


##custim function for calculating f1
get_f1<-function(PRECISION, RECALL){
  return (2*PRECISION*RECALL/((PRECISION+RECALL)*100))
}


## create x and y dataset separation function from the meth_data_mean_cond_red preprocessed
create_datasets<-function(meth_data_mean_cond_red, ded_ref, name_tag, TO_SAVE=TRUE, N=2000, SEED = TRUE){

  #minimum features of either group --> maximum number of features for each group
  Nsel=min(table(meth_data_mean_cond_red$category))
  if (Nsel > 2000) {
    Nsel = N
  }
  print(Nsel)
  #set up features (sequences) and labels
  if(SEED) set.seed("1234")
  
  x=c(ded_ref[sample(meth_data_mean_cond_red[category==1]$meta,size=Nsel)],
      ded_ref[sample(meth_data_mean_cond_red[category==-1]$meta,size=Nsel)])
  y=c(sample(meth_data_mean_cond_red[category==1]$category,size=Nsel),
      sample(meth_data_mean_cond_red[category==-1]$category,size=Nsel))
  stopifnot(length(x) == length(y))
  
  comb=data.table(id=names(x),seq=as.data.frame(x),lab=y)
  
  if(TO_SAVE){
    write.table(comb,paste0(subdir,"/sequences/sequence_table_", name_tag, ".tsv"),
                quote=FALSE, sep="\t", row.names=FALSE )
    
    frag_meth_count=table(meth_data_mean_cond_red$category)
    
    stats=data.table(Nmeth=frag_meth_count["1"],Nunmeth=frag_meth_count["-1"])
    
    write.table(stats,paste0(subdir,"/stats/", name_tag, "_meth_stats.tsv"),
                quote=FALSE,sep="\t",row.names=FALSE )
    
  }
  return(list(x=x, y=y))
  
}

##perform test and train for the initial analysis 
train_test<-function(x,y,type,ifRand, run,k, subdir, SAVE_TRAIN_IDS=FALSE){
  set.seed("1234")
  
  fold=createDataPartition(y=y,times=1)[[1]]
  x_train=x[-fold,,drop=FALSE]
  y_train=y[-fold]
  x_test=x[fold,,drop=FALSE]
  y_test=y[fold]
  
  print(length(y_test))
  print(length(y_train))
  
  if(SAVE_TRAIN_IDS){
    train_ids = names(x_train)
  }
  # Model selection
  specK = spectrumKernel(k=k)
  ms = kbsvm(x=x_train, y=y_train, kernel=specK,pkg="e1071", svm="C-svc", C=c(0.01,0.1,1,10), 
             explicit="yes",showProgress=TRUE, cross=10, noCross=1, nestedCross=10)
    ## show best parameter settings
  ms_res = modelSelResult(ms)
  best_ks = unlist(lapply(ms_res@selGridRow,function(x){return(x@k)}))
  best_Cs = ms_res@selGridCol$C
  sel_k=as.numeric(names(sort(table(best_ks),decreasing=TRUE)[1]))
  sel_C=as.numeric(names(sort(table(best_Cs[best_ks==sel_k]),decreasing=TRUE)[1]))
  specK_sel <- spectrumKernel(k=sel_k)
  fit <- kbsvm(x=x_train, y=y_train, kernel=specK_sel, pkg="e1071", svm="C-svc", C=sel_C, explicit="yes",
               featureWeights="yes")
  
  options( warn = -1 )
  feature_weights=melt(sort(as.data.table(featureWeights(fit))))
  options( warn = 0 )
  
  #predict
  preds <- predict( fit, x_test)
  pref<-evaluatePrediction(preds, y_test,allLabels=unique(y), decValues=preds, print=FALSE)
  print(pref)
  f1<-get_f1(pref$SENS, pref$PREC)
  preddec <- predict(fit, x_test, predictionType="decision")
  rocdata <- computeROCandAUC(preddec, y_test, allLabels=unique(y))
  #plot(rocdata)
  if(ifRand=="noRandTrain"){
    #save the key data: PWM and stats
    print("saving stuff")
    if (sel_k>1){
      pwm_low=makePWM(consensusMatrix(DNAStringSet(as.character(feature_weights[order(value,decreasing=FALSE)]$variable[1:10])),as.prob=TRUE)[1:4,])
      pwm_high=makePWM(consensusMatrix(DNAStringSet(as.character(feature_weights[order(value,decreasing=TRUE)]$variable[1:10])),as.prob=TRUE)[1:4,])
      print(subdir)
      
      pdf(paste0(subdir,"/seqLogo_low_", type, ".pdf"),height=3,width=4.5)
      seqLogo(pwm_low)
      dev.off()
      
      pdf(paste0(subdir,"/seqLogo_high_", type, ".pdf"),height=3,width=4.5)
      seqLogo(pwm_high)
      dev.off()
    }
    #create stats overview record to combine with other species
    k_freq=table(unlist(lapply(ms_res@selGridRow,function(x){x@k})))

    #!! No methylation statisctics here, saving this in the data readout
    stats=data.table(Species=type,k=fit@svmInfo@reqKernel@k,k_freq=max(k_freq/sum(k_freq)),c=fit@svmInfo@selSVMPar$cost,numSequences=fit@numSequences,AUC=rocdata@AUC, f1=f1)
    print(paste0(subdir,"/", type, "_stats.tsv"))
    write.table(stats,paste0(subdir,"/", type, "_stats.tsv"), quote=FALSE,sep="\t",row.names=FALSE )
  }
  
  roc_dt=data.table(fdr=unlist(rocdata@FPR),tpr=unlist(rocdata@TPR),auc=unlist(rocdata@AUC), f1=f1,run=run,type=type,ifRand=ifRand, k=sel_k, C=sel_C,min_motif=feature_weights[which.min(value)]$variable,max_motif=feature_weights[which.max(value)]$variable)
  
  if(SAVE_TRAIN_IDS){
    
    return(list(train_ids=train_ids, roc_dt=roc_dt,roc=rocdata,feature_weights=feature_weights,
                param=list(k=sel_k,C=sel_C),model=list(fit,ms_res)))
  }else{
    return(list(roc_dt=roc_dt,roc=rocdata,feature_weights=feature_weights,
                param=list(k=sel_k,C=sel_C),model=list(fit,ms_res)))
  }
}

