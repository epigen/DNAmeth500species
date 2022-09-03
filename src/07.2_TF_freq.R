#!/usr/bin/env Rscript
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

args = commandArgs(trailingOnly=TRUE)
species=args[1]

#load ???he meta table (we need to know the full number of methylated and unmethylated fragments)
wd=file.path(processed_dir, species)
print(wd)
setwd(wd)

meth_data_mean=fread(paste0("toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/",species,"_mean_meth.tsv"))
meth_data_mean_long=melt(meth_data_mean,id.vars=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),measure=patterns(".cov", ".meth"),variable.factor=TRUE,variable.name="sample",value.name=c("cov","meth"),na.rm=TRUE)
meth_data_mean_cond=meth_data_mean_long[,list(Nsamples=.N,mean_cov=mean(cov),min_cov=min(cov),max_cov=max(cov),mean_meth=mean(meth),min_meth=min(meth),max_meth=max(meth)),by=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
meth_data_mean_cond_red=meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000]
meth_data_mean_cond_red[,category:=ifelse(min_meth>80,"meth",ifelse(max_meth<20,"unmeth","amb")),]

#defining the number of methylated/unmethylated fragments
n <- list(meth=NROW(filter(meth_data_mean_cond_red, category=="meth")), unmeth=NROW(filter(meth_data_mean_cond_red, category=="unmeth")))

N <- NROW(meth_data_mean_cond_red)
rm(meth_data_mean)


wd=file.path(analysis_dir, "07_motifAnalysis/07.1_TFBS_detection_2020", species)
setwd(wd)

subdir=file.path(analysis_dir, "07_motifAnalysis/07.2_TF_freq_2020", species)
dir.create(subdir, recursive = TRUE)

#load the results table
TF_files = system("ls *_profile_table.tsv", intern = T)
TFs <- sapply(TF_files, function(x) strsplit(x, "_")[[1]][1])

get_counts <- function(res_path, TF){
  res <- read.csv(res_path, sep = "\t")
  print(paste0(TF, " ", NROW(res)))
  
  if( NROW(res) == 0 ){
    total_count <- data.frame(matrix(nrow = 1, ncol = 3, 0)) ##change to the dataframe
    colnames(total_count) <- c("TF", "n", "n_norm") 
    total_count$TF <- TF
    res_counts <- data.frame(matrix(nrow = 1, ncol = 6, 0)) ##change to the dataframe
    colnames(res_counts) <- c("TF", "meth", "unmeth", "meth_norm", "unmeth_norm", "f" ) 
    res_counts$TF <- TF
  }else{
    total_count<-res %>%
      count(TF) %>%
      mutate(n_norm=n/N)
    if(NROW(res[res$category!="amb",]) == 0){
      res_counts <- data.frame(matrix(nrow = 1, ncol = 3, 0)) ##change to the dataframe
      colnames(res_counts) <- c("TF", "meth", "unmeth") 
      res_counts$TF <- TF 
    }else{
    res_counts<-res %>%
        filter(category!="amb") %>%
        count(TF, category) %>%
        spread(category, n)
    if(!"meth" %in% colnames(res_counts)) res_counts$meth <-0 
    if(!"unmeth" %in% colnames(res_counts)) res_counts$unmeth <-0 
    }
    res_counts <- res_counts  %>%
      mutate(meth_norm=meth/n$meth, unmeth_norm=unmeth/n$unmeth, 
             f=meth_norm/unmeth_norm)
  }
  return(full_join(total_count, res_counts))
}

TF_fr_list <- lapply(seq_along(TF_files), function(i) get_counts(TF_files[[i]], TFs[[i]]))
TF_fr <- rbindlist(TF_fr_list)

TF_fr$N_ded_ref <- N
TF_fr$N_ded_ref_meth <- n$meth
TF_fr$N_ded_ref_unmeth <- n$unmeth

pdf(paste0(subdir, "/motif_frequency.pdf"))
ggplot(TF_fr, aes(x=meth_norm, y=unmeth_norm))+geom_point()+geom_abline(linetype="dotted")
dev.off()

write.csv(TF_fr, paste0(subdir, "/counts.csv"))
