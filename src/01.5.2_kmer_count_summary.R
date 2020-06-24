#!/bin/env Rscript
#init
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

##function to read in all the files and generate the full table
##input: list of files, id to save and list of species, corresponding to the file paths
#saves to wd+summary/kmer0{i}.csv

make_table<-function(x,i, sp_list){
  print(i)
  file_list<-mapply(function(x, y) read.table(x, skip = 1, row.names = 1, col.names = c("kmer",y)), x, sp_list, SIMPLIFY = FALSE)
  t<-Reduce(bind_cols, file_list)
  row.names(t)<-row.names(file_list[[1]])
  write.csv(t, paste0("summary/kmer0", i, ".csv"), quote = FALSE)
}


wd=file.path(analysis_dir, "99.4_kmer_count_filtered")
setwd(wd)

##collecting files 
stats_files<-lapply(c(1,2,3,4,5,6,7,8,9), function(x) system(paste0("ls */kmer0",x),intern=TRUE))
sp_list<-unlist(lapply(stats_files[[1]], function(x) strsplit(x, "/")[[1]][1] ))
lapply(seq_along(stats_files), function(i) make_table(stats_files[[i]], i, sp_list))

##update: order 10 in folders SP_10
stats_files_10<-system("ls *_10/kmer10",intern=TRUE)
sp_list<-as.character(sapply(stats_files_10, function(x) strsplit(x, "_")[[1]][1]))
make_table(stats_files_10, 10, sp_list )

##and now the same on the 6 to 9
stats_files<-lapply(c(1,2,3,4,5,6,7,8,9), function(x) system(paste0("ls *_10/kmer0",x),intern=TRUE))
sp_list<-as.character(sapply(stats_files[[1]], function(x) strsplit(x, "_")[[1]][1]))
lapply(seq_along(stats_files), function(i) make_table(stats_files[[i]], i, sp_list))

