source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(tibble)
wd <- paste0(analysis_dir, "/02_predict_meth/02.8_why_inv/repeats/repeat_frequencies_6")
setwd(wd)

read_files <- function(path, k){
  if(k<10) files_high <- system(paste0("ls */kmer_fr_", path, "/kmer0", k),intern=TRUE)
  else files_high <- system(paste0("ls */kmer_fr_", path, "/kmer", k),intern=TRUE)
  sp_list<-as.character(sapply(files_high, function(x) strsplit(x, "/")[[1]][1] ))
  stats_high<-lapply(files_high, read.csv, sep = " ")
  stats_new <- lapply(seq_along(stats_high), function(i){ a <- stats_high[[i]][, c(-3)]; 
                                            colnames(a) <- c("kmer", sp_list[[i]]); return(a)})
  stats_high <- Reduce(full_join, stats_new )
  return(stats_high)
}

stats_high <- read_files("high", 6)
write.csv(stats_high, paste0(wd,"/kmer_freq_high_6.csv"), quote = F)

stats_low <- read_files("low", 6)
write.csv(stats_low, paste0(wd,"/kmer_freq_low_6.csv"), quote = F)

##same for k = 9
stats_high <- read_files("high", 9)
write.csv(stats_high, paste0(wd,"/kmer_freq_high_9.csv"), quote = F)

stats_low <- read_files("low", 9)
write.csv(stats_low, paste0(wd,"/kmer_freq_low_9.csv"), quote = F)

##same for k = 12
stats_high <- read_files("high", 12)
write.csv(stats_high, paste0(wd,"/kmer_freq_high_12.csv"), quote = F)

stats_low <- read_files("low", 12)
write.csv(stats_low, paste0(wd,"/kmer_freq_low_12.csv"), quote = F)

files_high <- system("ls */kmer_fr_high/kmer06",intern=TRUE)
sp_list<-as.character(sapply(files_high, function(x) strsplit(x, "/")[[1]][1] ))
#f <- read.csv(files_high[[1]])
stats_high<-lapply(files_high, read.csv, sep = " ")
stats_new <- lapply(seq_along(stats_high), function(i){ a <- stats_high[[i]][, c(-3)]; 
          colnames(a) <- c("kmer", sp_list[[i]]); return(a)})
stats_high <- Reduce(full_join, stats_new )


files_low <- system("ls */kmer_fr_low/kmer06",intern=TRUE)
sp_list<-as.character(sapply(files_low, function(x) strsplit(x, "/")[[1]][1] ))
#f <- read.csv(files_high[[1]])
stats_low<-lapply(files_low, read.csv, sep = " ")
stats_new <- lapply(seq_along(stats_low), function(i){ a <- stats_low[[i]][, c(-3)]; 
          colnames(a) <- c("kmer", sp_list[[i]]); return(a)})
stats_low <- Reduce(full_join, stats_new )

write.csv(stats_high, paste0(wd,"/kmer_freq_high.csv"), quote = F)
write.csv(stats_low, paste0(wd,"/kmer_freq_low.csv"), quote = F)


##same for the frequency of 9

NCOL(stats_high[[1]])
stats_high[[1]] <- stats_high[[1]][, c(-3)]
colnames(stats_high[[1]]) <- c("kmer", sp_list[[1]])
head(rownames_to_column(stats_high[[1]]))

Reduce(full_join, stats_high)

stats_files=system("ls */stats.csv",intern=TRUE)
stats_files
stats_fr<-lapply(stats_files, read.csv,row.names = 1)
sp_list<-unlist(lapply(stats_files, function(x) strsplit(x, "/")[[1]][1] ))

freq_df <- rbindlist(stats_fr)

freq_df$delta <- freq_df$high - freq_df$low

ggplot(freq_df, aes(x = species, y = delta)) + geom_bar(stat = "identity") + facet_wrap(~kmer)

med_auc <- read.csv(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/median_auc_Actinopteri.csv"),
                    row.names = 1)
freq_df[freq_df$species == "WHH", ]

freq_and_med <- left_join(freq_df, med_auc, by = c("species" = "train_species"))

ggplot(freq_and_med, aes(x = delta, y = median_auc, color = kmer)) + geom_point()

pdf("summary_differences.pdf", width = 20, height = 5)
ggplot(freq_and_med, aes(x = species, y = delta, fill = kmer)) + geom_bar(stat = "identity", position = "dodge") + rotate_labels()
dev.off()