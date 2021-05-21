library(ggplot2)
library(data.table)
library(dplyr)


files_list <- paste0("../../resources/WGBS_public/", list.files("../../resources/WGBS_public/", pattern = "*.CGmap", recursive = T))
annot <- data.frame(path = files_list, 
           organism = as.character(sapply(files_list, function(x) strsplit(x, "/")[[1]][5])),
           tissue = as.character(sapply(files_list, function(x) strsplit(x, "/")[[1]][6])),
           replica = 1)
head(annot)

mean_ratio = data.table()
for(i in c(1:NROW(annot))){
    print(annot$path[[i]])
    df <- fread(as.character(annot$path[[i]]))
    print(NROW(df))
    colnames(df) <- c("chr", "nucleotide", "position", "context", "dinucl_context", "meth_level", "meth_read_count", "all_C_count")
    df <- df[context=="CG"]
    df <- df[all_C_count > 5]
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = mean(df$meth_level)))
}

mean_ratio <- left_join(mean_ratio, annot)

write.table(mean_ratio, "../validation/WGBS_public/Combined_study_GSE141609.csv", sep = ";", quote = F, row.names = F)