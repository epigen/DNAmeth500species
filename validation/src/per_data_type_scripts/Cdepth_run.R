library(ggplot2)
library(data.table)
library(dplyr)


files_list <- paste0("../../resources/WGBS_public/Crassostrea_gigas/", list.files("../../resources/WGBS_public/Crassostrea_gigas/", pattern = "*Cdepth.txt", recursive = T))

annot <- data.frame(path = files_list, 
           tissue = "mantle",
           replica = as.character(sapply(files_list, function(x) strsplit(x, "_")[[1]][7] )))
head(annot)
write.table(annot, "../../resources/WGBS_public/Crassostrea_gigas/WGBS_annot.csv", row.names = FALSE, quote = FALSE, sep = ";")

mean_ratio = data.table()
for(i in c(1:NROW(annot))){
    print(annot$path[[i]])
    df <- fread(as.character(annot$path[[i]]), skip = 26)
    print(NROW(df))
    df <- df[type == "CG"]
    
    colnames(df)[6] <- "mC_dep_rep1"
    colnames(df)[7] <- "non_mC_dep_rep1"
    
    colnames(df)[9] <- "mC_dep_rep2"
    colnames(df)[10] <- "non_mC_dep_rep2"

    df[, total_rep1 := mC_dep_rep1 + non_mC_dep_rep1,]
    df[, total_rep2 := mC_dep_rep2 + non_mC_dep_rep2,]
    print(head(df))
    
    ## mean in repl 1:
    df_r1 <- df[total_rep1 > 5]
    print(NROW(df_r1))
    df_r1$meth_ratio = df_r1$mC_dep_rep1/df_r1$total_rep1
    print(head(df_r1))
    print(df_r1$meth_ratio)
    mean_r1 = mean(df_r1$meth_ratio)
    
    
    
    ## mean in repl 2:
    df_r2 <- df[total_rep2 > 5]
    df_r2$meth_ratio = df_r2$mC_dep_rep2/df_r2$total_rep2
   
    mean_r2 = mean(df_r2$meth_ratio)
    
    
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = c(mean_r1, mean_r2), repl = c(1, 2)))
}

mean_ratio <- left_join(mean_ratio, annot)
mean_ratio$replica <- paste0(mean_ratio$replica, "_", mean_ratio$repl)
                                         
write.table(mean_ratio[,-3], "../validation/WGBS_public/Crassostrea_gigas_GSE40302.csv", sep = ";", quote = F, row.names = F)