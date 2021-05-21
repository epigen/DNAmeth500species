library(ggplot2)
library(data.table)
library(dplyr)
library(ggrepel)

COV=5

save_bed <- function(df, path){
    colnames(df) <- c("chr", "start", "end","strand", "coverage", "perc_meth_CpG")
    write.table(df, path, sep = "\t", row.names = F, quote = F)
}

mean_ratio = data.table()

annot_double <- fread("../../resources/WGBS_public/Danio_rerio_GSE134055/WGBS_annot.csv")

### saving the blood (one file only)
df <- fread(annot_double$path[[1]], skip = 1)
save_bed(df[, c(1,2,3,6,10,11)], annot_double$path_unif[[1]])
print("saved first")
## filter
mean_r = mean(df[V10 > 5]$V11)
mean_ratio <- rbind(mean_ratio, data.frame(path = annot_double$path[[1]], mean_ratio = mean_r))


### for loop over the paired samples
for(i in seq(2, NROW(annot_double), 2)){
    print(annot_double$path[[i]])
    df1 <- fread(annot_double$path[[i]], skip = 1)
    df2 <- fread(annot_double$path[[i+1]], skip = 1)
    ## transforming
    colnames(df1)[c(1:3)] <- c("chr", "start", "end")
    colnames(df2)[c(1:3)] <- c("chr", "start", "end")
    
    colnames(df1)[6] <- "strand"
    colnames(df2)[6] <- "strand"
    
    colnames(df1)[c(10,11)] <- c("coverage_1", "perc_1")
    colnames(df2)[c(10,11)] <- c("coverage_2", "perc_2")
    
    df <- full_join(df1[, c("chr", "start", "end","strand", "coverage_1", "perc_1" )], df2[, c("chr", "start", "end","strand", "coverage_2", "perc_2" )])
    
    df[is.na(coverage_1), coverage_1 :=0, ]
    df[is.na(coverage_2), coverage_2 :=0, ]
    
    df[, cov_total := coverage_1 + coverage_2, ]
    
    df[,perc_total:= weighted.mean(c(perc_1, perc_2), c(coverage_1, coverage_2)), by = list(chr, start, end)]
    save_bed(df[,c("chr", "start", "end", "strand", "cov_total", "perc_total")], annot_double$path_unif[i])
    
    
    mean_r = mean(df$perc_total)
    
    mean_ratio <- rbind(mean_ratio, data.frame(path = annot_double$path[[i]], mean_ratio = mean_r))
}

mean_ratio <- left_join(mean_ratio, annot_double)

write.table(mean_ratio,  "../validation/WGBS_public/Danio_rerio_GSE134055.csv", sep = ";", quote = F, row.names = F)