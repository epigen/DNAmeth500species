## ---------------------------
##
## WGBS-based validation
##
## Downloading the avaliable WGBS data
##
## Authors: Daria Romanovskaia
##
## Date Created: 2021-05-20
##
##
## ---------------------------

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

##needed coverage threshold
COV=5

save_bed <- function(df, path){
    colnames(df) <- c("chr", "start", "end","strand", "coverage", "perc_meth_CpG")
    write.table(df, path, sep = "\t", row.names = F, quote = F)
}

save_bed_nostrand <- function(df, path){
    colnames(df) <- c("chr", "start", "end","coverage", "perc_meth_CpG")
    write.table(df, path, sep = "\t", row.names = F, quote = F)
}

WGBS_dir <- file.path(data_dir, "resources", "WGBS_public")

WGBS_outdir <- file.path(analysis_dir, "validation", "03_WGBS")

wd <- file.path(WGBS_outdir, "03.2_mean_meth")
dir.create(wd)
setwd(wd)

### allC for the Xenopus:
files_list <- list.files(file.path(WGBS_dir, "Xenopus_laevis/"), 
                         pattern = "*allC.txt", recursive = T, full.names = T)
files_list

annot <-  data.frame(path = files_list[[1]],tissue = "embryo",
           replica = 2, path_unif = gsub("allC.txt", ".bed",files_list[[1]] ))


files_list[[1]]

df <- fread(files_list[[1]]) #change to fread, way faster
head(df)

colnames(df) <- c("chr", "position", "strand", "threenucl_context", "mC_count", "coverage")

#Filtering out non-CpG methylation

df <- df[threenucl_context=="CG"]

df$meth_ratio = df$mC_count/df$coverage

df[,end:=position + 1,]

save_bed(df[,c("chr", "position", "end", "strand", "coverage", "meth_ratio")], as.character(annot$path_unif[1]))

stop()

WGBS_RRBS <- fread("../WGBS_RRBS_match.txt") #hand-created file with available WGBS data
head(WGBS_RRBS)

unique(WGBS_RRBS$`scientific name`)

WGBS_RRBS[`scientific name` == "Crassostrea gigas"]

## data type 1 - coverage & ratio (Gallus gallus all tissues except brain)

files_list <- list.files(file.path(WGBS_dir, "galGal5"), pattern = "*.tsv", full.names = TRUE, recursive = T)

annot <- data.frame(path = files_list,
            path_unif= as.character(sapply(files_list, function(x) gsub(".tsv", ".bed",x))),
           tissue = as.character(sapply(files_list, function(x) gsub("_tissue", "", strsplit(x, "/")[[1]][10]))),
                    replica = as.character(sapply(files_list, function(x) strsplit(strsplit(x, "/")[[1]][11], "_")[[1]][2])), stringsAsFactors = FALSE)

write.table(annot, file.path(WGBS_dir, "galGal5", "WGBS_annot.csv"), row.names = FALSE, quote = FALSE, sep = ";")

meth_file <- "galGal5_GSE146620.csv"

if(!file.exists(meth_file)){
mean_ratio = data.table()
for(i in c(1:NROW(annot))){
    print(annot$path[[i]])
    df <- fread(annot$path[[i]])
    print(NROW(df))
    ### saving the dataframe in the universal format (coverage threshold 0)
    save_bed_nostrand(df[cov>0, c("chr", "start", "end", "cov", "ratio")], annot$path_unif[[i]])
    ## filtering by the COV threshold
    df_sub <- df[df$cov > COV, ]
    print(NROW(df_sub))
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = mean(df_sub$ratio)))
}
mean_ratio <- left_join(mean_ratio, annot)
write.table(mean_ratio, meth_file, sep = ";", quote = F, row.names = F)
}else mean_ratio <- fread(meth_file)



head(mean_ratio)



ggplot(mean_ratio, aes(x = tissue, y = mean_ratio, color = replica)) + geom_point(size = 5) + theme_bw() + 
theme(text = element_text(size = 20))
ggsave(gsub("csv", "pdf", meth_file), width = 6, height = 4)

files_list <- list.files(file.path(WGBS_dir, "phaCin_unsw_v4.1"), pattern = "*CpG_report.txt", recursive = T, full.names = T)

annot <- data.frame(path = files_list, 
            path_unif= as.character(sapply(files_list, function(x) gsub("_CpG_report.txt", ".bed",x))),
           tissue = as.character(sapply(files_list, function(x) strsplit(x, "/")[[1]][10])),
           replica = as.character(sapply(files_list, function(x) ifelse(length(grep("PC", x))==1, 1,2))), stringsAsFactors = FALSE)
head(annot)

write.table(annot, file.path(WGBS_dir, "phaCin_unsw_v4.1", "WGBS_annot.csv"), 
            row.names = FALSE, quote = FALSE, sep = ";")

meth_file <- "phaCin_unsw_v4.1_GSE149600.csv"

if(!file.exists(meth_file)){
mean_ratio = data.table()
for(i in c(1:NROW(annot))){
    print(annot$path[[i]])
    df <- fread(annot$path[[i]])
    print(NROW(df))
    colnames(df) <- c("chr", "position", "strand", "meth_count", "unmeth_count", "cytosine_context", "trinucleotide_context")
    df <- df[df$cytosine_context == "CG", ]
    df_f <- df %>% mutate(total_count = meth_count + unmeth_count) %>%
        filter(total_count > 0) %>%
        mutate(meth_ratio = meth_count/total_count, end = position + 1)
    save_bed(df_f[,c("chr", "position", "end", "strand", "total_count", "meth_ratio")], annot$path_unif[[i]])
    df_f <- df_f %>% filter(total_count > COV)
    print(NROW(df_f))
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = mean(df_f$meth_ratio, na.rm = T)))
}
mean_ratio <- left_join(mean_ratio, annot)
write.table(mean_ratio, meth_file, sep = ";", quote = F, row.names = F)
}else mean_ratio <- fread(meth_file)

options(repr.plot.width = 10)
ggplot(mean_ratio, aes(x = tissue, y = mean_ratio, color = as.factor(replica))) + geom_point(size = 5) + 
                            theme_bw() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(gsub("csv", "pdf", meth_file), width = 6, height = 4)

files_list <- list.files(file.path(WGBS_dir, "bosTau9"), pattern = "*CpG_report.txt", recursive = T, full.names = T)

files_list[1]

annot <- data.frame(path = files_list, 
            path_unif= as.character(sapply(files_list, function(x) gsub(".CpG_report.txt", ".bed",x))),
           tissue = as.character(sapply(files_list, function(x) strsplit(x, "/")[[1]][10])),
           replica = as.character(sapply(files_list, function(x) strsplit(x, "_")[[1]][4] )), stringsAsFactors = FALSE)
head(annot)

write.table(annot, file.path(WGBS_dir, "bosTau9","WGBS_annot.csv"), row.names = FALSE, quote = FALSE, sep = ";")

meth_file <- "bosTau9_GSE147087.csv"

if(!file.exists(meth_file)){
mean_ratio = data.table()
for(i in c(1:NROW(annot))){
    print(annot$path[[i]])
    df <- fread(annot$path[[i]])
    print(NROW(df))
    colnames(df) <- c("chr", "position", "strand", "meth_count", "unmeth_count", "cytosine_context", "trinucleotide_context")
    df <- df[df$cytosine_context == "CG", ]
    df_f <- df %>% mutate(total_count = meth_count + unmeth_count) %>%
        filter(total_count > 0) %>% ## theoretically can be filtered by 5
        mutate(meth_ratio = meth_count/total_count, end = position + 1)
    save_bed(df_f[,c("chr", "position", "end", "strand", "total_count", "meth_ratio")], annot$path_unif[[i]])
    df_f <- df_f %>% filter(total_count > COV)
    print(NROW(df_f))
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = mean(df_f$meth_ratio, na.rm = T)))
}
mean_ratio <- left_join(mean_ratio, annot)
write.table(mean_ratio, meth_file, sep = ";", quote = F, row.names = F)
}else mean_ratio <- fread(meth_file)

gsub("csv", "pdf", meth_file)

options(repr.plot.width = 10)
ggplot(mean_ratio, aes(x = tissue, y = mean_ratio, color = replica)) + geom_point(size = 5) + 
                            theme_bw() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(gsub("csv", "pdf", meth_file), width = 6, height = 4)

files_list <- list.files(file.path(WGBS_dir, "mm9"), full.names = T,
                         pattern = "*CpG.calls.txt", recursive = T)

annot <- data.frame(path = files_list, 
             path_unif= as.character(sapply(files_list, function(x) gsub(".CpG.calls.txt", ".bed",x))),
           tissue = as.character(sapply(files_list, function(x) gsub("mouse_", "", strsplit(x, "/")[[1]][10]))),
                    replica = 1, stringsAsFactors = FALSE)

head(annot)

write.table(annot, file.path(WGBS_dir, "mm9", "WGBS_annot.csv"), row.names = FALSE, quote = FALSE, sep = ";")

meth_file <- "mm9_GSE42836.csv"

if(!file.exists(meth_file)){
mean_ratio = data.table()
for(i in c(1:NROW(annot))){
    print(annot$path[[i]])
    df <- fread(as.character(annot$path[[i]]))

    print(NROW(df))
    colnames(df)[4] <- "CpG_meth"
    colnames(df)[5] <- "CpG_seq"
    colnames(df)[6] <- "CpG_meth_prc"
    df <- df[CpG_seq > 0]
    save_bed_nostrand(df[,c("chromosome", "CpG location left", "CpG location right", "CpG_seq", "CpG_meth_prc")], annot$path_unif[[i]])
    
    df_f <- df[CpG_seq > 5]
    
    print(NROW(df_f))
    
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = mean(df_f$CpG_meth_prc, na.rm = T)/100))
}
mean_ratio <- left_join(mean_ratio, annot)
write.table(mean_ratio, meth_file, sep = ";", quote = F, row.names = F)
}else mean_ratio <- fread(meth_file)

options(repr.plot.width = 10)
ggplot(mean_ratio, aes(x = tissue, y = mean_ratio, color = replica)) + geom_point(size = 5) + 
                            theme_bw() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(gsub("csv", "pdf", meth_file), width = 7, height = 4)

files_list <- list.files(file.path(WGBS_dir, "Bl71nemr_1"), pattern = "*allC.txt", recursive = T, full.names = T)

df <- fread(files_list[[1]], header = FALSE) #change to fread, way faster
head(df)

colnames(df) <- c("chr", "position", "strand", "threenucl_context", "mC_count", "coverage")

#Filtering out non-CpG methylation

df <- df[threenucl_context=="CG"]

ggplot(df[coverage>2], aes(x = mC_count)) + geom_histogram(bins = 50)

df$meth_ratio = df$mC_count/df$coverage

df[,end:=position + 1,]

ggplot(df[coverage>2], aes(x = meth_ratio)) + geom_histogram(bins = 50)

save_bed(df[,c("chr", "position", "end", "strand", "coverage", "meth_ratio")], file.path(WGBS_dir, "Bl71nemr", "GSM2728830_liver.bed")) ## moving the file to the rest of the Bl17nmr 

annot <- data.frame(path = files_list, 
           tissue = "liver",
           replica = "GSE102144")

annot

write.table(annot,file.path(WGBS_dir, "Bl71nemr_1", "WGBS_annot.csv"))

df <- df[coverage > 5, ]

mean_ratio = data.table(path = files_list[[1]], mean_ratio = mean(df$meth_ratio), replica = 1, tissue = "liver", genomeid = "Bl71nemr_1")

write.table(mean_ratio, "Bl71nemr_GSE102144.csv", sep = ";", quote = F, row.names = F)

files_list <- list.files(file.path(WGBS_dir, "oyster.v9"), pattern = "*Cdepth.txt", recursive = T, full.names = T)

files_list[[1]]

annot <- data.frame(path = files_list, 
           tissue = "mantle",
           replica = gsub(".Cdepth.txt", "", as.character(sapply(files_list, function(x) strsplit(x, "_")[[1]][8] ))))
head(annot)

meth_file <- "oyster.v9_GSE40302.csv"

### extremly long - running as a script on the side
if(!file.exists(meth_file)){
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
    
    df$total_rep1 = df$mC_dep_rep1 + df$non_mC_dep_rep1
    df$total_rep2 = df$mC_dep_rep2 + df$non_mC_dep_rep2
 
    
    ## mean in repl 1:
    df_r1 <- df[total_rep1 > 5]
    df_r1$meth_ratio = df_r1$mC_dep_rep1/df_r1$total_rep1
   
    mean_r1 = mean(df_r1$meth_ratio)
    
    
    
    ## mean in repl 2:
    df_r2 <- df[total_rep2 > 5]
    df_r2$meth_ratio = df_r2$mC_dep_rep2/df_r2$total_rep2
   
    mean_r2 = mean(df_r2$meth_ratio)
    
    
    mean_ratio <- rbind(mean_ratio, data.frame(path = as.character(annot$path[[i]]), mean_ratio = c(mean_r1, mean_r2), repl = c(1, 2)))
}

mean_ratio <- left_join(mean_ratio, annot)

mean_ratio$replica <- paste0(mean_ratio$replica, "_", mean_ratio$repl)
write.table(mean_ratio[,-3], meth_file, sep = ";", quote = F, row.names = F)
}else mean_ratio <- fread(meth_file)


options(repr.plot.width = 5)
ggplot(mean_ratio, aes(x = tissue, y = mean_ratio, color = replica)) + geom_point(size = 5) + 
                            theme_bw() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(gsub("csv", "pdf", meth_file), width = 3, height = 4)

files_list <- list.files(WGBS_dir, pattern = "*.CGmap", recursive = T, full.names = T) #here we collect across genome folders

files_list

annot <- data.frame(path = files_list, 
           organism = as.character(sapply(files_list, function(x) strsplit(x, "/")[[1]][9])),
           tissue = as.character(sapply(files_list, function(x) strsplit(x, "/")[[1]][10])),
           replica = 1)
head(annot)

                                 
 ##saving lancet  
#annot_tmp <- rbind(annot, data.frame(path = files_list[[2]], 
 #          tissue = "Neural_tube",
  #         replica = "GSE141609"))

w#rite.table(annot_tmp, "../../resources/WGBS_public/Bl71nemr/WGBS_annot.csv", row.names = FALSE, quote = FALSE, sep = ";")

                              
i <- 2
df <- fread(as.character(files_list[[i]]))
colnames(df) <- c("chr", "nucleotide", "position", "context", "dinucl_context", "meth_level", "meth_read_count", "all_C_count")
df <- df[context=="CG"]
df[, end:=position + 1, ]

save_bed_nostrand(df[,c("chr", "position", "end", "all_C_count", "meth_level")], file.path(WGBS_dir, "Bl71nemr", "GSM4209500_mC_NeuralTube.bed")

meth_file <- "Combined_study_GSE141609.csv"

if(!file.exists(meth_file)){
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
write.table(mean_ratio, meth_file, sep = ";", quote = F, row.names = F)
}else mean_ratio <- fread(meth_file)

head(mean_ratio)

options(repr.plot.width = 10)
ggplot(mean_ratio, aes(x = organism, y = mean_ratio, color = tissue)) + geom_point(size = 5) + 
                            theme_bw() + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(gsub("csv", "pdf", meth_file), width = 10, height = 4)

annot <- fread(file.path(WGBS_dir, "danRer11_gemBS", "BSX_0021_Zebrafish_WGBS_wgbs_gembs_samples.csv"))
head(annot)

annot[, path:=paste0(file.path(WGBS_dir, "danRer11_gemBS","extract/"), dataset, "_cpg.bed.gz"),]

annot[, path_bed:=paste0(WGBS_dir, "/Danio_rerio_", library,"/", sample_name, ".bed"),]

colnames(annot)[16] <- "path_unif" 

head(annot)

write.table(annot[library == "GSE149416",c("path","path_unif","tissue","replica")], paste0(WGBS_dir, "/Danio_rerio_GSE149416/WGBS_annot.csv")), row.names = FALSE, quote = FALSE, sep = ";")

write.table(annot[library == "GSE134055",c("path","path_unif","tissue","replica")], paste0(WGBS_dir, "/Danio_rerio_GSE134055/WGBS_annot.csv")), row.names = FALSE, quote = FALSE, sep = ";")

meth.file <- "danRer11_GSE149416.csv"

if(!file.exists(meth.file)){
mean_ratio = data.table()
for(i in c(1:NROW(annot[library == "GSE149416"]) + 1)){  ##processed the first & blood from the second
    print(annot$path[[i]])
    df <- fread(annot$path[[i]], skip = 1)
    save_bed(df[, c(1,2,3,6,10,11)], annot$path_bed[[i]])
    mean_r = mean(df[V10 > 5]$V11)
    mean_ratio <- rbind(mean_ratio, data.frame(path = annot$path[[i]], mean_ratio = mean_r))
    mean_ratio <- left_join(mean_ratio, annot)
    write.table(mean_ratio[library == "GSE149416",],  "../validation/WGBS_public/Danio_rerio_GSE149416.csv", sep = ";", quote = F, row.names = F)
    }   
}else mean_ratio <- fread(meth.file)





mean_ratio <- mean_ratio[library == "GSE134055",c(path, mean_ratio)]

mean_ratio

annot_double <- annot[c((NROW(annot[library == "GSE149416"])+2):NROW(annot))]

df1 <- fread(annot_double$path[[1]], skip =1)

df2 <- fread(annot_double$path[[2]], skip = 1)

head(df1)

head(df2)

colnames(df1)[c(1:3)] <- c("chr", "start", "end")
colnames(df2)[c(1:3)] <- c("chr", "start", "end")

colnames(df1)[6] <- "strand"
colnames(df2)[6] <- "strand"

colnames(df1)[c(10,11)] <- c("coverage_1", "perc_1")
colnames(df2)[c(10,11)] <- c("coverage_2", "perc_2")

df <- full_join(df1[, c("chr", "start", "end","strand", "coverage_1", "perc_1" )], df2[, c("chr", "start", "end","strand", "coverage_2", "perc_2" )])

head(df)

weighted.mean(c(100,100), c(3,4))

df[is.na(coverage_1), coverage_1 :=0, ]

df[is.na(coverage_2), coverage_2 :=0, ]

df[, cov_total := coverage_1 + coverage_2, ]

df[,perc_total:= weighted.mean(c(perc_1, perc_2), c(coverage_1, coverage_2)), by = list(chr, start, end)]

save_bed(df[,c("chr", "start", "end", "strand", "cov_total", "perc_total")], annot_double$path_bed[1])

mean_r = mean(df$perc_total)

mean_r

mean_ratio

mean_ratio <- rbind(mean_ratio, data.frame(path = annot_double$path[[1]], mean_ratio = mean_r))

for(i in seq(3, NROW(annot_double), 2)){
    print(annot_double$path[[i]])
    df1 <- fread(annot$path[[i]], skip = 1)
    df2 <- fread(annot$path[[i+1]], skip = 1)
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
    save_bed(df[,c("chr", "start", "end", "strand", "cov_total", "perc_total")], annot_double$path_bed[1])
    
    
    mean_r = mean(df$perc_total)
    
    mean_ratio <- rbind(mean_ratio, data.frame(path = annot_double$path[[i]], mean_ratio = mean_r))
}
