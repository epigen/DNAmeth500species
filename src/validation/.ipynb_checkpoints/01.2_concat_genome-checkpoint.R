source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

library(rtracklayer)
library(Biostrings)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

db <- args[1]

ref_dir <- paste0(data_dir, "/resources/reference_genomes/", db, "/")
dir.create(ref_dir, recursive = T)
setwd(ref_dir)


###downloading the data
### genome (soft-mask or unmasked)
if(file.exists(paste0(ref_dir, db, ".fa"))){
print("already downloaded!")
}else{

link_to_download <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", db, "/bigZips/", db, ".fa.gz")
cmd = paste0("wget ", link_to_download, " -P ", ref_dir)
print(cmd)
system(cmd)

## chrom sizes:
link_to_download <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", db, "/bigZips/", db, ".chrom.sizes")
cmd = paste0("wget ", link_to_download, " -P ", ref_dir)
print(cmd)
system(cmd)


## unzipping the flie:
cmd = paste0("gunzip ",ref_dir, db, ".fa.gz")
print(cmd)
system(cmd)
}

###uploading the files to work with
chrom_sizes = read.csv(paste0(db, ".chrom.sizes"), sep = "\t", header = F)
head(chrom_sizes)

if(NROW(chrom_sizes) < 30){
    print("Genome assembly good enough")
    cmd = paste0("cp ", ref_dir, db, ".fa ", ref_dir, db, "_concatinated.fa")
    system(cmd)
    chrom.map = data.table(new = chrom_sizes$V1, old = chrom_sizes$V1, width = chrom_sizes$V2)
    write.table(chrom.map, paste0( db, "_chrom_mapping.tsv"), sep = "\t", row.names = F, quote = F)
    stop()
}
##expected size of the artificial chromosome
tile_size = round(sum(chrom_sizes$V2)/20)
print(tile_size)

setDT(chrom_sizes)

##sorting by the chromosome size
chrom_sizes <- chrom_sizes[order(V2, decreasing = T)]


##uploading fasta
fa_file <- readDNAStringSet(paste0(db, ".fa"))

##first saving the first one (no append=True on purpose!) - chromosome_id:
if (file.exists(paste0(db, "_concatinated.fa"))){
    print("removing the previous version")
    file.remove(paste0(db, "_concatinated.fa"))
}

##data frame that knows which one maps to what
chrom.map = data.table()

##id count
chrom_id = 0

##saving chromosomes with no size change
if(chrom_sizes$V2[[1]] > 0.8*tile_size){
for(chrom_id in c(1:NROW(chrom_sizes[V2>0.8*tile_size]))){
    print(chrom_id) 
    
    chrom.map <- rbind(chrom.map, list("new" = paste0("chrArt",chrom_id), 
                                       "old" = chrom_sizes$V1[[chrom_id]], 
                                       "width" = chrom_sizes$V2[[chrom_id]]))
    
    seq_obj <- DNAStringSet(as.character(fa_file[[chrom_sizes$V1[[chrom_id]]]]))
    names(seq_obj) <- paste0("chrArt",chrom_id)
    writeXStringSet(seq_obj,paste0(db, "_more_chr_concatinated.fa") , append=TRUE)
}}

print(chrom.map)

print("...now saving joined chromosomes")
## joined chromosomes should be separated by a linker
linker = paste(rep("N", 100), collapse = "")

##saving joined chromosomes
chrom_id = chrom_id + 1
j = chrom_id ##first row of the leftover

print(j)
while(j < NROW(chrom_sizes)){
    i = j ## getting the moving count
    used_length = chrom_sizes[i]$V2
    print(as.character(chrom_sizes[i]$V1))
    
    while(used_length < tile_size && i  < NROW(chrom_sizes) && i <1000){
        i = i+1
        used_length = used_length + chrom_sizes[i]$V2 ## next chromosome to add
    }
    
    chrom.map <- rbind(chrom.map, list("new" = paste0("chrArt",chrom_id), "old" = paste(as.character(chrom_sizes[c(j:i)]$V1), collapse = ";"), width = used_length ))
    
    seq = paste(as.character(fa_file[chrom_sizes[c(j:i)]$V1]), collapse = linker)
    seq_obj <- DNAStringSet(seq, use.names = T)
    names(seq_obj) <- paste0("chrArt",chrom_id)
    writeXStringSet(seq_obj,paste0(db, "_more_chr_concatinated.fa") , append=TRUE)
    
    j = i + 1 ## moving on to the next chromosome(that is not included yet)
    print(j)
    
    chrom_id <- chrom_id + 1  ## starting the new artificial chromosome
}

##save the mapping: 
write.table(chrom.map, paste0( db, "_chrom_mapping.tsv"), sep = "\t", row.names = F, quote = F)
    
print("Done")

