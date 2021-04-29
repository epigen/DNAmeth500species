library(rtracklayer)
library(Biostrings)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

db <- args[1]
genome_id <- args[2]

reference_path <- paste0("../../resources/reference_genomes/genomes/", genome_id, "/concat/")
dir.create(reference_path)



###downloading the data
### genome (soft-mask or unmasked)
if(file.exists(paste0(reference_path, db, ".fa"))){
	print("already downloaded!")
}else{

link_to_download <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", db, "/bigZips/", db, ".fa.gz")
cmd = paste0("wget ", link_to_download, " -P ", reference_path)
print(cmd)
system(cmd)

## chrom sizes:
link_to_download <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", db, "/bigZips/", db, ".chrom.sizes")
cmd = paste0("wget ", link_to_download, " -P ", reference_path)
print(cmd)
system(cmd)


## unzipping the flie:
cmd = paste0("gunzip ", reference_path, db, ".fa.gz")
print(cmd)
system(cmd)
}

###uploading the files to work with
chrom_sizes = read.csv(paste0(reference_path, db, ".chrom.sizes"), sep = "\t", header = F)
head(chrom_sizes)
##expected size of the artificial chromosome
tile_size = round(sum(chrom_sizes$V2)/20)
print(tile_size)

setDT(chrom_sizes)
##sorting by the chromosome size
chrom_sizes <- chrom_sizes[order(V2, decreasing = T)]


##uploading fasta
fa_file <- readDNAStringSet(paste0(reference_path, db, ".fa"))

##first saving the first one (no append=True on purpose!) - chromosome_id:
cat(">chrArt1", file = paste0(reference_path, db, "_concatinated.fa"), sep = "\n")
##and the seuqence
cat(as.character(fa_file[[chrom_sizes$V1[[1]]]]), file = paste0(reference_path, db, "_concatinated.fa"), sep = "\n", append = T)

##data frame that knows which one maps to what
chrom.map = data.table(new = "chrArt1", old = as.character(chrom_sizes$V1[[1]]), width = chrom_sizes$V2[[1]])


##saving chromosomes with no size change
for(chrom_id in c(2:NROW(chrom_sizes[V2>0.8*tile_size]))){
    print(chrom_id)
   cat(paste0(">chrArt",chrom_id), file = paste0(reference_path, db, "_concatinated.fa"), sep = "\n", append = T) 
    cat(as.character(fa_file[[chrom_sizes$V1[[chrom_id]]]]), file = paste0(reference_path, db, "_concatinated.fa"), sep = "\n", append = T)
     chrom.map <- rbind(chrom.map, list("new" = paste0("chrArt",chrom_id), "old" = as.character(chrom_sizes$V1[[chrom_id]]), "width" = chrom_sizes$V2[[chrom_id]]))
}

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
    while(used_length < tile_size && i  < NROW(chrom_sizes)){
        i = i+1
        used_length = used_length + chrom_sizes[i]$V2 ## next chromosome to add
    }
    chrom.map <- rbind(chrom.map, list("new" = paste0("chrArt",chrom_id), "old" = paste(as.character(chrom_sizes[c(j:i)]$V1), collapse = ";"), width = used_length ))
    cat(paste0(">chrArt",chrom_id), file =paste0(reference_path, db, "_concatinated.fa"), sep = "\n", append = T) 
    cat(paste(as.character(fa_file[chrom_sizes[c(j:i)]$V1]), collapse = linker), file = paste0(reference_path, db, "_concatinated.fa"), sep = "\n", append = T)
    j = i + 1 ## moving on to the next chromosome(that is not included yet)
    print(j)
    chrom_id <- chrom_id + 1  ## starting the new artificial chromosome
}

##save the mapping: 
write.table(chrom.map, paste0(reference_path, db, "_chrom_mapping.tsv"), sep = "\t", row.names = F, quote = F) 
print("Done")

