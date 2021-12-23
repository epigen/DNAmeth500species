## ---------------------------
##
## Modelling the methylation prediction in WGBS
##
## Part 1: data preparation.
## The script parses through pre-downloaded genomes and CpG-calls in the bed format and generates 
## the freffreedma like fragments and the methylation table for each of them 
## The genomes should be pre-downlaoded/linked from previous ecxperiments (see -03.1_download_WGBS.ipynb)
## Data from GEO is predownloaded (03.1) and parsed/reformatted (03.2_process_WGBS_data.ipynb)
##
## Authors: dromanovskaia
##
## Date Created: 2021-05-12, upd. 2021-11-04
## ---------------------------

## uploading needed libraries
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

suppressMessages(library(Biostrings))
suppressMessages(library(rtracklayer)) ## for reading - exporting
suppressMessages(library(BSgenome))
suppressMessages(library(GenomicRanges))


##############################################
## PARSE CpG DATA IN CALLS FORMAT           ##
## path to read the table from              ##
##tissue name                               ##
## replica                                  ##
## the bed file of the genome split in tiles##
##############################################

parse_CpG_calls <- function(path, tissue, replica, bed_tiles){
    
    print(path)

    df <- fread(path)
    df$id <- as.numeric(row.names(df)) ## row ids that we can match the subject hits
    print(head(df))
    
    bed_tiles_df <- as.data.table(bed_tiles) ## df with tiles and their names (for the end we just need the tile name and the sequence level there)
    bed_tiles_df$id_tile <- as.numeric(row.names(bed_tiles_df)) ## row numbers to match the CpGs
    bed_tiles_df[, name:=paste0(seqnames, "_", start, "_", end),] ## names of the sequences
    
    
    ##checking if the chromosome identification matches
    if (length(grep("chr", bed_tiles_df$seqnames[1])) != length(grep("chr", df$chr[1]))) {
    if(length(grep("chr", bed_tiles_df$seqnames[1])) == 1) df[,chr:=paste0("chr", chr),]
    else df[,chr:=gsub("chr","", chr),]}
    print(setdiff(unique(df$chr), unique(tiles_50bp_df$seqname)))
    
    if("strand" %in% colnames(df)){
    df <- df[,-c("strand")] ## otherwise you can't make a GRanges object
    }
    
    ##creating the GRanges
    df_granges <- makeGRangesFromDataFrame(df,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("chr"),
                         start.field=c("start"),
                         end.field=c( "end"),
                         starts.in.df.are.0based=FALSE)
    ##overlapping the regions (aka fake fragments) with the CpG coordinates
    ovlaps <- as.data.table(findOverlaps(bed_tiles, df_granges))
    
    
    ## joining by the number of regions in both datasets
    ovlaps_annot <- left_join(left_join(ovlaps, df, by = c("subjectHits" = "id")), bed_tiles_df, by = c("queryHits" = "id_tile"))
    
    if (max(ovlaps_annot$perc_meth_CpG)<=1) ovlaps_annot[,perc_meth_CpG:=100*perc_meth_CpG,] ##so we always have the same range
    
    summary_df <- ovlaps_annot %>% group_by(name) %>% 
            summarize(cov = mean(coverage), mean_meth = weighted.mean(perc_meth_CpG, coverage, na.rm = TRUE))
    colnames(summary_df)[[2]] <- paste0("cov", "_", tissue, "_", replica)
    colnames(summary_df)[[3]] <- paste0("mean_meth", "_", tissue, "_", replica)
    return(summary_df)
}


##########MAIN##########

##reading in the info about the species
args <- commandArgs(trailingOnly=TRUE)
print(args)

##reading which one to use from the annotation
species <- args[[1]]
genomeid <- args[[2]]


## setup the paths
path_to_folder <- file.path(data_dir, "resources", "WGBS_public", species)

## output
outfolder <- file.path(analysis_dir, "validation", "03_WGBS","03.3_fragments", species)
dir.create(outfolder, recursive = T)
setwd(outfolder)

### reading in the fasta file:
fasta_file = file.path(path_to_folder, "refgenome", paste0(genomeid, ".fa"))
genome_fa <- readDNAStringSet(fasta_file)


## splitting the genome into tiles
tiles_50bp <- unlist(tileGenome(seqinfo(genome_fa), tilewidth = 50))

#tiles_50bp <- unlist(tileGenome(seqinfo(genome_sequence), tilewidth = 50))
print(head(tiles_50bp))

## setting up the sequence names
tiles_50bp_df <- as.data.table(tiles_50bp)
tiles_50bp_df[, name:=paste0(seqnames, "_", start, "_", end),]



if(!file.exists(paste0(genomeid, "_fragments.fa"))){
## extracting the sequences
seq = getSeq(genome_fa, tiles_50bp)
## for the seq:
names(seq) = tiles_50bp_df$name
## saving the fasta file
writeXStringSet(seq, paste0(genomeid, "_fragments.fa"))
print("sequences saved")
}

## upload the per sample annotation table
annotation_file = file.path(path_to_folder, "WGBS_annot.csv")
annot_samples <- fread(annotation_file)

##temp path issue fix - REMOVE
annot_samples$path_unif <- gsub("/nobackup/lab_bock/projects/DNAmeth500species", data_dir, annot_samples$path_unif) 

df_list <- mclapply(seq_along(annot_samples$path_unif), 
                       function(x) parse_CpG_calls(annot_samples$path_unif[[x]], 
                                                   annot_samples$tissue[[x]], annot_samples$replica[[x]], tiles_50bp))
saveRDS(file = "df_list.RDS",object = df_list)                  
df <- Reduce(full_join, df_list)
write.table(df, "mean_meth_per_fragment.tsv", sep = ";", quote = F, row.names = F)
