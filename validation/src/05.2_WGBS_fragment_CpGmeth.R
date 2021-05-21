## ---------------------------
##
## Modelling the methylation prediction in WGBS
##
## Part 1: data preparation.
## The script parses through pre-downloaded genomes and CpG-calls in the bed format and generates 
## the freffreedma like fragments and the methylation table for each of them 
## The genomes should be pre-installed in the BSgenome format (see download_WGBS.ipynb)
## Data from GEO is predownloaded and parsed/reformatted (get_mean_meth.ipynb)
##
## Authors: dromanovskaia
##
## Date Created: 2021-05-12
## ---------------------------

## uploading needed libraries

suppressMessages(library(Biostrings))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer)) ## for reading - exporting
suppressMessages(library(BSgenome))
suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))


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

annot_full <- fread("../meta/WGBS_prediction_selection.tsv")

##reading which one to use from the annotation
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
idx <- as.numeric(args[[1]]) - 1
print(idx)
print(annot_full$BSlibrary)
#upload the library with the genome
library(as.character(annot_full$BSlibrary[idx]), character.only = TRUE)

## setup the paths
path_to_folder <- paste0("../../resources/WGBS_public/", paste(strsplit(annot_full$Species[[idx]], " ")[[1]], collapse="_"), "/")
outfolder <- paste0("../validation/WGBS_public/prediction/", paste(strsplit(annot_full$Species[[idx]], " ")[[1]], collapse="_"), "/",annot_full$assembly[[idx]] ,"/")
dir.create(outfolder, recursive = T)

genome_sequence <- get(annot_full$genomeId[[idx]])

## splitting the genome into tiles
tiles_50bp <- unlist(tileGenome(seqinfo(genome_sequence), tilewidth = 50))
print(head(tiles_50bp))
## setting up the sequence names
tiles_50bp_df <- as.data.table(tiles_50bp)
tiles_50bp_df[, name:=paste0(seqnames, "_", start, "_", end),]


if(!file.exists(paste0(outfolder, "/",annot_full$assembly[idx], "_fragments.fa"))){
## extracting the sequences
seq_mmm9 = getSeq(genome_sequence, tiles_50bp)
## for the seq:

names(seq_mmm9) = tiles_50bp_df$name

## reading the fasta file
writeXStringSet(seq_mmm9, paste0(outfolder, "/",annot_full$assembly[idx], "_fragments.fa"))
print("sequences saved")
}
## upload the per sample annotation table
annot_samples <- fread(paste0(path_to_folder, "WGBS_annot.csv"))

df_list <- mclapply(seq_along(annot_samples$path_unif), 
                       function(x) parse_CpG_calls(annot_samples$path_unif[[x]], annot_samples$tissue[[x]], annot_samples$replica[[x]], tiles_50bp))
                       
saveRDS(df_list,paste0(annot_full$assembly[idx], 'error.RDS'))

df <- Reduce(full_join, df_list)

head(df)

write.table(df, paste0(outfolder, "mean_meth_per_fragment.tsv"), sep = ";", quote = F, row.names = F)