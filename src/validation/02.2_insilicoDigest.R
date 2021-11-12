## ---------------------------
##
## Insilico digestion analysis
##
## Simulation of the RRBS prototcol insilico, evaluation of the percentage of CpGs covered in regulatory elements of interest
##
## Authors: Johanna Klughammer, Daria Romanovskaia
##
## Date Created: 2020-12-17, upd. 2021-10-18
##
##
## ---------------------------

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))


## uploading needed libraries

suppressMessages(library(Biostrings))
suppressMessages(library(rtracklayer))## for reading - exporting


##parcing the arguments:
options <- commandArgs(trailingOnly = TRUE)
genome_id = options[1]

output_dir <- file.path(analysis_dir,"validation", "02_insilico_digest", genome_id)
dir.create(output_dir,  showWarnings=FALSE, recursive = TRUE)
print(output_dir)

resource_dir <- file.path(data_dir, "resources","reference_genomes", genome_id)


### prepare the conversion into old coordinates
###going back to old coordinates:
chrom_mapping <-  fread(file.path(data_dir, "resources", "reference_genomes", genome_id, paste0(genome_id,"_chrom_mapping.tsv")))

##filter out tiny chromosomes
chrom_mapping[,n:=length(unlist(strsplit(old,";"))),by=1:nrow(chrom_mapping)]

if(max(chrom_mapping$n) > 10) chrom_mapping=chrom_mapping[n<quantile(chrom_mapping$n,0.8) & n < 3000]

###create long chromosome mapping + sizes
chrom_mapping_long=chrom_mapping[,.(old = as.character(unlist(tstrsplit(old, ";", type.convert = TRUE)))), by = "new"]

chrom_sizes <- fread(file.path(data_dir, "resources", "reference_genomes", genome_id, paste0(genome_id,".chrom.sizes")))

chrom_mapping_long[,old:=ifelse(!is.na(as.numeric(old)),chrom_sizes[as.numeric(old)]$V1,old),] # this is to use it as index "old" as index if necessary

chroms=merge(chrom_mapping_long, chrom_sizes, by.x='old',by.y='V1',sort=FALSE)

chroms[,new_pos:=ifelse(.N!=1,0,c(0,cumsum(V2 + 100)[1:length(V2)-1])),by=new]
chroms[,new_pos:=c(0,cumsum(V2 + 100)[1:length(V2)-1]),by=new]


unconcat_coordinates <- function(chroms, bam_dt){
    print(NROW(bam_dt))
    bam_dt=bam_dt[!is.na(chr)]
    bam_dt=bam_dt[chr%in%chroms$new]
    print(NROW(bam_dt))
    bam_dt[,name:=paste(chr, start, end, sep = "_"), by = row.names(bam_dt)]
    
#find correct chrom and pos
    bam_dt_cor = data.table()
    for(chrArt in unique(bam_dt$chr)){
        print(chrArt)
        bam_dt_sub = bam_dt[chr == chrArt]
        chroms_sub = chroms[new == chrArt]
        print(NROW(bam_dt_sub))
        bam_dt_cor_sub = merge(bam_dt_sub, chroms_sub, by.x='chr',by.y='new',allow.cartesian=TRUE)
        bam_dt_cor_sub[,dist:=start-new_pos,]
        bam_dt_cor_sub=bam_dt_cor_sub[dist>0]
        bam_dt_cor_sub[,min_dist:=min(dist),by='name']
        bam_dt_cor_sub=bam_dt_cor_sub[dist==min_dist]
        ##combinig for each art chromosome
        bam_dt_cor = rbind(bam_dt_cor, bam_dt_cor_sub)
    }
   
    #bam_dt_cor=merge(bam_dt,chroms, by.x='chr',by.y='new',allow.cartesian=TRUE)
    
    bam_dt_cor[,chr_cor:=old,]
    bam_dt_cor[,start_cor:=start-new_pos,]
    bam_dt_cor[,end_cor:=end - new_pos,]   

    bam_dt_cor <- bam_dt_cor %>% 
      dplyr::rename(
        chr_conc = chr,
        start_conc = start,
        end_conc = end)
    
    bam_dt_cor <- bam_dt_cor %>% 
      dplyr::rename(
        chr = chr_cor,
        start = start_cor,
        end = end_cor)
    
    return(bam_dt_cor)
}



# Identify the MspI recongnition sites for each chromosomal entry
# Generate a dataframe with the length of MspI digested fragments
insilicoRRBS <- function(genome, output_dir_genome){
    
    #genome_set <- get(genome)
    ##reading in the full data
    resource_dir <- file.path(data_dir, "resources", "reference_genomes", genome)
    if(file.exists(file.path(resource_dir, paste0(genome, "_concatinated.fa")))){
        genome_set <- readDNAStringSet(file.path(resource_dir, paste0(genome, "_concatinated.fa")))
    }else{
        print("no fasta file avaliable!")
        stop()
    }
    
    
    if(!file.exists(paste0(output_dir_genome, "/uniqFragments_full.RData"))){

    mdf=data.frame()
    uniqFragments=data.table()
    
    ## looking for start positions of RRBS-sp enzymes
    for (pattern_seq in c("CCGG", "TCGA")){
    print(paste0("running on ", pattern_seq))
    for (i in seq_along(genome_set)){
     #for(i in c(1,2)){
              if ("DNAStringSet"%in%is(genome_set[[i]])) {
                        print(length(genome_set[[i]]))
                                                    next}
      print(paste0("running on chromosome", i))    
      ##matching the pattern to the genome
      m <- matchPattern(pattern_seq, genome_set[[i]])
      ##getting everything, EXCEPT the pattern      
      frags = gaps(m)

      if(length(frags) > 0){      
      ##extracting geomic coordinates
      starts<-start(frags)
      ends<-end(frags)
      widths<-width(frags)

      ##generating a df and performing coordinate shifts to include that reads to start with CGG/TCG and end with the reverse complement
      temp_df<-data.frame(start=starts-3,end=ends+3,width=widths+6,chr=names(genome_set)[i],genome=genome) #actually end = ends
          
      ##as the first becomes -2 => shifting to 0
      temp_df$start<-replace(temp_df$start, temp_df$start == -2, 0)
      temp_df$name <- pattern_seq
      temp_df<-temp_df[c("genome","chr","start","end","width", "name")]
      mdf<-rbind(mdf,temp_df)
  
      #get coordinates of the start and the end for the generated sequences:
      
      #50 - size-selection(shortest you get via the library prep), but: including the adapters, which we don't account for
      #500 - upper threshold for the width of the fragment we are interested in
      sub=frags[width(frags)>=50]
      if(length(sub) > 0){
      #fragment ends
      tail_start=width(sub)-49
      tail_end=width(sub)
      tail_end_pos=end(sub)+3
      tail_start_pos=end(sub)-49+3
  
      ##fragment starts
      front_end_pos=start(sub)+49-3
      front_start_pos=start(sub)-3
      
      ##summary dataframe    
      uniqFragments=rbindlist(list(uniqFragments,data.table(genome=genome, pattern = pattern_seq, chr=names(genome_set)[i],
                                                            maxlength=length(genome_set[[i]]),
                                                            start=c(front_start_pos,tail_start_pos),
                                                            end=c(front_end_pos,tail_end_pos),
                                                            width=c(width(sub),width(sub)),
                                                            origin=c(rep("front",length(front_start_pos)),rep("tail",length(tail_start_pos))))))}
      }
    }
    }
    save(uniqFragments, file = paste0(output_dir_genome, "/uniqFragments_full.RData"))
    }else{

      print("already extracted pattern coordinates, loading...")
      load(paste0(output_dir_genome, "/uniqFragments_full.RData"))
    }
    #uniqFragments[,dupl:=duplicated(seq),by=list(genome, pattern)]
    ## adjusting the start of the first fragment to 1:
    print("adjusting chromosome borders...")
    uniqFragments[uniqFragments$start <= 0, ]$start <- 1
    
    ##adjusting the end of the chromosomes after correction to the chr length
    if(NROW(uniqFragments[uniqFragments$end>uniqFragments$maxlength, ])> 0){
    uniqFragments[uniqFragments$end>uniqFragments$maxlength, ]$end <- uniqFragments[uniqFragments$end>uniqFragments$maxlength, ]$maxlength
    }
    ##filtering by width
    print("filtering the sequences: ")
    print(NROW(uniqFragments))
    uniqFragments <- uniqFragments[uniqFragments$width > 50 & uniqFragments$width < 1000,]
    print(NROW(uniqFragments))

    ##extracting the sequence:
   # print("extracting sequences...")

  #  uniqFragments[, seq_exact :=as.character(subseq(genome_set[[chr]], start, end)), by = 1:nrow(uniqFragments)]
    
    ##!! WIDTH IS NOT ACTUAL, (NOT NEEDED)

    ## we are only interested in the sequences, 
    ##saving output
    
   
    print("saving length distribution...")
    ##saving stats of fragment length:
    mdf=as.data.table(mdf)
    print(head(mdf))
    
    print(unique(mdf$name))
    
    count=mdf[,.N,by=c("width","genome", "name")]
    
    ggplot(count[width>50&width<1000],aes(x=width,y=N, color = name)) + 
                                geom_line() +
                                facet_wrap(~genome,scale="free_y") +
                                    xlab("fragment length") + theme_bw()
    
    ggsave(file.path(output_dir_genome, "fragmentLengths.pdf"),width=8,height=5)
    
    write.csv(mdf, file.path(output_dir_genome, "fragmentLengths.csv"))
    
     #saving Rdata
    save(uniqFragments,file=paste0(output_dir_genome, "/uniqFragments.RData"))
    
    return(uniqFragments)

}

### creating the insilico RRBS simulation

if(file.exists(paste0(output_dir, "/uniqFragments.RData"))){
    print("insilico simulation already performed")
    load(paste0(output_dir, "/uniqFragments.RData"))
}else{
    print("running insilico RRBS simulation")
    uniqFragments <- insilicoRRBS(genome_id, output_dir)
}



if(!file.exists(paste0(output_dir, "/uniqFragments_gRanges.RDS"))){
    ##transforming the coordinates
    uniqFragments_coord <- unconcat_coordinates(chroms, uniqFragments)
    ##transforming the UnF to GRanges:
    unF_granges <- makeGRangesFromDataFrame(uniqFragments_coord, keep.extra.columns=TRUE)
    saveRDS(unF_granges,file=paste0(output_dir, "/uniqFragments_gRanges.RDS"))
}else unF_granges <- readRDS(paste0(output_dir, "/uniqFragments_gRanges.RDS"))

generate_CpG_map <- function(genome_id){
###  creating a GRanges obect with all the possible CpGs:
    genome_set <- readDNAStringSet(file.path(data_dir, "resources","reference_genomes", genome_id, paste0(genome_id, "_concatinated.fa")))

  CpGs_table = data.table()

#collecting the CpGs from each chromosome:
  for (i in seq_along(genome_set)){
   # for(i in c(1,2)){
    if ("DNAStringSet"%in%is(genome_set[[i]])) {
                        print(length(genome_set[[i]]))
                                                    next}
      CpGs = matchPattern("CG", genome_set[[i]])
      if(length(CpGs) > 0){
      CpGs_df <- data.frame(chr = names(genome_set)[[i]], start =  start(CpGs), end = end(CpGs), width = width(CpGs))
      CpGs_table <- rbind(CpGs_table, CpGs_df)
      }
      
  }

  print("coordinates of CpGs in the genome extracted")
    return(CpGs_table)
 # 
   # return(CpGs_GRanges)
}
    


if(!file.exists(paste0(output_dir, "/CpGs_genomewide.RDS"))){
    CpGs_table <- generate_CpG_map(genome_id)
    saveRDS(object = CpGs_table,file = paste0(output_dir, "/CpGs_genomewide.RDS"))
}else{
  print("coordinates of CpGs in the genome already extracted")
  CpGs_table <- readRDS(paste0(output_dir, "/CpGs_genomewide.RDS"))
}

if(!file.exists(paste0(output_dir, "/CpGs_genomewide_granges.RDS"))){
    CpGs_table_coord <- unconcat_coordinates(chroms, CpGs_table)
    CpGs_GRanges <- makeGRangesFromDataFrame(CpGs_table_coord)
    saveRDS(object = CpGs_GRanges,file = paste0(output_dir, "/CpGs_genomewide_granges.RDS"))
}else CpGs_GRanges <- readRDS(paste0(output_dir, "/CpGs_genomewide_granges.RDS"))


### looking for overlaps:
## total overlap:
CpGs_in_RRBS <- GenomicRanges::intersect(unF_granges, CpGs_GRanges)

##creating the output file:
output_stat_file <- file.path(output_dir, paste0(genome_id, "_overlap_stats_final.tsv"))

cat("type\tCpGs\tCpGs_in_RRBS", file=output_stat_file, sep="\n")

##saving first overlap:
cat(paste0(c("total", NROW(CpGs_GRanges), NROW(CpGs_in_RRBS)), sep = "\t"),
            file=output_stat_file, append = TRUE)

###annotations:

annot_files <- list.files(file.path(resource_dir, "tracks"))

for (file_path in annot_files){
  annot_id <- gsub("\\.bed", "", gsub(".*_","", file_path))
  print(annot_id)
  ##reading the annotation:
  annot <- import.bed(file.path(resource_dir, "tracks", file_path))
  ##CpGs within the annotation
  CpGs_in_annot <- GenomicRanges::intersect(CpGs_GRanges, GRanges(annot), ignore.strand = TRUE)
  ##which of them are in RRBS
  CpGs_in_annot_in_RRBS <- GenomicRanges::intersect(CpGs_in_annot, unF_granges, ignore.strand = TRUE)

  ##saving output numbers:
  if (annot_id %in% c("refSeqComposite", "ensGene", "xenoRefGene", "refGene", "genscan", "lampreyGene")){

  annot_id <- "transcripts"
  cat(paste0(c(paste0("\n",annot_id), NROW(CpGs_in_annot), NROW(CpGs_in_annot_in_RRBS)), sep = "\t"),
                                                             file=output_stat_file, append = TRUE)
  annot_id <- "promoters"
  print(annot_id)
  ##extracting promoters
  promoters_annot <- promoters(annot, upstream = 1000, downstream = 500)

  CpGs_in_promoters <- GenomicRanges::intersect(CpGs_GRanges, GRanges(promoters_annot), ignore.strand = TRUE)
  CpGs_in_promoters_in_RRBS <- GenomicRanges::intersect(CpGs_in_promoters, unF_granges, ignore.strand = TRUE)

  cat(paste0(c(paste0("\n",annot_id, "1000_500"), NROW(CpGs_in_promoters), NROW(CpGs_in_promoters_in_RRBS)), sep = "\t"),
                                                             file=output_stat_file, append = TRUE)


}else{

  cat(paste0(c(paste0("\n",annot_id), NROW(CpGs_in_annot), NROW(CpGs_in_annot_in_RRBS)), sep = "\t"),
                                                             file=output_stat_file, append = TRUE)
}

}

print("Done")

