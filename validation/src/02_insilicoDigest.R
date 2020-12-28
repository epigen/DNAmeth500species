## ---------------------------
##
## Insilico digestion analysis
##
## Simulation of the RRBS prototcol insilico, evaluation of the percentage of CpGs covered in regulatory elements of interest
##
## Authors: Johanna Klughammer, Daria Romanovskaia
##
## Date Created: 2020-12-17
##
##
## ---------------------------



## uploading needed libraries

suppressMessages(library(Biostrings))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer)) ## for reading - exporting


##parcing the arguments:
options <- commandArgs(trailingOnly = TRUE)  ## as an argument read in the row in the annotation matrix

idx = options[1]

annot_table <- read.csv("meta/validation_species_resources.csv", sep = ";", stringsAsFactors = FALSE)

print(paste0("working with ", annot_table[idx,]$genome_id))
print(annot_table[idx,])

genome_id = as.character(annot_table[idx,]$genome_id)
species = as.character(annot_table[idx,]$species)

##identifying where to save the output
output_dir <- file.path( "/scratch/lab_bock/shared/projects/compEpi/validation/insilicoDigest/model", genome_id)
resource_dir <- file.path("/scratch/lab_bock/shared/projects/compEpi/validation/insilicoDigest/data/tracks/", species)

dir.create(output_dir,  showWarnings=FALSE, recursive = TRUE)
print(output_dir)


library(annot_table[idx,]$BSlibrary, character.only = TRUE)
#library(BSgenome.Cmilii.UCSC.calMil1)
##extracting sequences from a chromosome
get_seq <- function(chr, start, end, genome){
    
    tryCatch({s <- subseq(genome[[chr]], start, end)}, 
             error = function(e) {
                 print(genome[[chr]])
                 print(start)
                 print(end)
             })
    return(as.character(s))
}

# Identify the MspI recongnition sites for each chromosomal entry
# Generate a dataframe with the length of MspI digested fragments
insilicoRRBS <- function(genome, output_dir_genome){
    
    ##redirecting output
    print(output_dir_genome)

    genome_set <- get(genome)
    
    if(!file.exists(paste0(output_dir_genome, "/uniqFragments_full.RData"))){

    mdf=data.frame()
    uniqFragments=data.table()
    
    ## looking for start positions of RRBS-sp enzymes
    for (pattern_seq in c("CCGG", "TCGA")){
      for (i in seq_along(genome_set)){
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
      temp_df<-data.frame(start=starts-3,end=ends+3,width=widths+6,chr=seqnames(genome_set)[i],genome=genome) #actually end = ends
          
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
      uniqFragments=rbindlist(list(uniqFragments,data.table(genome=genome, pattern = pattern_seq, chr=seqnames(genome_set)[i],
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



###  creating a GRanges obect with all the possible CpGs:
genome_set <- get(genome_id)

if(!file.exists(paste0(output_dir, "/CpGs_genomewide.RData"))){

  CpGs_table = data.table()

#collecting the CpGs from each chromosome:
  for (i in seq_along(genome_set)){
      if ("DNAStringSet"%in%is(genome_set[[i]])) {
                        print(length(genome_set[[i]]))
                                                    next}
      CpGs = matchPattern("CG", genome_set[[i]])
      if(length(CpGs) > 0){
      CpGs_df <- data.frame(chr = seqnames(genome_set)[[i]], start =  start(CpGs), end = end(CpGs), width = width(CpGs))
      CpGs_table <- rbind(CpGs_table, CpGs_df)
      }
      
  }

  print("coordinates of CpGs in the genome extracted")

  CpGs_GRanges <- makeGRangesFromDataFrame(CpGs_table)

  save(CpGs_GRanges,file=paste0(output_dir, "/CpGs_genomewide.RData"))
}else{
  print("coordinates of CpGs in the genome already extracted")
  load(paste0(output_dir, "/CpGs_genomewide.RData"))

}


### looking for overlaps:

##transforming the UnF to GRanges:
unF_granges <- makeGRangesFromDataFrame(uniqFragments, keep.extra.columns=TRUE)

## total overlap:
CpGs_in_RRBS <- GenomicRanges::intersect(unF_granges, CpGs_GRanges)

##creating the output file:
output_stat_file <- file.path(output_dir, paste0(as.character(annot_table[idx,]$ucsc_genome), "_overlap_stats_final.tsv"))

cat("type\tCpGs\tCpGs_in_RRBS", file=output_stat_file, sep="\n")

##saving first overlap:
cat(paste0(c("total", NROW(CpGs_GRanges), NROW(CpGs_in_RRBS)), sep = "\t"),file=output_stat_file, append = TRUE)

###annotations:

annot_files <- list.files(resource_dir)

for (file_path in annot_files){
  annot_id <- gsub("\\.bed", "", gsub(".*_","", file_path))
  print(annot_id)
  ##reading the annotation:
  annot <- import.bed(file.path(resource_dir, file_path))
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
  ## coord2
  #promoters_annot <- promoters(annot, upstream = 2000, downstream = 1000)

  #CpGs_in_promoters <- GenomicRanges::intersect(CpGs_GRanges, GRanges(promoters_annot), ignore.strand = TRUE)
  #CpGs_in_promoters_in_RRBS <- GenomicRanges::intersect(CpGs_in_promoters, unF_granges, ignore.strand = TRUE)

#  cat(paste0(c(paste0("\n",annot_id, "2000_1000"), NROW(CpGs_in_promoters), NROW(CpGs_in_promoters_in_RRBS)), sep = "\t"),
  #                                                           file=output_stat_file, append = TRUE)

  # coord 3
 # promoters_annot <- promoters(annot, upstream = 500, downstream = 200)

  #CpGs_in_promoters <- GenomicRanges::intersect(CpGs_GRanges, GRanges(promoters_annot), ignore.strand = TRUE)
  #CpGs_in_promoters_in_RRBS <- GenomicRanges::intersect(CpGs_in_promoters, unF_granges, ignore.strand = TRUE)

#  cat(paste0(c(paste0("\n",annot_id, "500_200"), NROW(CpGs_in_promoters), NROW(CpGs_in_promoters_in_RRBS)), sep = "\t"),
 #                                                            file=output_stat_file, append = TRUE)


}else{

  cat(paste0(c(paste0("\n",annot_id), NROW(CpGs_in_annot), NROW(CpGs_in_annot_in_RRBS)), sep = "\t"),
                                                             file=output_stat_file, append = TRUE)
}

}

print("Done")

