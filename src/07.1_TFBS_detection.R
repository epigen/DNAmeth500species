#!/usr/bin/env Rscript

#classical initiation
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))


args = commandArgs(trailingOnly=TRUE)
species=args[1]
print(species)

#specific libraries
library(Biostrings)
library(TFBSTools)
#library(JASPAR2014)
library(GenomicRanges)
library(genomeIntervals)

wd=file.path(processed_dir, species)
setwd(wd)

subdir=file.path(analysis_dir, "07_motifAnalysis/07.1_TFBS_detection_2020", species)
dir.create(subdir, recursive = T)

#load data
meth_data_mean=fread(paste0("toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/",species,"_mean_meth.tsv"))
ded_ref_concat=readDNAStringSet("reduced/consensus//toSelf_filtered_0.08mm_final_concat/toSelf_filtered_0.08mm_final_concat.fa")
ded_ref_concat_bed=fread("reduced/consensus/toSelf_filtered_0.08mm_final_concat.bed")
ded_ref_concat_gr=with(ded_ref_concat_bed,GRanges(seqnames = Rle(V1), IRanges(start=V2, end=V3),strand=Rle("*"),name=V4))

#select features (methylated/unmethylated)
meth_data_mean_long=melt(meth_data_mean,id.vars=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts"),measure=patterns(".cov", ".meth"),variable.factor=TRUE,variable.name="sample",value.name=c("cov","meth"),na.rm=TRUE)
sample_names=sub(".meth","",grep(".meth",names(meth_data_mean),value=TRUE))

meth_data_mean_long[,sample_name:=sample_names[sample],]
meth_data_mean_long[,species:=sub("_.*","",sample_name),]
meth_data_mean_long[,tissue:=sub(".*_","",sample_name),]
meth_data_mean_cond=meth_data_mean_long[,list(Nsamples=.N,mean_cov=mean(cov),min_cov=min(cov),max_cov=max(cov),mean_meth=mean(meth),min_meth=min(meth),max_meth=max(meth)),by=c("meta","score","consN","mean_diffNts","max_diffNts","min_diffNts")]
meth_data_mean_cond_red=meth_data_mean_cond[Nsamples>(max(Nsamples)*0.5)&mean_cov>10&mean_cov<1000]
meth_data_mean_cond_red[,category:=ifelse(min_meth>80,"meth",ifelse(max_meth<20,"unmeth","amb")),]

#cleanup
rm(meth_data_mean)

#scan for all TFs that are available in JASPAR
#opts <- list()
#opts[["tax_group"]] <- "vertebrates"
#opts[["collection"]] <- "CORE"
#PFMatrixList <- getMatrixSet(JASPAR2014, opts)


PFMatrixList <- readJASPARMatrix(fn = file.path(data_dir, "resources","JASPAR", "JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"))#, type = "all")

modify_set<-function(siteset){
  gff=as.data.table(TFBSTools::writeGFF3(siteset))
  if (NROW(gff)>0){
    siteset_gr=with(gff,GRanges(seqnames = Rle(seqname), IRanges(start=start, end=end),strand=Rle(strand)))
    ol=as.data.table(findOverlaps(siteset_gr,ded_ref_concat_gr,type="within"))
    siteset_annot=cbind(gff[ol$queryHits],ded_ref_concat_bed[ol$subjectHits])
    siteset_annot[,rel_start:=start-V2+1,]
    siteset_annot[,rel_end:=end-V2+1,] 
  }
  else{
    siteset_annot=gff
  }
  return(siteset_annot)
}

pwm_list=lapply(PFMatrixList,function(x)
{return(toPWM(x, type="log2probratio", pseudocounts=0.8,
              bg=c(A=0.25, C=0.25, G=0.25, T=0.25)))})

full_searchSeq <- function(pwmat){
    filepath <- paste0(subdir, "/", pwmat@name, "_profile_table.tsv")
    if(file.exists(filepath)){
        print(paste0(pwmat@name, " done"))
    }else{
        siteset = searchSeq(pwmat, ded_ref_concat, min.score='90%', strand='*')

  siteset_annot=modify_set(siteset)
    
  if(NROW(siteset_annot) > 0){
  siteset_annot[,c("TF","class","sequence") := tstrsplit(gsub("TF=|class=|sequence=","",attributes), ";")]  
  siteset_annot_red=merge(siteset_annot,meth_data_mean_cond_red,
                          by.x="V4",by.y="meta",all=FALSE)
  write.table(siteset_annot_red,
              paste0(subdir, "/", pwmat@name, "_profile_table.tsv"),sep="\t",
              quote=FALSE,row.names=FALSE) }
  else(print(paste0("NONE ", pwmat@name)))
    }
  }

lapply(pwm_list, full_searchSeq)

