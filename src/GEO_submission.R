source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library("data.table")
library("gtools")

wd=data_dir
setwd(wd)

#sampleAnnot=fread("/data/groups/lab_bock/jklughammer/gitRepos/RefFreeDMA/meta/sampleAnnotation_combined.csv")
#sampleAnnot[,results:=ifelse(Organism=="Human","Hum_full_9",ifelse(Organism=="Bovine","Cow_full_6","Carp_full_3")),]
#sampleAnnot[,ref_genome:=ifelse(Organism=="Human","hg19",ifelse(Organism=="Bovine","bosTau6","carp_ena")),]
#sampleAnnot[,ded_genome:="toSelf_filtered_0.08mm_final_concat",]


GEO_dir="forGEO/"
processed_data_dir="processed_data"
raw_dir="raw_data"
system(paste0("mkdir -p ",GEO_dir,processed_data_dir))
system(paste0("mkdir -p ",GEO_dir,raw_dir))

ref_genome <- "toSelf_filtered_0.08mm_final_concat"
## checking if all the raw bam files are named according to last changes:
for (i in 1:nrow(stats_annot)){
    bam_file <- file.path("raw_data", stats_annot[i]$species, stats_annot[i]$Sample_Name, paste0(stats_annot[i]$Sample_Name, ".bam"))
    if(!file.exists(bam_file)) print(stats_annot[i]$Sample_Name)
}

### creating a sample annotation table(no file copying):
sample_annotation = data.table()
for (i in 1:nrow(stats_annot)){
    sample_name=stats_annot[i]$Sample_Name
    sample_annotation <- rbind(sample_annotation, data.frame(sample = sample_name, species = stats_annot[i]$species, raw_file=paste0(sample_name,"_unmapped.bam"),
  processed_file=paste0(sample_name,"_cpgMeth.bed"), 
   concat_fasta = paste0(stats_annot[i]$species, "_consensus_reference.fa")))
    }
sample_annotation <- unique(sample_annotation)
sample_annotation <- left_join(sample_annotation, stats_annot[, c("Sample_Name", "Tissue", "ncbi_name", "sex", "replicate", "age")], by = c("sample" = "Sample_Name"))

##
sample_annotation$title <- sample_annotation$sample
sample_annotation$molecule <- "genomicDNA"

write.table(sample_annotation[,c("sample", "title", "Tissue", "ncbi_name","sex", "replicate", "age", "molecule", "processed_file", "concat_fasta", "raw_file")],paste0(GEO_dir,"/sampleAnnotation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    

for (i in 2198:nrow(stats_annot)){
  sample_name=stats_annot[i]$Sample_Name
  print(sample_name)
  header=fread(file.path(processed_dir, stats_annot[i]$species, ref_genome, sample_name,"biseqMethcalling", "header_for_RRBS_cpgMethylation_file.txt"))
  
    meth_file <- file.path(processed_dir, stats_annot[i]$species, ref_genome, sample_name,"biseqMethcalling", paste0("RRBS_cpgMethylation_", sample_name, ".bed"))
   
   if(!file.exists(meth_file)){
       meth_file <- paste0(meth_file, ".gz")
   }
    
    bam_file <- file.path(processed_dir, stats_annot[i]$species,"unmapped_bam",paste0(sample_name,".bam"))
    
    fasta_file <- file.path(processed_dir,stats_annot[i]$species, "reduced", "consensus",ref_genome, paste0(ref_genome, ".fa"))
    
    if(!(file.exists(fasta_file) & file.exists(bam_file) & file.exists(meth_file))){
        print(paste0("one of the files missing for ", sample_name))
    }else{
    #meth file    
    meth_ded <- fread(meth_file)
    setnames(meth_ded,names(meth_ded),names(header))
write.table(meth_ded,paste0(GEO_dir,processed_data_dir,"/",sample_name,"_cpgMeth.bed"),sep="\t",quote=FALSE,row.names=FALSE)
   
    #bam file
   # system(paste0("cp ",bam_file," ",GEO_dir,"/",raw_dir,"/",sample_name,"_unmapped.bam"))
  
     #fasta file
    fasta_file_dest <- file.path(GEO_dir, processed_data_dir, paste0(stats_annot[i]$species, "_consensus_reference.fa"))
    if(!file.exists(fasta_file_dest)){
        system(paste0("cp ",fasta_file, " ", fasta_file_dest))
    }
     
   }
}

raw_md5=system(paste0("cd ", GEO_dir,"/",raw_dir," ; md5sum ","*" ),intern=TRUE)
proc_md5_fa=system(paste0("cd ", GEO_dir,"/",processed_data_dir," ; md5sum ","*.fa" ),intern=TRUE)

spl=unlist(strsplit(raw_md5, " "))
idx=seq(1,length(spl),3)
raw_md5_dt=data.table(md5_raw=spl[idx],raw_file=spl[idx+2])

spl=unlist(strsplit(proc_md5_fa, " "))
idx=seq(1,length(spl),3)
proc_md5_fa_dt=data.table(md5_proc=spl[idx],processed_file=spl[idx+2])

write.table(raw_md5_dt,paste0(GEO_dir,"/raw_files_md5.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(proc_md5_fa_dt,paste0(GEO_dir,"/processed_files_fa_md5.tsv"),sep="\t",quote=FALSE,row.names=FALSE)


proc_md5_bed=system(paste0("cd ", GEO_dir,"/",processed_data_dir," ; md5sum ","*.bed" ),intern=TRUE)

spl=unlist(strsplit(proc_md5_bed, " "))
idx=seq(1,length(spl),3)
proc_md5_bed_dt=data.table(md5_proc=spl[idx],processed_file=spl[idx+2])


raw_md5=system(paste0("cd raw_data ; md5sum ","*/*/*.bam" ),intern=TRUE)
proc_md5_fa=system(paste0("cd ", GEO_dir,"/",processed_data_dir," ; md5sum ","*.fa" ),intern=TRUE)