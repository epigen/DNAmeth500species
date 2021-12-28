library("data.table")
library("gtools")

wd="/data/scratch/lab_bock/jklughammer/projects/RefFreeDMA/"
setwd(wd)

sampleAnnot=fread("/data/groups/lab_bock/jklughammer/gitRepos/RefFreeDMA/meta/sampleAnnotation_combined.csv")
sampleAnnot[,results:=ifelse(Organism=="Human","Hum_full_9",ifelse(Organism=="Bovine","Cow_full_6","Carp_full_3")),]
sampleAnnot[,ref_genome:=ifelse(Organism=="Human","hg19",ifelse(Organism=="Bovine","bosTau6","carp_ena")),]
sampleAnnot[,ded_genome:="toSelf_filtered_0.08mm_final_concat",]


GEO_dir="forGEO/"
processed_dir="processed_data"
raw_dir="raw_data"
system(paste0("mkdir -p ",GEO_dir,processed_dir))
system(paste0("mkdir -p ",GEO_dir,raw_dir))


for (i in 1:nrow(sampleAnnot)){

  sample_name=sampleAnnot[i]$Sample_Name
  print(sample_name)
  header=fread(system(paste0("ls ",sampleAnnot[i]$results,"/",sampleAnnot[i]$ref_genome,"/*__",sample_name,"/biseqMethcalling/header_for_RRBS_cpgMethylation_file.txt"),intern=TRUE))
  meth_ref=fread(system(paste0("ls ",sampleAnnot[i]$results,"/",sampleAnnot[i]$ref_genome,"/*__",sample_name,"/biseqMethcalling/RRBS_cpgMethylation*__",sample_name,".bed"),intern=TRUE))
  meth_ded=fread(system(paste0("ls ",sampleAnnot[i]$results,"/",sampleAnnot[i]$ded_genome,"/*__",sample_name,"/biseqMethcalling/RRBS_cpgMethylation*__",sample_name,".bed"),intern=TRUE))
  
  setnames(meth_ref,names(meth_ref),names(header))
  setnames(meth_ded,names(meth_ded),names(header))
  
  write.table(meth_ref,paste0(GEO_dir,processed_dir,"/",sample_name,"_cpgMeth_ref.bed"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(meth_ded,paste0(GEO_dir,processed_dir,"/",sample_name,"_cpgMeth_ded.bed"),sep="\t",quote=FALSE,row.names=FALSE)

  system(paste0("cp ",sampleAnnot[i]$results,"/unmapped_bam/*",sample_name,".bam ",GEO_dir,"/",raw_dir,"/",sample_name,"_unmapped.bam"))

  sampleAnnot[i,raw_file:=paste0(sample_name,"_unmapped.bam")]
  sampleAnnot[i,processed_file_ref:=paste0(sample_name,"_cpgMeth_ref.bed")]
  sampleAnnot[i,processed_file_ded:=paste0(sample_name,"_cpgMeth_ded.bed")]
}

write.table(sampleAnnot[mixedorder(Sample_Name)],paste0(GEO_dir,"/sampleAnnotation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

raw_md5=system(paste0("cd ", GEO_dir,"/",raw_dir," ; md5sum ","*" ),intern=TRUE)
proc_md5=system(paste0("cd ", GEO_dir,"/",processed_dir," ; md5sum ","*" ),intern=TRUE)

spl=unlist(strsplit(raw_md5, " "))
idx=seq(1,length(spl),3)
raw_md5_dt=data.table(md5_raw=spl[idx],raw_file=spl[idx+2])

spl=unlist(strsplit(proc_md5, " "))
idx=seq(1,length(spl),3)
proc_md5_dt=data.table(md5_proc=spl[idx],processed_file=spl[idx+2])

write.table(raw_md5_dt[mixedorder(raw_file)],paste0(GEO_dir,"/raw_files_md5.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(proc_md5_dt[mixedorder(processed_file)],paste0(GEO_dir,"/processed_files_md5.tsv"),sep="\t",quote=FALSE,row.names=FALSE)




