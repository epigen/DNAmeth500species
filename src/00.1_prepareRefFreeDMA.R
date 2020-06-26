source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=processed_dir

dir.create(wd,recursive=TRUE)
setwd(wd)
dir.create("00_RefFreeDMA_log")

#set configurations for RefFreeDMA
config="
tool_path=/data/groups/lab_bock/jklughammer/resources/RefFreeDMA_resources/tools/
picard_path=$tool_path/picard-tools_1.118/
trim_galore_path=$tool_path/trim_galore_0.3.3/
cutadapt_path=$tool_path/cutadapt_1.8.3/
bowtie2_path=$tool_path/bowtie_2.2.4/bin/
bsmap_path=$tool_path/bsmap_2.90/
samtools_path=$tool_path/samtools_1.2/bin/
bedtools_path=$tool_path/bedtools_2.26.0/bin/
bwa_path=$tool_path/bwa_0.7.8/bin/
bwameth_path=$tool_path/bwa-meth-0.10/
decon_reference=/data/groups/lab_bock/jklughammer/resources/RefFreeDMA_resources/decon_reference/bacterial_extracted_add
cross_genome_fa=-
sample_annotation=$working_dir/meta/sampleAnnotation.tsv
compCol=comp
groupsCol=Tissue
parallel=TRUE
nonCpG=TRUE
unconv_tag=\"_uc\"
decon=TRUE
restrictionSites=\"CGG|CGA\"
wait_time=10
nProcesses=4
nameSeparator=\"#\"
maxReadLen=51
maxSamples=`ls $working_dir/unmapped_bam/*.bam|wc -l`
cLimit=0.05
mapToSelf_filter=0.08
consensus_dist=0.05
crossMap_mismatchRate=0.2
nTopDiffMeth=500
"

##Functions
getBamPath=function(FlowCell,lane,Abbreviation,alternative){

  pattern=paste0(FlowCell,"_.*_",lane,"\\#",Abbreviation,"(?:_S[0-9]+)?.bam$")
  
  path=list.files(dir(BSF_dir,pattern=paste0("BSF_[0]*",FlowCell,"_.*"),full.names=TRUE),pattern=pattern,full.names=TRUE,recursive = TRUE)
  if (length(path)==0){
    if(alternative!="skip"|!is.na(FlowCell)){
    path=list.files(dir(BSF_dir,pattern=paste0("BSF_[0]*",FlowCell,"_.*"),full.names=TRUE),pattern=paste0(FlowCell,"_.*_",lane,"\\#",alternative,"(?:_S[0-9]+)?.bam$"),full.names=TRUE,recursive = TRUE)
    }else{return("NA");stop("Skipping!")}
  }
  if (length(path)==0){message(paste0(Abbreviation," bam not found!"))}
  return(path)
}

makeLink=function(bam_path,run_species,Abbreviation){
  link=paste0(run_species,"/unmapped_bam/",Abbreviation,".bam")
  command=paste0("ln -f -s ",bam_path," ",link)
  if(grepl("NA ",command)){
    entry=readline("INFO: Link failed. Press ENTER to continue.")
  }
  system(command)
  print(command)
}

#fill all blanks with NA
sampleAnnot[sampleAnnot==""]=NA

#reconstruct speciess abbreviation and replicate
sampleAnnot[,abbreviation_sp:=unlist(lapply(Abbreviation,function(x){unlist(strsplit(x,"_"))[1]}))]
sampleAnnot[,replicate:=unlist(lapply(Abbreviation,function(x){unlist(strsplit(x,"_"))[2]}))]
sampleAnnot[,abbreviation_ti:=unlist(lapply(Abbreviation,function(x){unlist(strsplit(x,"_"))[3]}))]

#check for which species all samples have been sequenced
sampleAnnot[,complete:=ifelse(all(!is.na(FlowCell))&all(!is.na(lane)),TRUE,FALSE),by=scientific_name]


#check for unwanted name duplications (should be empty)
sampleAnnot[!duplicated(sampleAnnot[,c("abbreviation_sp","replicate","Patho-Nr"),with=FALSE])&duplicated(sampleAnnot[,c("abbreviation_sp","replicate"),with=FALSE])]
#also include tissue
sampleAnnot[!duplicated(sampleAnnot[,c("abbreviation_sp","replicate","Patho-Nr","abbreviation_ti","Tissue"),with=FALSE])&duplicated(sampleAnnot[,c("abbreviation_sp","replicate","abbreviation_ti"),with=FALSE])]

#make stats of listed species/samples
#unconverted species/samples sequenced
sub=sampleAnnot[!is.na(Flowcell_uc)&!is.na(lane_uc)]
message(paste0("Lanes used (unconv): ",nrow(sub[,.N,by=c("lane_uc","Flowcell_uc")]),"\nSpecies sequenced (unconv): ",length(unique(sub$scientific_name)),"\nSamples sequenced (unconv): ",length(unique(sub$Abbreviation))))
sort(sub$scientific_name)[duplicated(sort(sub$scientific_name))]

#prepped samples/ species
sub=sampleAnnot[!is.na(FlowCell)&!is.na(lane)]
message(paste0("Lanes used: ",nrow(sub[,.N,by=c("lane","FlowCell")]),"\nSpecies sequenced: ",length(unique(sub$scientific_name)),"\nSamples sequenced: ",length(unique(sub$Abbreviation))))
#some sub species obtained the same species abbreviation --> this is why there are more "sub species" than species with unconverted samples
#check that all species abbreviations have a unconverted sample (should be empty)
sub[!abbreviation_sp%in%sub[!is.na(Flowcell_uc)]$abbreviation_sp]

#prepped samples/ species
sub=sampleAnnot[!is.na(Adapter)&!is.na(Pool)]
message(paste0("Lanes used: ",nrow(sub[,.N,by=c("lane","FlowCell")]),"\nSpecies prepped: ",length(unique(sub$scientific_name)),"\nSamples prepped: ",length(unique(sub$Abbreviation))))

#listed samples/ species
sub=sampleAnnot[!is.na(`DNA (ng/µl)`)&(comment!="DNA failed"|is.na(comment))&`DNA (ng/µl)`>0]
message(paste0("Lanes used: ",nrow(sub[,.N,by=c("lane","FlowCell")]),"\nSpecies DNA available: ",length(unique(sub$scientific_name)),"\nSamples available: ",length(unique(sub$Abbreviation))))

#samples per species
samples_per_species=sub[,length(unique(Abbreviation)),by="scientific_name"]
table(samples_per_species$V1)

#only include complete sets
sampleAnnot_complete=sampleAnnot[as.logical(complete)&FlowCell!="failed"&lane!="failed"]

#include only selected species
#sampleAnnot_complete=sampleAnnot[abbreviation_sp=="EAB"&!is.na(FlowCell)&!is.na(lane)]
#sampleAnnot_complete=sampleAnnot[`Experiment ID`=="AK41_2"&!is.na(lane)]

#calculate number of used lanes and sequnced species (for keeping track)
message(paste0("Lanes used: ",nrow(sampleAnnot_complete[,.N,by=c("lane","FlowCell")]),"\nSpecies sequenced: ",length(unique(sampleAnnot_complete$scientific_name)),"\nSamples sequenced: ",length(unique(sampleAnnot_complete$Abbreviation))))

#get paths to bams
sampleAnnot_complete[,bam_path:=getBamPath(FlowCell,lane,ifelse(!is.na(renamed),gsub("__","_",renamed),Abbreviation),`Fortlaufende Nr`),by=1:nrow(sampleAnnot_complete)]
sampleAnnot_complete[,bam_path_uc:=getBamPath(Flowcell_uc,lane_uc,ifelse(!is.na(renamed_uc),renamed_uc,unconverted),"skip"),by=1:nrow(sampleAnnot_complete)]
sampleAnnot_complete[sampleAnnot_complete=="NA"]=NA

#only select subset
#sampleAnnot_complete=sampleAnnot_complete[abbreviation_sp=="TD"]
#sampleAnnot_complete=sampleAnnot_complete[abbreviation_sp%in%abbreviation_sp[!is.na(bam_path_uc)]]

#delete preexisting command file
overwrite=FALSE
commands_only=FALSE
run_wait=10
max_pending=5

dir.create("00_commands")

run_file_name="00_commands/RefFreeDMA_run_1.sh"
convCtr_file_name="00_commands/RefFreeDMA_convCtr_1.sh"
parse_file_name="00_commands/RefFreeDMA_parse_1.sh"

#make extra scripts for the new data (makes running easier)
run_file_name_new="00_commands/RefFreeDMA_run_new.sh"
convCtr_file_name_new="00_commands/RefFreeDMA_convCtr_new.sh"
parse_file_name_new="00_commands/RefFreeDMA_parse_new.sh"
cat("#!/bin/bash\n",file=run_file_name_new,append=FALSE)
cat("#!/bin/bash\n",file=convCtr_file_name_new,append=FALSE)
cat("#!/bin/bash\n",file=parse_file_name_new,append=FALSE)

if (overwrite==TRUE){
  cat("#!/bin/bash\n",file=run_file_name,append=FALSE)
  cat("#!/bin/bash\n",file=convCtr_file_name,append=FALSE)
  cat("#!/bin/bash\n",file=parse_file_name,append=FALSE)
  }
counter=0
for(run_species in unique(sampleAnnot_complete$abbreviation_sp)){
  message(run_species)
  counter=counter+1  
  
  #create directories
  dir.create(paste0(run_species,"/meta"),recursive=TRUE,showWarnings=FALSE)
  dir.create(paste0(run_species,"/unmapped_bam"),recursive=TRUE,showWarnings=FALSE)  
  
  
  #create sample annotation
  annot=sampleAnnot_complete[abbreviation_sp==run_species,c("Abbreviation","abbreviation_sp","Tissue","English"),with=FALSE]
  dupl=annot[duplicated(annot)]
  if (nrow(dupl)>0){
    entry=readline("INFO: Merging duplicate entries. Press ENTER to continue.")
    dupl
    annot=unique(annot)
  }
  sel=annot[,.N,by="Tissue"][order(N,decreasing=TRUE)]$Tissue[c(1:2)]
  annot[,comp:=ifelse(Tissue%in%sel,Tissue,NA),]
  setnames(annot,names(annot),c("Sample_Name","Organism","Tissue","Species","comp"))
  annot[,Species:=gsub(" ","",Species),]
  annot[,Select:=1,]
  
  #create config
  working_dir=paste0("working_dir=",wd,"/",run_species)
  species=paste0("species=",run_species)
  if (nrow(annot)<4){filtLim=paste0("filtLim=",nrow(annot)+1)}else{filtLim="filtLim=4"}
  config_add=paste(working_dir,species,filtLim,config,sep="\n")

  
  #create and store excecute command
  command_run=paste0("#",run_species,"\t",as.character(nrow(annot)),"\necho ",run_species,"\n",RefFreeDMA_dir,"/RefFreeDMA.sh ",wd,"/",run_species,"/meta/",run_species,"_RefFreeDMA.cfg > ",wd,"/00_RefFreeDMA_log/",run_species,".log&\nwhile [  `sacct --parsable --long|grep '|PENDING|'|wc -l` -gt ",max_pending," ]; do echo Too many pending!; sleep 5m; done\n")


  command_convCtr=paste0("#",run_species,"\t",as.character(nrow(annot)),"\necho ",run_species,"\n",convCtr_dir,"submit_conversionCtr.sh ",wd,"/",run_species," > ",wd,"/00_RefFreeDMA_log/",run_species,"_convCrt.log","&\nsleep 5m\n")
  command_parse=paste0("#",run_species,"\t",as.character(nrow(annot)),"\necho ",run_species,"\n",RefFreeDMA_dir,"/scripts/parse_stats.sh ",wd,"/",run_species,"/meta/",run_species,"_RefFreeDMA.cfg\n")
  
  
  #write files
  if (!file.exists(paste0(run_species,"/fastq/"))|overwrite==TRUE){
    if(commands_only!=TRUE){
    #link bams
    sub_annot=sampleAnnot_complete[abbreviation_sp==run_species][!duplicated(Abbreviation,fromLast=TRUE)]
    sub_annot[,makeLink(bam_path,run_species,Abbreviation),by=1:nrow(sub_annot)]
    if (any(!is.na(sub_annot$bam_path_uc))){
    sub_annot[!is.na(unconverted),makeLink(bam_path_uc,run_species,unconverted),by=1:nrow(sub_annot[!is.na(unconverted)])]
    uc_annot=sub_annot[!is.na(unconverted),list(Sample_Name=unconverted,Organism=abbreviation_sp,Tissue=Tissue,Species=gsub(" ","",English),comp=NA,Select=1),]
    annot=rbindlist(list(annot,uc_annot))  
    
    }
    
    #config
    cat(config_add,file=paste0(run_species,"/meta/",run_species,"_RefFreeDMA.cfg"))
    #sample annotation
    write.table(annot,paste0(run_species,"/meta/sampleAnnotation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    }
    #command
    #write to the existing scripts at the proper spot Attention: this makes them contain more than 30 commands
    cat(command_run,file=run_file_name,append=TRUE)
    cat(command_convCtr,file=convCtr_file_name,append=TRUE)
    cat(command_parse,file=parse_file_name,append=TRUE)
    #write to the "new" scripts
    cat(command_run,file=run_file_name_new,append=TRUE)
    cat(command_convCtr,file=convCtr_file_name_new,append=TRUE)
    cat(command_parse,file=parse_file_name_new,append=TRUE) 
  }else{print(paste0(run_species," Already setup!"))}

  
  if (counter%%30==0){
    run_file_name=paste0("00_commands/RefFreeDMA_run_",as.character(counter/30+1),".sh")
    convCtr_file_name=paste0("00_commands/RefFreeDMA_convCtr_",as.character(counter/30+1),".sh")
    parse_file_name=paste0("00_commands/RefFreeDMA_parse_",as.character(counter/30+1),".sh")
    
    if (overwrite==TRUE){
      cat("#!/bin/bash\n",file=run_file_name,append=FALSE)
      cat("#!/bin/bash\n",file=convCtr_file_name,append=FALSE)
      cat("#!/bin/bash\n",file=parse_file_name,append=FALSE)
    }
  }
}
