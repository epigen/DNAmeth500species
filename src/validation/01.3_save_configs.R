source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

## setting the working directory 
wd=file.path(analysis_dir,"validation", "01_crossMapping")
dir.create(wd,recursive=TRUE)
setwd(wd)

mapped_species <- fread("species_ucsc_matches.tsv")

mapped_species <- mapped_species[!is.na(ncbi_id_map)]

ucscdb = mapped_species$ucsc_db[[1]]

max_pending = 200
cat("#cross-maping mode\n",file=file.path(processed_dir, "00_commands", "RefFreeDMA_concat.sh"))
write_mapped_config <- function(run_species, ucscdb){
    print(run_species)
#main config body  
config=paste0("
tool_path=/home/lv71484/droman/reffreedma/tools/
picard_path=$tool_path/picard-tools_1.118/
trim_galore_path=$tool_path/trim_galore_0.3.3/
cutadapt_path=$tool_path/cutadapt_1.8.3/
bowtie2_path=$tool_path/bowtie_2.2.4/bin/
bsmap_path=$tool_path/bsmap_2.90/
samtools_path=$tool_path/samtools_1.2/bin/
bedtools_path=$tool_path/bedtools_2.26.0/bin/
bwa_path=$tool_path/bwa_0.7.8/bin/
bwameth_path=$tool_path/bwa-meth-0.10/
decon_reference=/binfl/lv71484/droman/DNAmeth500species/resources/decon_reference/bacterial_extracted_add
cross_genome_fa=/binfl/lv71484/droman/DNAmeth500species/resources/reference_genomes/",ucscdb,"/", ucscdb,"_concatinated.fa
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
conversionCtr_path=/home/lv71484/droman/reffreedma/conversionCtr/
bisulfiteBlast_path=/home/lv71484/droman/reffreedma/bisulfiteBlast/
")
#create sample annotation
annot=stats_annot[abbreviation_sp==run_species,c("Sample_Name","abbreviation_sp","Tissue","English"),with=FALSE]
dupl=annot[duplicated(annot)]
if (nrow(dupl)>0){
   entry=readline(paste0("INFO: Merging duplicate entries for ",run_species," . Press ENTER to continue."))
   dupl
   annot=unique(annot)
}
sel=annot[,.N,by="Tissue"][order(N,decreasing=TRUE)]$Tissue[c(1:2)]
annot[,comp:=ifelse(Tissue%in%sel,Tissue,NA),]
setnames(annot,names(annot),c("Sample_Name","Organism","Tissue","Species","comp"))
annot[,Species:=gsub(" ","",Species),]
annot[,Select:=1,]    
#create config
working_dir=paste0("working_dir=",processed_dir,"/",run_species)
species=paste0("species=",run_species)
if (nrow(annot)<4){filtLim=paste0("filtLim=",nrow(annot)+1)}else{filtLim="filtLim=4"}
config_add=paste(working_dir,species,filtLim,config,sep="\n")

#create and store excecute command
command_run=paste0("#",run_species,"\t",as.character(nrow(annot)),"\necho ",run_species,"\n",RefFreeDMA_dir,
                   "/RefFreeDMA_VSC_expanded.sh ",processed_dir,"/",run_species,"/meta/",run_species,"_RefFreeDMA_",ucscdb, ".cfg -c > ",
                   processed_dir,"/00_RefFreeDMA_log/",run_species,"_", ucscdb, ".log&\nwhile [  `sacct --parsable --long|grep '|PENDING|'|wc -l` -gt ",
                   max_pending," ]; do echo Too many pending!; sleep 5m; done\n")

#write cmd
#cat("sleep 5m\n", file=file.path(processed_dir, "00_commands", "RefFreeDMA_concat.sh"),append=TRUE)
cat(command_run,file=file.path(processed_dir, "00_commands", "RefFreeDMA_concat.sh"),append=TRUE)
#write config
cat(config_add,file=file.path(processed_dir, paste0(run_species,"/meta/",run_species,"_RefFreeDMA_", ucscdb, ".cfg")))
}

for (i in seq_along(mapped_species$species)){
write_mapped_config( mapped_species$species[[i]],  mapped_species$ucsc_db[[i]])    
}
