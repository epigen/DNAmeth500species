#!/bin/env Rscript

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(tableHTML)

wd=file.path(analysis_dir,"01_basicStats")
dir.create(wd,recursive=TRUE)
setwd(wd)


#functions
getStats=function(files){
  all_stats=data.table()
  for (file in files){
    stats=fread(file)
    
    if ("all_bases"%in%names(stats)){
    stats[,all_bases:=as.numeric(all_bases)/1000000,]
    setnames(stats,"all_bases","all_bases_Mio")}
    all_stats=rbindlist(list(all_stats,stats),use.names=TRUE,fill=TRUE)
  }
  return(all_stats)
}


#remove species to be excluded
sampleAnnot=sampleAnnot[!grepl("GPA_",Abbreviation)]


#fill all blanks with NA
sampleAnnot[sampleAnnot==""]=NA
sampleAnnot=sampleAnnot[FlowCell!="failed"&lane!="failed"]
sampleAnnot=sampleAnnot[!duplicated(Abbreviation,fromLast=TRUE)] #fromLast=TRUE to match how samples were selected for analysis in prepareRefFreeDMA
setnames(sampleAnnot,"Abbreviation","Sample_Name")
sampleAnnot[,`pre-BC_CT`:=as.numeric(gsub(",",".",`pre-BC_CT`))]

#reconstruct speciess abbreviation and replicate
sampleAnnot[,abbreviation_sp:=unlist(lapply(Sample_Name,function(x){unlist(strsplit(x,"_"))[1]}))]
sampleAnnot[,replicate:=unlist(lapply(Sample_Name,function(x){unlist(strsplit(x,"_"))[2]}))]
sampleAnnot[,Sample_Name_unif:=gsub("_uc","",Sample_Name),]

stats_files=list.files(processed_dir,pattern="^summary.txt",recursive=TRUE,full.names=TRUE)
overlap_files=system(paste0("ls ",processed_dir,"/*/*/diffMeth_cpg/*overlap*.tsv"),intern=TRUE)
ref_files=list.files(processed_dir,pattern="^ref_summary.txt",recursive=TRUE,full.names=TRUE)


stats=getStats(stats_files)
setnames(stats,"sample","Sample_Name")
overlap=getStats(overlap_files)
stats=merge(stats,overlap[,-grep("sample",names(overlap)),with=FALSE],by="Sample_Name",all.x=TRUE)
ref=getStats(ref_files)

#correct some sample names
stats[,Sample_Name:=sub("\\+","P",Sample_Name),]
stats[,Sample_Name:=sub("\\-","N",Sample_Name),]
stats[,Sample_Name:=sub("_macro","_Macro",Sample_Name),]

#add unified sample names (to also match the unconverted)
stats[,Sample_Name_unif:=gsub("_uc","",Sample_Name),]
#mark converted and unconverted samples
stats[,conversion_type:=ifelse(grepl("_uc",Sample_Name),"unconverted","converted"),]

#add some correlations (distance to diagonal)
stats[,res_mapRate_nonMotifReads:=(others+mapping_efficiency-100)/sqrt(2),]

#add overlap percentages
stats[,min_overlap_perc:=min_overlap/coveredCpGs*100,]
stats[,max_overlap_perc:=max_overlap/coveredCpGs*100,]

stats_annot=merge(stats,sampleAnnot[,-c("Sample_Name"),],by="Sample_Name_unif")
my_wt(stats_annot,"all_stats.tsv")


stats_long=melt(stats[,-c("max_cont_sp","blast_species1","blast_species2")],id.vars=c("Sample_Name","Sample_Name_unif","conversion_type","species"))
stats_long[,average:=mean(value,na.rm=TRUE),by="variable"]
stats_long[,sub_average:=ifelse(value<average,TRUE,FALSE)]
stats_long_annot=merge(stats_long,sampleAnnot[,-c("Sample_Name"),],by="Sample_Name_unif")

stats_long_annot[,value:=round(value,2),]
stats_long_annot[,average:=round(average,2),]

my_wt(stats_long_annot,"all_stats_long.tsv")

#only needed for sharing
#write.table(stats_long_annot,"/data/groups/lab_bock/public_html/jklughammer/compEpi_auto/patho_seq_stats.tsv",sep="\t",quote=FALSE,row.names=FALSE)
#cat(paste0("<a href='patho_seq_stats.tsv'>download table</a>\n",tableHTML(stats_long_annot,collapse="separate",theme="rshiny-blue",spacing = "4px")),file="/data/groups/lab_bock/public_html/jklughammer/compEpi_auto/patho_seq_stats.html")


ref_long=melt(ref[,-c("motifs"),with=FALSE],id.vars=c("species"))
ref_long[,average:=mean(value,na.rm=TRUE),by="variable"]
ref_long[,sub_average:=ifelse(value<average,TRUE,FALSE)]
my_wt(ref_long,"ref_stats_long.tsv")

