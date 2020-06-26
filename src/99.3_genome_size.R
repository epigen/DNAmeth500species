source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(stringr)

wd=file.path(analysis_dir,"01_basicStats/01.6_stats_detail")
dir.create(wd)
setwd(wd)

#genome sizes downloaded from: https://www.ncbi.nlm.nih.gov/genome/browse/#!/eukaryotes/
genome_sizes=fread(paste0(meta_dir,"/annotations_ext/eukaryotes.csv"))

#simplify organism names for better matches
genome_sizes[,`#Organism Name`:=word(`#Organism Name`,1,2," "),by=1:nrow(genome_sizes)]

#some organisms have multiple records --> only keep the largest genome size
genome_sizes=genome_sizes[order(`Size(Mb)`,decreasing=TRUE)][!duplicated(`#Organism Name`)]


stats_annot_gs=merge(stats_annot,genome_sizes,by.x="ncbi_name",by.y="#Organism Name",all.x=TRUE)

length(unique(stats_annot_gs[!is.na(`Size(Mb)`)]$ncbi_name))

stats_annot_gs_red=unique(stats_annot_gs[!is.na(`Size(Mb)`),c("ncbi_name","fragments_ref","abbreviation_sp","color_class","Size(Mb)","GC%","Level"),])


ggplot(stats_annot_gs_red,aes(x=log10(`Size(Mb)`),y=log10(fragments_ref),col=color_class))+geom_point()+geom_text(aes(label=abbreviation_sp))+scale_color_manual(values = class_colors)
