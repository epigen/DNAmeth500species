source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))


#different species as analyzed
length(unique(stats_annot$abbreviation_sp))

#different species as annotated by scientific name
length(unique(stats_annot$scientific_name))

#number of samples
stats_annot[,length(unique(Sample_Name)), by=conversion_type]

#average samples per species
mean(stats_annot[conversion_type=="converted",.N,by=abbreviation_sp]$N)

#average tissues per species
mean(stats_annot[conversion_type=="converted",length(unique(Tissue)),by=abbreviation_sp]$V1)

#average individuals per species
mean(stats_annot[conversion_type=="converted",length(unique(replicate)),by=abbreviation_sp]$V1)
