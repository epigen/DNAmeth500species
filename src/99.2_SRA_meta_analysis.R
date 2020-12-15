source(file.path(Sys.getenv("CODEBASE"),"compEpi/src/00.0_init.R"))

dat=fread(paste0(meta_dir,"/sra_result.csv"))
dat[,nSpecies:=length(unique(`Organism Name`)),by="Study Accession"]
dat=dat[order(nSpecies,decreasing=TRUE)]

unique(dat[,c("Study Title","Study Accession","nSpecies")])[1:10]

#Study Title Study Accession nSpecies
#1: The evolution of CHROMOMETHYLTRANSFERASES and gene body DNA methylation in plants       SRP075487       34
#2:                      Sociality and DNA methylation are not evolutionary dependent       SRP076799       30  Insects
#3:                Widespread natural variation of DNA methylation within angiosperms       SRP072226       23
#4:                   Genome-wide evolutionary analysis of eukaryotic DNA methylation       SRP002417       17  Across all eucaryotes  10 Animals, 2 Vertebrates
#5:                                    Evolution of DNA methylation across arthropods       SRP230024       11
#6:                                    The DNA Methylation Landscape of Giant Viruses       SRP252852       10
#7:    Single-molecule detection of N6-methyladenine in microbial reference materials       SRP151172       10
#8:                         Genome size and DNA methylation covary across Land Plants       SRP059572        9
#9:                                               DNA methylation in cotton evolution       SRP071640        8
#10:                                            DNA methylation in mammalian placentas       SRP049936        8

# paste Study Accession here to see species
unique(dat[`Study Accession`=="SRP049936"]$`Organism Name`)
