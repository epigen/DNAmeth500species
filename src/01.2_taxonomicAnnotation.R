source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(taxize)


#load annotation data
anage=fread(file.path(meta_dir,"annotations_ext/anage_data.txt"))
patho_anno=fread(file.path(meta_dir,"Patholist_annotation.tsv"))

wd=file.path(analysis_dir,"01_basicStats")
dir.create(wd,recursive=TRUE)
setwd(wd)

#load stats data
stats=fread("all_stats.tsv")

#####functions#######

get_db_anno=function(sci_name,db){
  print(sci_name)
  Sys.sleep(2)  
  tax=try(as.data.table(classification(gsub(" NBRC.*","",trimws(sci_name)),db, rows = 1)[[1]]),silent=TRUE) # if double UID -  we need an extra parameter, selected rows = 1 as the latest):
  print(tax)
    #for those species still without annotation, try only the first part of the scientific name (=genus)
    
   if(is.na(tax[1, 1])){
      tax=try(as.data.table(classification(gsub(" .*","",trimws(sci_name)),db)[[1]]),silent=TRUE)}
    
  if(is.na(tax[1, 1])){return()}
  
  group=ifelse("Metazoa"%in%tax$name,ifelse("Vertebrata"%in%tax$name,"Vertebrata","Invertebrata"),as.character(NA))
  species=tail(tax,n=1)
  
  species_name=species$name
  species_id=species$id
  class=head(tax[rank%in%c("class")]$name,n=1)
  order=head(tax[rank%in%c("order")]$name,n=1)

  if (length(order)==0){
    sel=tax[which(rank=="superorder")+1,]
    order=sel$name
    print(paste0("No order. Using ", sel$rank))
  }

  return(list(species_name,species_id,order,class,group))
}

## addition to fix the non-ambiguity of clades
split_clades <- function(ncbi_id){
    print("specifying clade")
    df <- classification(ncbi_id,db='ncbi')[[1]]
    m <- max(which(df$rank !="clade"))
    new_id <- paste0("below", (NROW(df) -m), "_",  df[m,]$rank)
    return(new_id)
}


lowest_common_custom=function(ids,db="ncbi"){
  lc=try(lowest_common(ids,db=db),silent=TRUE)
   
  if (any(grepl("Error",lc))){
    lc=as.character(NA)
    return(lc)
  }
  lcr=lc$rank
  if (lcr=="below-no rank"){
    lcr=lc$name
    return(lcr)
  }else{
      if(lcr == "clade") lcr = split_clades(lc$id) 
    return(lcr)
  } 
}


##################
#annotate species
species_annot=data.table(scientific_name=unique(c(stats$scientific_name,stats$blast_species1,stats$blast_species2)),ncbi_name=as.character(NA))


#should not fail with error. If it does (probably random server connection loss or similar) just rerun on yet unannotated using "is.na(ncbi_name)"
## sometimes needs to be restarted (server interaction mistakes...)
simpleCache(cacheName="taxonomic_annotation",instruction="species_annot[is.na(ncbi_name),c('ncbi_name','ncbi_id','ncbi_order','ncbi_class','ncbi_group'):=get_db_anno(scientific_name,'ncbi'),by='scientific_name']",recreate=FALSE,assignToVariable="species_annot")

#check which are not yet annotated from ncbi (~1)
species_annot[is.na(ncbi_name)]

#add ncbi_name by hand (use higher ranking)
species_annot[scientific_name=="Colurodontis paxmani",ncbi_name:="Monacanthidae",]
species_annot[scientific_name=="Colurodontis paxmani",c("ncbi_name","ncbi_id","ncbi_order","ncbi_class","ncbi_group"):=get_db_anno(ncbi_name,"ncbi"),]

## Xenopus (Silurana) tropicalis -> Xenopus tropicalis
species_annot[scientific_name=="Xenopus (Silurana) tropicalis",ncbi_name:="Xenopus tropicalis",]
species_annot[scientific_name=="Xenopus (Silurana) tropicalis",c("ncbi_name","ncbi_id","ncbi_order","ncbi_class","ncbi_group"):=get_db_anno(ncbi_name,"ncbi"),]


#add Reptilia, Hyperoartia and Leptocardii as class (although not in ncbi taxonomy)
species_annot[ncbi_order%in%c("Crocodylia","Testudines","Squamata","Sphenodontia"),ncbi_class:="Reptilia",]
species_annot[ncbi_order%in%c("Petromyzontiformes"),ncbi_class:="Hyperoartia",]
species_annot[grep("Branchiostoma",ncbi_name),c("ncbi_order", "ncbi_class"):=list("Amphioxiformes","Leptocardii"),]

#add missing orders for included species
species_annot[is.na(ncbi_order)&!is.na(scientific_name)&!grepl("virus|[0-9]",scientific_name),ncbi_order:=get_db_anno(scientific_name,'ncbi')[3],by='scientific_name',]



#add still missing relevant orders by hand
manual_orders=c("Batoidea"="Batoidea","Amblyglyphidodon leucogaster"="Perciformes",
                "Dascyllus aruanus"="Perciformes","Pomacentrus coelestis"="Perciformes",
                "Embiotoca jacksoni"="Perciformes","Chromis punctipinnis"="Perciformes",
                "Centropomus robalito"="Perciformes","Saccoglossus kowalevskii"="Enteropneusta",
                "Oxycercichthys veliferus"="Perciformes","Centropyge heraldi"="Perciformes",
                "Morone saxatilis"="Perciformes","Cynoscion regalis"="Perciformes")

species_annot[,ncbi_order:=ifelse(scientific_name%in%names(manual_orders),manual_orders[scientific_name],ncbi_order),by=1:nrow(species_annot)]

#checks
species_annot[is.na(ncbi_id)]
species_annot[is.na(ncbi_order)]
species_annot[is.na(ncbi_class)] #14
## all other NAs are coming from the blast species - we can ignore them
#now save annotation


##now fixing the two names that were assigned wrong
species_annot[scientific_name == "Ostorhinchus monospilus"]$ncbi_name <- "Ostorhinchus moluccensis"
species_annot[scientific_name == "Ostorhinchus monospilus",c("ncbi_name","ncbi_id","ncbi_order","ncbi_class","ncbi_group"):=get_db_anno(ncbi_name,"ncbi"),]


species_annot[scientific_name == "Ostorhinchus gracilis"]$ncbi_name <- "Rhabdamia gracilis"
species_annot[scientific_name == "Ostorhinchus gracilis",c("ncbi_name","ncbi_id","ncbi_order","ncbi_class","ncbi_group"):=get_db_anno(ncbi_name,"ncbi"),]

my_wt(species_annot,"species_annot.tsv")

#Add missing orders to stats_annot if it already exists (fixin order retrospectively without recreating everything)
if(exists("stats_annot")){
  species_annot_c=copy(species_annot)
  setnames(species_annot_c,names(species_annot_c),paste0(names(species_annot_c),"_a"))

  stats_annot[,ncbi_order:=ifelse(is.na(ncbi_order),species_annot_c[scientific_name_a==scientific_name]$ncbi_order_a,ncbi_order),by=1:nrow(stats_annot)]

my_wt(stats_annot,"stats_annot.tsv")
save(stats_annot,file="stats_annot.RData")

}
###########################only continue here if running de novo#############################################
#now merge with stats
stats_annot=merge(stats,species_annot,by="scientific_name",all.x=TRUE)
stats_annot=merge(stats_annot,setnames(species_annot[,c("scientific_name","ncbi_id","ncbi_name"),],c("ncbi_id","ncbi_name"),c("ncbi_id_blastS1","ncbi_name_blastS1")),by.x="blast_species1",by.y="scientific_name",all.x=TRUE)
stats_annot=merge(stats_annot,setnames(species_annot[,c("scientific_name","ncbi_id","ncbi_name"),],c("ncbi_id","ncbi_name"),c("ncbi_id_blastS2","ncbi_name_blastS2")),by.x="blast_species2",by.y="scientific_name",all.x=TRUE)



### fixing mismatches in annotation:
 ## hand-replacing the annotations (where species were merged into one higher-level):
stats_annot[species == "WFA" ]$ncbi_id <- "12929"
stats_annot[species == "WFA" ]$ncbi_name <- "Amazona"

stats_annot[species == "GRE" ]$ncbi_id <- "30389"
stats_annot[species == "GRE" ]$ncbi_name <- "Ardea"

stats_annot[species == "SQM"]$ncbi_id <- "378850"
stats_annot[species == "SQM"]$ncbi_name <- "Saimiriinae"

##in case of missing samples:
stats_annot_backup <- stats_annot
#check species (comparison with blast species through lowest_common())
simpleCache(cacheName="taxonomic_lowest_common",instruction={stats_annot[,lowest_common_blastS1:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS1)),by=c("ncbi_id","ncbi_id_blastS1")]},recreate=FALSE,assignToVariable="stats_annot") 
simpleCache(cacheName="taxonomic_lowest_commonS2",instruction={stats_annot[,lowest_common_blastS2:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS2)),by=c("ncbi_id","ncbi_id_blastS2")]},recreate=FALSE,assignToVariable="stats_annot")

## for NA:
stats_annot <- rbind(stats_annot, stats_annot_backup[!Sample_Name %in% stats_annot$Sample_Name], fill = TRUE)

print(head(stats_annot[, c("lowest_common_blastS1", "lowest_common_blastS2")]))

print(NROW(stats_annot[is.na(lowest_common_blastS1)]))

stats_annot[is.na(lowest_common_blastS1),lowest_common_blastS1:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS1)),by=c("ncbi_id","ncbi_id_blastS1")]


print(NROW(stats_annot[is.na(lowest_common_blastS2)]))
      
stats_annot[is.na(lowest_common_blastS2),lowest_common_blastS2:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS2)),by=c("ncbi_id","ncbi_id_blastS2")]  


stopifnot(nrow(stats_annot[is.na(lowest_common_blastS1)])<=1) #1 NA row
stopifnot(nrow(stats_annot[is.na(lowest_common_blastS2)])<=4) #4 NA rows

## fixing by hand annotations that went wrong (ncbi mistakes) - depricated:
## GCF and BCL mapped to the same actinopteri nbci id

#create rank hierarchy
sort(unique(c(stats_annot$lowest_common_blastS1,stats_annot$lowest_common_blastS2)))

tax_hierarchy=c("cellular organisms","superkingdom","below1_superkingdom", "kingdom","below1_kingdom","below2_kingdom","below3_kingdom", "below4_kingdom","below5_kingdom",
                "phylum","below1_phylum","below2_phylum","subphylum","below1_subphylum", "below2_subphylum","below3_subphylum","below4_subphylum",
                "superclass","below1-superclass", "below2_superclass", "below3_superclass", "below4_superclass", "below5_superclass", "below6_superclass","below7_superclass",
                "class","below1_class","below2_class","below3_class","below4_class","below5_class","below7_class","subclass", "infraclass", 'below1_infraclass',"below-infraclass", 
                "cohort","below1_cohort","below2_cohort","below3_cohort","below4_cohort","below5_cohort","below6_cohort","below7_cohort","below-cohort", 
                "subcohort","below1_subcohort", "below2_subcohort","below3_subcohort","below4_subcohort","below5_subcohort","below6_subcohort","below7_subcohort",
                "superorder", "order", "below1_order","below2_order","below3_order", "below4_order", 
                "suborder","below1_suborder", "infraorder","parvorder",
                "superfamily", "family", "subfamily","tribe","genus","subgenus",
                "species", "subspecies")

stats_annot[,lowest_common_blastS1:=factor(lowest_common_blastS1,levels=tax_hierarchy),]
stats_annot[,lowest_common_blastS2:=factor(lowest_common_blastS2,levels=tax_hierarchy),]

stats_annot[,max_lowest_common:=max(c(as.numeric(lowest_common_blastS1),as.numeric(lowest_common_blastS2)),na.rm=TRUE),by=1:nrow(stats_annot)]

table(stats_annot$lowest_common_blastS1)
stats_annot[as.numeric(lowest_common_blastS1)<9&as.numeric(lowest_common_blastS2)<9]

#annotate with anAge
anage[,scientific_name:=paste0(Genus," ",Species),by=1:nrow(anage)]
stats_annot=merge(stats_annot,anage[,-c("HAGRID","Kingdom","Phylum","Class","Order","Family","Genus","Species", "Common name","Source","References"),],by="scientific_name",all.x=TRUE)

#annotate with patho info
#patho_anno
stats_annot=merge(stats_annot,unique(patho_anno),by.x=c("Fortlaufende Nr","Patho-Nr"),by.y=c("Nr","Patho-Nr"),all.x=TRUE)

print(stats_annot[species == "WFA" ]$ncbi_id )
my_wt(stats_annot,"stats_annot.tsv")
save(stats_annot,file="stats_annot.RData")
