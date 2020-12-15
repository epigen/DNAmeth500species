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
  tax=try(as.data.table(classification(gsub(" NBRC.*","",trimws(sci_name)),db)[[1]]),silent=TRUE)
  #for those species still without annotation, try only the first part of the scientific name (=genus)
   if(is.null(ncol(tax))){
      tax=try(as.data.table(classification(gsub(" .*","",trimws(sci_name)),db)[[1]]),silent=TRUE)}
  if(is.null(ncol(tax))){return()}
  
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
    return(lcr)
  } 
}


##################
#annotate species
species_annot=data.table(scientific_name=unique(c(stats$scientific_name,stats$blast_species1,stats$blast_species2)),ncbi_name=as.character(NA))

#should not fail with error. If it does (probably random server connection loss or similar) just rerun on yet unannotated using "is.na(ncbi_name)"
simpleCache(cacheName="taxonomic_annotation",instruction="species_annot[is.na(ncbi_name),c('ncbi_name','ncbi_id','ncbi_order','ncbi_class','ncbi_group'):=get_db_anno(scientific_name,'ncbi'),by='scientific_name']",recreate=FALSE,assignToVariable="species_annot")

#check which are not yet annotated from ncbi (~1)
species_annot[is.na(ncbi_name)]

#add ncbi_name by hand (use higher ranking)
species_annot[scientific_name=="Colurodontis paxmani",ncbi_name:="Monacanthidae",]
species_annot[scientific_name=="Colurodontis paxmani",c("ncbi_name","ncbi_id","ncbi_order","ncbi_class","ncbi_group"):=get_db_anno(ncbi_name,"ncbi"),]


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

#now save annotation
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

#check species (comparison with blast species through lowest_common())
stats_annot[,lowest_common_blastS1:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS1)),by=c("ncbi_id","ncbi_id_blastS1")]
stats_annot[,lowest_common_blastS2:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS2)),by=c("ncbi_id","ncbi_id_blastS2")]

#IMPORTANT:repeat to fix NAs due to server timeout
stats_annot[is.na(lowest_common_blastS1),lowest_common_blastS1:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS1)),by=c("ncbi_id","ncbi_id_blastS1")]
stats_annot[is.na(lowest_common_blastS2),lowest_common_blastS2:=lowest_common_custom(c(ncbi_id,ncbi_id_blastS2)),by=c("ncbi_id","ncbi_id_blastS2")]

stopifnot(nrow(stats_annot[is.na(lowest_common_blastS1)])<=1) #1 NA row
stopifnot(nrow(stats_annot[is.na(lowest_common_blastS2)])<=4) #4 NA rows

#create rank hierarchy
sort(unique(c(stats_annot$lowest_common_blastS1,stats_annot$lowest_common_blastS2)))

tax_hierarchy=c("cellular organisms","superkingdom","kingdom","below-kingdom","phylum","below-phylum","below-subphylum","superclass","below-superclass","class","below-class","subclass","infraclass","below-infraclass","superorder","below-superorder","order","below-order","suborder","below-suborder","infraorder","parvorder","superfamily","family","subfamily","tribe","genus","subgenus","species","subspecies")
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

my_wt(stats_annot,"stats_annot.tsv")
save(stats_annot,file="stats_annot.RData")
