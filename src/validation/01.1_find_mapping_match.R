source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

suppressMessages(library(rtracklayer))
suppressMessages(library(lubridate))
suppressMessages(library(taxize))

## setting the working directory 
wd=file.path(analysis_dir,"validation", "01_crossMapping")
dir.create(wd,recursive=TRUE)
setwd(wd)

##### functions #######

#query ncbi for an annotation
get_db_anno=function(sci_name,db){
  print(sci_name)
  Sys.sleep(2)  
  tax=try(as.data.table(classification(gsub(" NBRC.*","",trimws(sci_name)),db, rows = 1)[[1]]),silent=TRUE) # if double UID -  we need an extra parameter, selected rows = 1 as the latest):
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
    df <- classification(ncbi_id,db='ncbi')[[1]]
    m <- max(which(df$rank !="clade"))
    new_id <- paste0("below", (NROW(df) -m), "_",  df[m,]$rank)
    return(new_id)
}
## extract lowest common for two ncbi ids
lowest_common_custom=function(ids,db="ncbi"){
  lc=try(lowest_common(ids,db=db))
  if (any(grepl("Error",lc))){
    lc=as.character(NA)
    return(lc)
  }
  lcr=lc$rank
  if (lcr=="below-no rank"){
    lcr=lc$name
    return(lcr)
  }else{
    if(lcr == "clade"){
        d <- classification(lc$id, 'ncbi')[[1]]
        m <- max(which(d$rank!="clade"))
        lcr = paste0("below", (NROW(d) - m), "_", d[m, "rank"])
    }
    return(lcr)
  } 
}


## find a crossmapping match
find_match <- function(idncbi, orderncbi, classncbi, full_ucsc_unique,tax_hierarchy){
    print(idncbi)
    if(idncbi %in% full_ucsc_unique$ncbi_id) return(list(idncbi, "species"))
    else{
        print(paste0("looking within ", orderncbi))
        if(NROW(full_ucsc_unique[ncbi_order == orderncbi]) == 0){
            print(paste0("no genomes avaliable for this order, looking within class ", classncbi))
            if(NROW(full_ucsc_unique[ncbi_class == classncbi]) == 0){
                print("no genomes avaliable for this class either")
               return(list(NA, NA)) 
            }else ids <- full_ucsc_unique[ncbi_class == classncbi]$ncbi_id
            
        }else{
            ids <- full_ucsc_unique[ncbi_order == orderncbi]$ncbi_id}
        print(ids)
        list_ranks <- lapply(setNames(ids, ids), function(x) lowest_common_custom(c(as.integer(x), idncbi)))
        list_ranks <- factor(list_ranks, levels = tax_hierarchy)
                             
        best_fits <- list_ranks[which(as.numeric(list_ranks) == max(as.numeric(list_ranks)))]
                             
        if(length(best_fits) == 1){
                return(list(names(best_fits), as.character(best_fits[[1]])))
        }else return(list(paste0(names(best_fits), collapse = ","), as.character(best_fits[[1]])))
            
        }
    }

## Analyzing the genomes avaliable from the UCSC

## extracting the full avaliable list
full_ucsc <- ucscGenomes()



## keeping only the latest from several assemblies
full_ucsc$date_my <- my(full_ucsc$date)
full_ucsc_unique <- full_ucsc %>% arrange(desc(date_my)) %>% group_by(species) %>% summarize(name = first(name), db = first(db), date = first(date))

head(full_ucsc_unique)

setDT(full_ucsc_unique)

full_ucsc_unique[, ncbi_name:=as.character(NA),]

## getting the ncbi ids for the ucsc data
simpleCache(cacheName="taxonomic_annotation_ucsc",instruction="full_ucsc_unique[is.na(ncbi_name),c('ncbi_name','ncbi_id','ncbi_order','ncbi_class','ncbi_group'):=get_db_anno(species,'ncbi'),by='species']",recreate=FALSE,assignToVariable="full_ucsc_unique")

head(full_ucsc_unique)

to_drop = c("C. japonica", "Chinese hamster ovary cell line","Yeast" ) ## flower, cell line,

#name fix so they are parsible by taxize
ncbi_names <- list("A. gambiae" = "Anopheles gambiae","A. mellifera" = "Apis mellifera", "Bushbaby" = "Galago", 
                    "C. brenneri" = "Caenorhabditis brenneri", "C. briggsae" = "Caenorhabditis briggsae", "C. elegans"= "Caenorhabditis elegans",
                  "C. intestinalis" = "Ciona intestinalis","C. remanei" =  "Caenorhabditis remanei", "Chimp" = "Pan troglodytes",
                  "D. ananassae" = "Drosophila ananassae","D. erecta" = "Drosophila erecta", "D. grimshawi" = "Drosophila grimshawi", 
                  "D. melanogaster" = "Drosophila melanogaster", "D. mojavensis" = "Drosophila mojavensis","D. persimilis" = "Drosophila persimilis",
                  "D. pseudoobscura" = "Drosophila pseudoobscura", "D. sechellia" = "Drosophila sechellia", "D. simulans" = "Drosophila simulans", "D. virilis" = "Drosophila virilis",
                  "D. yakuba" = "Drosophila yakuba", "Dolphin"="Tursiops truncatus","Garter snake" = "Thamnophis sirtalis", "Hedgehog" = "Erinaceus europaeus",
                  "Lamprey" = "Petromyzon marinus", "Lancelet" = "Cephalochordata", "Lizard" = "Anolis carolinensis", "Manatee" = "Trichechus manatus",
                  "Marmoset" = "Callithrix jacchus", "Megabat" = "Pteropus vampyrus", "Microbat" = "Myotis lucifugus", "Mouse lemur" = "Microcebus murinus",
                  "Opossum" = "Monodelphis domestica", "P. pacificus" = "Pristionchus pacificus", "Pika" = "Ochotona princeps", "Rhesus" = "Macaca mulatta",
                  "S. purpuratus" = "Strongylocentrotus purpuratus", "Shrew" = "Sorex araneus", "Sloth" = "Choloepus hoffmanni", "Squirrel" = "Spermophilus tridecemlineatus",
                  "Squirrel monkey" = "Saimiri boliviensis boliviensis", "Stickleback" = "Gasterosteus aculeatus","Wallaby" = "Macropus eugenii", "X. tropicalis" = "Xenopus tropicalis")

full_ucsc_unique <- full_ucsc_unique[!species %in% to_drop]

for (nameid in names(ncbi_names)){
    full_ucsc_unique[species == nameid, ncbi_name:=as.character(ncbi_names[[nameid]]),]
}

saveRDS(full_ucsc_unique, "full_ucsc_unique.RDS")

simpleCache(cacheName="taxonomic_annotation_ucsc2",instruction="full_ucsc_unique[is.na(ncbi_id),c('ncbi_name','ncbi_id','ncbi_order','ncbi_class','ncbi_group'):=get_db_anno(ncbi_name,'ncbi'),by='species']",recreate=FALSE,assignToVariable="full_ucsc_unique")

write.table(full_ucsc_unique, "avail_genomes_ucsc.tsv", sep = "\t", quote = F, row.names = F)

full_ucsc_unique <- fread("avail_genomes_ucsc.tsv")

mapped_species = data.table(species = stats_annot$species, ncbi_id = stats_annot$ncbi_id, 
                            ncbi_name = stats_annot$ncbi_name,class = stats_annot$ncbi_class, 
                            ncbi_order = stats_annot$ncbi_order,
                            ncbi_id_map = as.character(NA), common_level = as.character(NA))

mapped_species = unique(mapped_species)

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

head(mapped_species)

##pairwise scoring withon order/class to find the best matches
simpleCache(cacheName="taxonomic_annotation_mapping",
            instruction="mapped_species[is.na(ncbi_id_map),c('ncbi_id_map','common_level'):=find_match(ncbi_id,ncbi_order, class, full_ucsc_unique, tax_hierarchy ),by=c('species','ncbi_id')]",
            recreate=FALSE,assignToVariable="mapped_species")

##checking if there is really nothing more avaliable
mapped_species[is.na(ncbi_id_map),c('ncbi_id_map','common_level'):=find_match(ncbi_id,ncbi_order, class, full_ucsc_unique, tax_hierarchy ),by=c('species','ncbi_id')]

head(mapped_species)

mapped_species[, common_level := factor(common_level, levels = tax_hierarchy)]

mapped_species <- left_join(mapped_species, unique((stats_annot[, c("species", "color_class")])))

options(repr.plot.width = 10)
ggplot(mapped_species, aes(x = common_level, fill = color_class)) + geom_bar(position = "dodge") + 
    scale_fill_manual(values = class_colors) + rotate_labels()
ggsave("level_matched.pdf", height = 6, width = 12)

head(full_ucsc_unique)

full_ucsc_unique[,ncbi_id := as.character(ncbi_id),]

colnames(full_ucsc_unique) <- paste0("ucsc_", colnames(full_ucsc_unique))

## separating into individual pairs those that have more then one match
morethanonemap <- c()
N = NROW(mapped_species)
for (i in c(1:N )){
    if(length( strsplit(mapped_species[i,]$ncbi_id_map, ",")[[1]]) > 1){
        print(i)
        mapped_species <- rbind(mapped_species, data.table(species = mapped_species[i,]$species,
        ncbi_id = mapped_species[i,]$ncbi_id,
        ncbi_name = mapped_species[i,]$ncbi_name,
        class = mapped_species[i,]$class,
        ncbi_order = mapped_species[i,]$ncbi_order,
        ncbi_id_map = strsplit(mapped_species[i,]$ncbi_id_map, ",")[[1]],
        common_level = mapped_species[i,]$common_level,
        color_class = mapped_species[i,]$color_class)) 
    morethanonemap <- c(morethanonemap, i)
    }}

mapped_species <- mapped_species[-morethanonemap,]

mapped_species <- left_join(mapped_species, full_ucsc_unique, by=c("ncbi_id_map" = "ucsc_ncbi_id"))

write.table(mapped_species, "species_ucsc_matches.tsv", quote = F, row.names = F, sep = "\t")

write.table(unique(mapped_species[!is.na(ucsc_db), c("ucsc_db", "ucsc_species")]), "genomes_ucsc_todownload.tsv", sep = "\t", row.names = F, col.names = F, quote = F)
