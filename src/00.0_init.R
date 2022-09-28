#libraries

##technical
library(simpleCache, quietly = TRUE)
library(parallel)
suppressMessages(library(Rsamtools))

##working with data table
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

##vizualization
suppressMessages(library(ggplot2))
library(ggrepel)
library(ggseqlogo)
library(pheatmap)
suppressMessages(library(ComplexHeatmap)) ##should be installed directly from gitHub
#suppressMessages(library(circlize))



#plotting settings
theme_set(theme_bw())
theme_update(axis.text = element_text(color="black"),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             text = element_text(size=12))

pdf.options(useDingbats=FALSE)

#plotting helper functions

rotate_labels <- function(angle = 60,vjust = 1,hjust = 1)
        {return(theme(axis.text.x = element_text(angle = angle, vjust = vjust,hjust = hjust)))}

#use like this: stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)
give.n <- function(x,y = NULL){
  return(data.frame(y = ifelse(is.null(y),0,y), label = paste0("N=",length(x)))) 
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

##abbreviations for plotting
class_short <- c("Invertebrata" = "Inv.",
                 "Jawless_vertebrate" = "Jl.vb.",
                 "Chondrichthyes" = "Ch.",
                 "Actinopteri" = "Act.", "Amphibia" = "Amp.", 
                 "Reptilia" = "Rep.", "Aves" = "Av.", 
                 "Marsupialia" = "Mar.", "Mammalia" = "Mam.")

#colors
class_colors <- c("Invertebrata" = "#cfcfcf",
                  "Jawless_vertebrate" = "#636363",
                  "Chondrichthyes" = "#F68383", 
                  "Actinopteri" = "#EA4142",
                  "Amphibia" = "#5AB349",
                  "Reptilia" = "#8761AC",
                  "Aves"="#FE9222",
                  "Marsupialia"="#8EBFDA",
                  "Mammalia"="#4892C2")

tissue_colors <- c("H" = "#005fb1", "L" = "#e6ab02", "S" = "#a9a900",
                   "M" = "#ae6600", "B" = "#D1B0FF", "G" = "#59ae00",
                   "F" = "#00a783", "X" = "#525252")

#Useful vectors

#classes
classes <- c("Invertebrata","Jawless_vertebrate",  "Chondrichthyes", "Actinopteri", "Amphibia",
           "Reptilia", "Aves", "Marsupialia", "Mammalia")

#core tissues
core_inv <- c("Muscle","Tube_feet","Tentacle","Arm","Gills","Gonad","Pharynx")
core_vert <- c("Liver","Heart","Brain","Spleen","Muscle","Gills","Fin")
core <- c("Liver","Heart","Lung", "Brain","Spleen","Muscle","Gills","Fin", "Kidney")

#paths (have to be updated)
RefFreeDMA_dir="/home/lv71484/droman/reffreedma/RefFreeDMA" 
convCtr_dir = "/home/lv71484/droman/reffreedma/conversionCtr"
bisulfiteBlast_dir = "/home/lv71484/droman/reffreedma/bisulfiteBlast"
BSF_dir = "/binfl/lv71484/droman/DNAmeth500species/raw_data"

data_dir = "/binfl/lv71484/droman/DNAmeth500species/"
source_dir = file.path(Sys.getenv("CODEBASE"),"DNAmeth500species")
rrbs_dir = "toSelf_filtered_0.08mm_final_concat"


#paths (inherited)
processed_dir = file.path(data_dir,"results_pipeline")
analysis_dir = file.path(data_dir,"results_analysis")
cache_dir = file.path(data_dir,"RCache")
dir.create(cache_dir,showWarnings = FALSE)
options(RCACHE.DIR = cache_dir)

meta_dir = file.path(source_dir,"meta")
scripts_dir = file.path(source_dir,"scripts")


#annotations
stats_annot_file <- file.path(analysis_dir,"01_basicStats/stats_annot.RData")

if (file.exists(stats_annot_file)){
  load(stats_annot_file)
  stats_annot[,color_class:=ifelse(ncbi_class %in% c("Mammalia","Amphibia","Reptilia","Actinopteri","Chondrichthyes","Aves"),
                                   ncbi_class,"Invertebrata"),]
  stats_annot[ncbi_order%in%c("Diprotodontia","Dasyuromorphia"),color_class:="Marsupialia",]
    stats_annot[species == "JL",color_class:="Jawless_vertebrate", ]
  stats_annot[,color_class:=factor(color_class,levels=classes),]
 
    
    stats_annot[, group := class_short[color_class],]
   stats_annot[, group:=factor(group, levels = class_short)]
    stats_annot[species == "JL", group := "Inv.",] #so that its visualised with inv.
    ### assigning to the JL color_class/group NA (as it is a jawless vertebrate and doesn't fit in any group)
    
    #
   
    ## fixing the 
  }else{
  print("stats_annot.RData doesn't exist. Run script 01.2 to create. Only loading basic sample annotation.")
}

sampleAnnot=fread(file.path(meta_dir,"Patholist_selection.tsv"))

##simple helpful dataframe species-class, ordered phylogenetically 

phylo_order_path <- file.path(Sys.getenv("CODEBASE"), "DNAmeth500species/meta/species_list_ordered_2021_2.txt")

if(file.exists(phylo_order_path) & file.exists(stats_annot_file)){
  
  sp_df <- read.table(phylo_order_path)
  colnames(sp_df) <- c("species")
  sp_df <- unique(left_join(sp_df, stats_annot[, c("species", "color_class", "group")]))
  row.names(sp_df) <- sp_df$species
  setDT(sp_df)

  
}else{
  print("species_list_ordered.txt doesn´t exist. Run 01.7 to create") ##
}

#helper functions
my_wt=function(tab,name){
  write.table(tab,name,sep="\t",quote=FALSE,row.names=FALSE)
}



