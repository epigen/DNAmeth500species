source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(phylobase)
library(tidytree)
library(ape)
wd=file.path(analysis_dir,"01_basicStats/01.6_ITOL")
dir.create(wd,recursive=TRUE)
setwd(wd)


#annotation settings
create_annot_bar=function(label,colors,fields,scale,width,legend_col,internal=0,data){
  settings=paste0("DATASET_MULTIBAR
                  SEPARATOR COMMA",
                  "\nDATASET_LABEL,",label,
                  "\nCOLOR,",legend_col,
                  "\nFIELD_COLORS,",colors,
                  "\nFIELD_LABELS,",fields,
                  "\nDATASET_SCALE,",scale,
                  "\nWIDTH,",width,
                  "\nSHOW_INTERNAL,",internal,
                  "\nDATA\n",paste(data,collapse="\n"))
  return(settings)
}


create_annot_popup=function(data){
  settings=paste0("POPUP_INFO
SEPARATOR COMMA
DATA\n",paste(data,collapse="\n"))
  return(settings)
}

create_annot_gradient=function(data,label,legend_col){
  settings=paste0("DATASET_GRADIENT
SEPARATOR COMMA
STRIP_WIDTH,70
BORDER_WIDTH,3
SHOW_INTERNAL,1
DATASET_LABEL,",label,
"\nCOLOR,",legend_col,
"\nDATA\n",paste(data,collapse="\n"))
  return(settings)
}

tree_colors="
TREE_COLORS
SEPARATOR COMMA
DATA\n
"
#create stats per species (summarize replicate)
make_stats_per_species=function(stats){

  stats_per_species=stats[!is.na(ncbi_id),list(abbreviation=abbreviation_sp[1],scientific_name=scientific_name[1],
                                               common_name=English[1],german_name=Tierart[1],class=ncbi_class[1],
                                               As=mean(perc_As),Ts=mean(perc_Ts),Gs=mean(perc_Gs),Cs=mean(perc_Cs),Ns=mean(perc_Ns),
                                               CGG=100*sum(total_reads_untrimmed*CGG)/sum(total_reads_untrimmed*(100-others)),
                                               TGG=100*sum(total_reads_untrimmed*TGG)/sum(total_reads_untrimmed*(100-others)),
                                               CGA=100*sum(total_reads_untrimmed*CGA)/sum(total_reads_untrimmed*(100-others)),
                                               TGA=100*sum(total_reads_untrimmed*TGA)/sum(total_reads_untrimmed*(100-others)),
                                               other_reads=100*sum(total_reads_untrimmed*others)/sum(total_reads_untrimmed*100),rand_frag=mean(others),
                                               avg_meth=mean(avg_meth),nonCpG_meth=mean(100-conversionRate),CpG_meth=mean(CpG_meth),samples=.N,CpGs=mean(coveredCpGs),
                                               tissues=length(unique(Tissue)),contamination=mean(cont_rat),mapping_efficiency=mean(mapping_efficiency),
                                               conversionRate=mean(1-k1_unmeth),fragments_ref=fragments_ref[1],
                                               lowest_rank=min(max_lowest_common,na.rm=TRUE)),by=c("ncbi_id","ncbi_id_int","ncbi_name","color_class")]
  return(stats_per_species)
}


#output species to create tree 
#Use https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
#include unranked taxa
#expand all (might need manual expanding to make sure all are expanded)
#cat(na.omit(unique(stats_annot$ncbi_id)),file="species_list_ncbiID.txt",sep="\n") #only do once
<<<<<<< Updated upstream
#save tree as species_tree_unranked.phy in the meta directory

#read tree to set NCBI ids and fix internal nodes
#tree=readLines(paste0(source_dir, "/meta/species_tree_unranked.phy"))

=======


#read tree to set NCBI ids and fix internal nodes
#tree=readLines(paste0("species_tree_unranked.phy"))
tree_id = "2021_2" ##adjust 
tree=readLines(paste0(meta_dir, "/species_tree_unranked_", tree_id, ".phy"))
>>>>>>> Stashed changes
tree=paste0(gsub("'","",tree),collapse="")

writeLines(tree,con=paste0("species_tree_unranked_names_", tree_id, ".phy"))

replace=unique(stats_annot[,c("ncbi_id","ncbi_name"),])

for(nameid in replace$ncbi_name){
    if(!grepl(paste0(nameid,":"), tree)){
        print(nameid)
    }
}


for (i in 1:nrow(replace)){
  tree=with(replace[i,],gsub(paste0(")",ncbi_name,":"),paste0(")INT",ncbi_id,":"),tree,))
  tree=with(replace[i,],gsub(paste0(ncbi_name,":"),paste0(ncbi_id,":"),tree,))
}
writeLines(tree,con=paste0("species_tree_unranked_id_", tree_id, ".phy"))


##checking which species are missing and/or which labels are not present
tree_df <- read.tree(paste0("species_tree_unranked_id_", tree_id, ".phy"))

labels<-tree_df$tip.label
nodes<-tree_df$node.label

df <- as_tibble(tree_df)

stopifnot(length(setdiff(labels, stats_annot$ncbi_id))==0)

int_mapped <- paste0("INT", setdiff(stats_annot$ncbi_id, labels))
stopifnot(length(setdiff(int_mapped, nodes))==0)


#needed if tree is built on ncbi ids (instead of names)
stats_annot[,ncbi_id_int:=ifelse(any(grepl(paste0("INT",ncbi_id[1],":"),tree)),paste0("INT",ncbi_id[1]),as.character(ncbi_id[1])),by=ncbi_id]

                            

                            
                            
#for unconverted
stats_per_species_unconv=make_stats_per_species(stats_annot[conversion_type=="unconverted"])
#for converted
stats_per_species_conv=make_stats_per_species(stats_annot[conversion_type=="converted"])



#general species annotation
legend_col="#888888"
#group color
stats_per_species_conv[,col:=class_colors[as.character(color_class)],by=1:nrow(stats_per_species_conv)]
stats_per_species_conv[,annot_col:=paste(ncbi_name,"range",paste0(col,"89"),color_class,sep=","),]
cat(tree_colors,file=paste0("colors_new.txt"),sep=",")
cat(stats_per_species_conv$annot_col,file=paste0("colors_new.txt"),sep="\n",append=TRUE)

#species info popup
cat(create_annot_popup(data=stats_per_species_conv[,paste(ncbi_name,scientific_name,paste0("<p>Common name: ",common_name,"</p>","<p>Abbreviation: ",abbreviation,"</p>","<p>NCBI ID: ",ncbi_id,"</p>","<p style='color:blue'>More info at <a target='_blank' href=",paste0("'https://www.google.at/search?q=",gsub(" ","+",scientific_name),"'> google."),"</a></p>"),sep=","),]),file=paste0("popup.txt"))

#Number of sample
cat(create_annot_bar(label=paste0("Number of samples"),colors="#75b1c7",fields="samples",scale="0,10,20,30,40",width=200*max(stats_per_species_conv$samples)/40,legend_col=legend_col,internal=1,data=stats_per_species_conv[,paste(ncbi_name,samples,sep=","),]),file=paste0("annot_samples.txt"))

#Number of tissues
cat(create_annot_bar(label=paste0("Number of tissues"),colors="#75b1c7",fields="tissues",scale="0,1,2,3,4,5,6,7,8,9,10",width=200*max(stats_per_species_conv$tissues)/10,legend_col=legend_col,internal=1,data=stats_per_species_conv[,paste(ncbi_name,tissues,sep=","),]),file=paste0("annot_tissues.txt"))

#Number of reference fragments
cat(create_annot_bar(label=paste0("Reference fragments"),colors="#75b1c7",fields="fragments_ref",scale="0,1*10^6,2*10^6,3*10^6,4*10^6,5*10^6,6*10^6,7*10^6,8*10^6,9*10^6,10^7",legend_col=legend_col,width=200*max(stats_per_species_conv$fragments_ref)/10^7,internal=1,data=stats_per_species_conv[,paste(ncbi_name,fragments_ref,sep=","),]),file=paste0("annot_fragments_ref.txt"))


#stats annotation (separate datasets for converted and unconverted)
set_list=list(c("stats_per_species_conv","converted","#ffffff"),c("stats_per_species_unconv","unconverted","#000000"))

for (sel_set in set_list) {

  stats_per_species=get(sel_set[1])
  suffix=sel_set[2]
  legend_col=sel_set[3]

 cat(create_annot_bar(label=paste0("Number of covered CpGs (1k) ",suffix),colors="#ff9e15",fields="CpG_count",scale="0,1000,2000,3000,4000,5000,6000,7000,8000",width=200,legend_col=legend_col,data=stats_per_species[,paste(ncbi_name,round(CpGs/1000,0),sep=","),]),file=paste0("annot_CpGs_count_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Restriction Fragemnt Distribution ",suffix),colors="#2b8cbe,#a6bddb,#31a354,#a1d99b",fields="CGG,TGG,CGA,TGA",scale="0,10,20,30,40,50,60,70,80,90,100",width=200,legend_col=legend_col,data=stats_per_species[,paste(ncbi_name,CGG,TGG,CGA,TGA,sep=","),]),file=paste0("annot_resFrag_",suffix,".txt"))
  
  cat(create_annot_bar(label=paste0("Base distribution ",suffix),colors="#66c2a5,#fc8d62,#8da0cb,#e78ac3,#888888",fields="As,Ts,Cs,Gs,Ns",scale="0,10,20,30,40,50,60,70,80,90,100",width=200,legend_col=legend_col,data=stats_per_species[,paste(ncbi_name,As,Ts,Cs,Gs,Ns,sep=","),]),file=paste0("annot_baseDistribution_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Random fragmentation ",suffix),colors="#10ea77",fields="rand_frag",scale="0,10,20,30,40,50,60,70,80,90,100",width=200*max(stats_per_species$rand_frag)/100,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,rand_frag,sep=","),]),file=paste0("annot_randFrag_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("CpG methylation ",suffix),colors="#ff9e15",fields="CpG_meth",scale="0,10,20,30,40,50,60,70,80,90,100",width=200*max(stats_per_species$CpG_meth)/100,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,CpG_meth,sep=","),]),file=paste0("annot_CpG_meth_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("NonCpG methylation ",suffix),colors="#ff9e15",fields="nonCpG_meth",scale="0,10,20,30,40,50,60,70,80,90,100",width=200*max(stats_per_species$nonCpG_meth)/100,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,nonCpG_meth,sep=","),]),file=paste0("annot_nonCpG_meth_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Average methylation ",suffix),colors="#ff9e15",fields="avg_meth",scale="0,10,20,30,40,50,60,70,80,90,100",width=200*max(stats_per_species$avg_meth)/100,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,avg_meth,sep=","),]),file=paste0("annot_avg_meth_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Contamination ratio ",suffix),colors="#10ea77",fields="contamination",scale="0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1",width=200*max(stats_per_species$contamination)/1,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,contamination,sep=","),]),file=paste0("annot_contamination_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Conversion rate ",suffix),colors="#10ea77",fields="conversionRate",scale="0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1",width=200*max(stats_per_species$conversionRate)/1,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,conversionRate,sep=","),]),file=paste0("annot_conversionRate_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Mapping effieciency ",suffix),colors="#10ea77",fields="mapping_efficiency",scale="0,10,20,30,40,50,60,70,80,90,100",width=200*max(stats_per_species$mapping_efficiency)/100,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,mapping_efficiency,sep=","),]),file=paste0("annot_mapping_efficiency_",suffix,".txt"))

  cat(create_annot_bar(label=paste0("Lowest common rank ",suffix),colors="#10ea77",fields="lowest_rank",scale="0,3,6,9,12,15,18,21,24,27,30,33",width=200*max(stats_per_species$lowest_rank,na.rm=TRUE)/33,legend_col=legend_col,internal=1,data=stats_per_species[,paste(ncbi_name,lowest_rank,sep=","),]),file=paste0("annot_lowest_rank_",suffix,".txt"))

  cat(create_annot_gradient(label=paste0("Lowest common rank ",suffix),legend_col=legend_col,data=stats_per_species[,paste(ncbi_name,lowest_rank,sep=","),]),file=paste0("annot_lowest_rank_gradient_",suffix,".txt"))
}
<<<<<<< Updated upstream
=======


>>>>>>> Stashed changes
