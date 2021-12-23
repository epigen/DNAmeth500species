## helper that extracts the required order of the the tips and 
##internal nodes, that are mapped to the species
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

setwd(analysis_dir)

library(phylobase)
library(tidytree)
library(ape)
library(ggtree) ##installation note - downgrade rvcheck to 0.1.8, if R < 4.1.1!!!
#treefile<-
tree <- read.tree("01_basicStats/01.6_ITOL/species_tree_unranked_id_2021.phy")

labels<-tree$tip.label
nodes<-tree$node.label

## retriving the internaly mapped to adiitioinal tips
ncbi_ids_tree<-unique(stats_annot$ncbi_id)

## if the ncbi_id is not a tip, that means that it is mapped internaly (and have the id INTXX)
internal_mapped<-paste0("INT", setdiff(ncbi_ids_tree, labels))

df <- as_tibble(tree)
df_internal <- df[df$label %in% internal_mapped, ]

##getting the true ncbi_id
df_internal["id"] <- sapply(df_internal$label, function(x) strsplit(x, "T")[[1]][2])
##tis new tips will have the currently assigned nodes as parents:
df_internal$parent <- df_internal$node
##new tip ids - start right after the last one, which is the length of tip.labels
df_internal$node <- as.integer(row.names(df_internal))+length(labels)
##than we shift all the node ids by the number of added tips
df_internal$parent <- df_internal$parent+NROW(df_internal)
df_internal$label <- df_internal$id
df_internal <- df_internal[, setdiff(colnames(df_internal), "id")]


##now we have to update the df of the initial tree: 
##for that we shift all the parents and all the internal node ids (ones, that are bigger that the length of )
df$parent <- df$parent+NROW(df_internal)

df[df$node>length(labels), ]$node <- df[df$node>length(labels), ]$node + NROW(df_internal)

df_new <- rbind(df, df_internal)




##now we relabel the nodes to species abbreviation
df_new$label_sp <- df_new$label
new_L <- length(labels)+NROW(df_internal) ## number of nodes

df_new[df_new$node<new_L+1,]$label_sp <- unlist(sapply(df_new[df_new$node<new_L+1,]$label,
                        function(x) as.character(stats_annot[stats_annot$ncbi_id==x, "species"][[1]][1])))


##checking if there are non -unique tips:
dupl_sp<-df_new[df_new$node<new_L,]$label_sp[duplicated(df_new[df_new$node<new_L,]$label_sp)]

df_new[df_new$label_sp %in% dupl_sp,]


######depricated, tree is based on a clean annotation now########
##saving those, we want to drop:
#217, 231, 576, 311, 490
#delete_tip<-function(df, tip_id){
#  df_new<-df[df$node!=tip_id, ]
#  df_new$node_updt<-df_new$node
#  df_new[df_new$node>tip_id, ]$node_updt<-df_new[df_new$node>tip_id, ]$node - 1
#  df_new$parent_updt<-df_new$parent
#  df_new[df_new$parent>tip_id, ]$parent_updt<-df_new[df_new$parent>tip_id, ]$parent - 1
#  df_new$parent<-df_new$parent_updt
#  df_new$node<-df_new$node_updt
#  return(df_new[, setdiff(colnames(df_new), c("node_updt", "parent_updt"))])
#}
#we need to start with the biggest number, so that nothing chages
#deleting node 
#df_new<-delete_tip(df_new, 576)
#df_new<-delete_tip(df_new, 490)
#df_new<-delete_tip(df_new, 311)
#df_new<-delete_tip(df_new, 231)
#df_new<-delete_tip(df_new, 217)


#now we have to check if we actually have ALL the species
missed <- setdiff(unique(stats_annot$species), df_new$label_sp)
missed
                                                       
amb_ids <- unique(stats_annot[stats_annot$species %in% missed, ]$ncbi_id)
amb_ids

unique(stats_annot[stats_annot$ncbi_id %in% amb_ids, ]$species)

######depricated, tree is based on a clean annotation now########
#add_tip<-function(df, parent, label, sp_label, N_of_tips){
  ##all that are bigger than the new node id should be shifted +1 
#  df[df$node>N_of_tips,]$node<-df[df$node>N_of_tips,]$node+1
  ##same for the parents (but we are shifting ALL, as they are all after the tips)
 # df$parent<-df$parent+1
  #new row to add
  #new_row<-c(parent+1, N_of_tips+1, 4, label, sp_label)
  #df_new<-rbind(df, new_row)
  #df_new$node<-as.integer(df_new$node)
  #df_new$parent<-as.integer(df_new$parent)
  #df_new$branch.length<-as.integer(df_new$branch.length)
  #return(df_new)
#}


## second is 
#split_tip<-function(df, tip_label, N_of_tips, sp_id){
  #<-df[df$label==tip_label, ]$parent
  #creating a new internal nofr
  #df_updt<-rbind(df, c(p, NROW(df)+1, 4, tip_label, paste0("INT",tip_label)))
  #df_updt[c(df$label==tip_label, F),"label"]<-paste0(tip_label, "_1")
  #df_updt[c(df$label==tip_label, F),"parent"]<-NROW(df)+1
  #df_updt$parent<-as.integer(df_updt$parent)
  #df_updt$node<-as.integer(df_updt$node)
  #df_updt<-add_tip(df_updt, NROW(df)+1, paste0(tip_label, "_2"), sp_id, N_of_tips) ##5 is the number of deleted tips
  #return(df_updt)
#}



#df_new<-split_tip(df_new, "9031", new_L-5, "CHK")
#df_new<-split_tip(df_new, "638290", new_L-5+1, "GCF")

df_new$label<-df_new$label_sp

df_new<-df_new[, c("parent", "node", "branch.length", "label")]


tree_clean<-as.phylo(df_new)

write.tree(tree_clean, paste0("01_basicStats/01.6_ITOL/tree_species_as_tips.phy"))


write.csv(df_new, "01_basicStats/01.6_ITOL/tree_species_as_tips.csv", row.names = F)


###saving the tree with the ordered species
d <- subset(ggtree::fortify(tree_clean), isTip)

sp_full_order <- with(d, label[order(y, decreasing=T)])

##adding the "jawless vertebrata"
stats_annot$color_class <- as.character(stats_annot$color_class)
                                                       
stats_annot[species == "JL"]$color_class <-"Jawless_vertebrate" 
                                                       
stats_annot$color_class <- factor(stats_annot$color_class, 
        levels = c('Invertebrata','Jawless_vertebrate', 'Actinopteri','Chondrichthyes','Amphibia','Reptilia', 'Aves', 'Marsupialia', 'Mammalia'))
                                                       
##now we need to rearrange them, so that we start with invertebrata
sp_full_order_colors <- sapply(sp_full_order, 
    function(x) stats_annot[stats_annot$species == x, ]$color_class[[1]])

sp_correct_order <- as.character(unlist(sapply(names(class_colors), 
function(x) names(sp_full_order_colors[sp_full_order_colors == x]))))



write.table(as.data.frame(sp_correct_order), file.path(meta_dir, 
                            "species_list_ordered_2021.txt"), 
                                            row.names = F, quote=F, col.names = F)

