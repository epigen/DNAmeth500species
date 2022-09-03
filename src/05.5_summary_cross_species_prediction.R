source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(treeio)
library(circlize)
library(ComplexHeatmap)
library(ggtree)

wd = file.path(analysis_dir, "05_predict_meth", "05.2_test_on_other_species")
setwd(wd)

dir.create("summary")
all_sp=unique(stats_annot$species)

all_id<-lapply(all_sp, function(x) stats_annot[stats_annot$species==x]$ncbi_id[[1]])
all_class<-sapply(all_sp, function(x) as.character(stats_annot[stats_annot$species==x]$color_class[[1]]))
all_sp_id<-setNames(all_sp, all_id)
all_class_sp<-setNames(all_class, all_sp)
##executing once - collecting prediction data

##aucs need to be load for the class - specific mss score

stats_files=system(paste0("ls ", analysis_dir, "/05_predict_meth/05.2_test_on_other_species/screen/*/all_aucs.csv"),intern=TRUE)
                  
sp_list <- unlist(lapply(stats_files, function(x) strsplit(x, "/")[[1]][11] ))

##check if there are some missing:
missed <- setdiff(sp_df$species, sp_list)
length(missed)
missed

if (length(missed) > 5){
  write.table(as.data.frame(missed), file.path(Sys.getenv("CODEBASE"), "compEpi/meta/to_run_cross.txt"), 
              col.names = F, quote = F, row.names = F)
}

##collecting the stats data
auc_file <- file.path("summary", "all_aucs_full.csv")
if(!(file.exists(auc_file))){
    
    all_aucs = data.table()
    for (i in 1:length(stats_files)){
        df <- fread(stats_files[i])
        df$train_species <- sp_list[i]
        all_aucs <- rbind(all_aucs, df)
    }
    
    all_aucs <- left_join(all_aucs, sp_df, by = c("type" = "species"))
    
    colnames(all_aucs)[5] <- "color_class_test"
    colnames(all_aucs)[6] <- "group_test"
    
    all_aucs <- left_join(all_aucs, sp_df, by = c("train_species" = "species"))
    colnames(all_aucs)[7] <- "color_class_train"
    colnames(all_aucs)[8] <- "group_train"
    all_aucs[ifRand == "noRandTrain", ifRand := "noRand",]
  all_aucs <- all_aucs[!all_aucs$type %in% missed,]
  my_wt(all_aucs, auc_file)
}

all_aucs <- fread(auc_file)

#tree with the species we explore, saving separately
tree <- read.tree(file.path(analysis_dir, "01_basicStats",
                            "01.6_ITOL", "tree_species_as_tips_2021_2.phy"))
p_tree <- ggtree(tree)
tre1e <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% unique(all_aucs$train_species)])
p_tree <- ggtree(tre1e, branch.length = "none" )

sp_df <- as.data.frame(sp_df)
row.names(sp_df) <- sp_df$species
#sapply(tre1e$tip.label, function(x) sp_df[x,]$color_class)

df <- data.frame(label = tre1e$tip.label,
                 color_class=sapply(tre1e$tip.label, 
                                    function(x) sp_df[x,]$color_class))
                                    
#row.names(df) <- df[,1]

p1 <- ggtree(tre1e, branch.length = "none") %<+% df + geom_tiplab(size = 1, aes(color=color_class)) 

pdf(paste0(analysis_dir,   "/05_predict_meth/05.2_test_on_other_species/summary/tree_full_w_labels.pdf"), height = 18, width = 9)
p1 + scale_color_manual(values = class_colors)
dev.off()

pdf(paste0(analysis_dir, 
"/05_predict_meth/05.2_test_on_other_species/summary/tree_full_w_labels_and_nodes.pdf"), height = 18, width = 9)
print(p1 +  geom_nodelab(size = 2) + scale_color_manual(values = class_colors))
dev.off()

## now we have to rotate, for this we also need the node labels##first, lets try rotating "Amniota"

tre1e_annot <- tre1e %>% as_tibble
tre1e_annot[tre1e_annot$label == "Amniota", ] #585 - big rotate
tre1e_annot[tre1e_annot$label == "Theria", ] #586 ##marsupiala and mammalia
tre1e_annot[tre1e_annot$label == "Sauria", ] #678

tre1e_rotated <- ggtree(tre1e, branch.length = "none") %>% rotate(585)  %>%
                                    rotate(586) %>% rotate(678) 

pdf(paste0(analysis_dir, "/05_predict_meth/05.2_test_on_other_species/summary/tree_full_right_order.pdf"), height = 18, width = 9)
tre1e_rotated %<+% df + geom_tiplab(size = 1, aes(color=color_class)) + scale_color_manual(values = class_colors)
dev.off()


### saving the correct order
td_out <- tre1e_rotated$data
td_out <- dplyr::arrange(td_out, y)
sp_full_order <- td_out %>% filter(isTip) %>% pull(label)
                                    
sp_full_order <- fread(file.path(meta_dir, "species_list_ordered_2021_2.txt"), header = F)$V1
                                    



all_auc_m<-all_aucs %>% 
  filter(ifRand=="noRandTest"|ifRand=="noRand" ) %>%
  select(c("train_species", "type", "auc")) %>% 
  spread(type, auc)

##transferring in the matrix
all_auc_m <- as.data.frame(all_auc_m)
row.names(all_auc_m) <- all_auc_m$train_species

sp_full_order <- sp_full_order[sp_full_order %in% all_auc_m$train_species]

##annotation the heatmap
col_annot <- data.frame( sapply(sp_full_order, function(x) as.character(all_class[x])))
colnames(col_annot) <- c("class")
                                
                                
all_auc_m <- as.matrix(all_auc_m[sp_full_order, sp_full_order])

write.table(all_auc_m, "summary/all_aucs_full_matrix.csv", sep = "\t", quote = FALSE) ##need to save rownames!

##plotting the annotation
                                
col_fun_num = colorRamp2(c(0, 2000), c("white", "darkgrey"))
ha = columnAnnotation(df = col_annot, col = list("class" = class_colors), 
                      show_legend=c(F), height = unit(0.2, "cm"))

ra = rowAnnotation(df = col_annot, col = list("class" = class_colors), show_legend=c(F), width = unit(0.2, "cm"))
                                
col_fun = colorRamp2(c(0, 0.4,0.5,0.6, 1), c("#2b83ba","#abd9e9", "#ffffbf", "#fdae61","#d7191c"))

Heatmap(all_auc_m,
        cluster_columns = F, cluster_rows = F)

pdf(paste0(analysis_dir, "/05_predict_meth/05.2_test_on_other_species/summary/aucs_full_small.pdf"), height=4, width=5)
ra + Heatmap(all_auc_m,
        cluster_columns = F, cluster_rows = F,
        column_order = sp_full_order, row_order = sp_full_order, 
        show_column_names = F, show_row_names = F,
        top_annotation = ha, name="ROC-AUC", col = col_fun, use_raster = T, raster_device = "png")
dev.off()

#summary boxplots

all_auc_test <- all_aucs[all_aucs$ifRand == "noRandTest",]

                                all_auc_test$color_class_test <- factor(all_auc_test$color_class_test, levels = names(class_colors))
                                
                                
all_auc_test$color_class_train <- factor(all_auc_test$color_class_train, levels = names(class_colors))
                                
all_auc_test[color_class_test=="Jawless_vertebrate", group_test := "Jl.vb.",] 
                                
all_auc_test$group_test <- factor(all_auc_test$group_test, levels = class_short)                                
                                
all_auc_test %>%
  group_by(.dots = c("color_class_test", "color_class_train")) %>%
  summarize(n=n())

#all_auc_test[color_class_test=="Jawless_vertebrate", color_class_test := "Invertebrata",]
#all_auc_test[color_class_train=="Jawless_vertebrate", color_class_train := "Invertebrata",]

#all_auc_test[group_test=="Jl.vb.", group_test := "Inv.",]
#all_auc_test[group_train=="Jl.vb.", group_train := "Inv.",]
                                
pdf(paste0(analysis_dir, "/05_predict_meth/05.2_test_on_other_species/summary/aucs_by_class_lines.pdf"), height=4, width=12)
ggplot(all_auc_test[color_class_train!="Jawless_vertebrate"], aes(x = group_test, y = auc, fill = color_class_test)) + geom_boxplot(outlier.shape = 21)  + 
  scale_fill_manual(values=c(class_colors)) + facet_wrap(~color_class_train, ncol = 4) + 
  theme(legend.position = "None") + labs(x = "class, tested on", y = "ROC-AUC") + 
                      rotate_labels()+  
stat_summary(fun.data = give.n,fun.args = c(y=0.1), geom = "text",size=4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", size= 0.1) + 
  geom_hline(yintercept = 0.25, linetype = "dashed", size= 0.1) + 
  geom_hline(yintercept = 0.75, linetype = "dashed", size= 0.1)
dev.off()

                        
all_aucs %>%
  filter(ifRand == "noRand") %>%
  group_by(color_class_train) %>%
  summarise(mean_auc = mean(auc))

all_aucs %>%
  group_by(ifRand) %>%
  summarise(median_auc = median(auc))


#####NOT USED FOR THE FINAL FIGURE#####################
#only_aucs<-lapply(aucs, function(a) {a_short<-a[(a$ifRand=="noRandTest"|a$ifRand=="noRandTrain" ), c("type", "auc")]; sp<-a[a$ifRand=="noRandTrain", "type"];colnames(a_short)<-c("species", as.character(sp)); return(a_short)})
#full_auc_m<-Reduce(full_join, only_aucs)
#write.csv(full_auc_m, paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/all_aucs_matrix.csv"), quote = F)

###reading in the 
full_auc_m<-read.csv(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/all_aucs_matrix.csv"), row.names = 1)
row.names(full_auc_m)<-full_auc_m$species
full_auc_mat<-t(as.matrix(full_auc_m[, -1]))


tree<-read.tree(paste0(analysis_dir, "/99.2_ITOL/tree_species_as_tips.phy"))

p_tree<-ggtree(tree)
tre1e<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% unique(all_aucs$train_species)])
p_tree<-ggtree(tre1e, branch.length = "none" )
p_tree

pdf(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/tree_full.pdf"), height = 18, width = 9)
print(p_tree)
dev.off()

d = fortify(tre1e)
d = subset(d, isTip)
sp_full_order<-with(d, label[order(y, decreasing=T)])

##now the subelection
tre1e<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% colnames(full_auc_mat)])
p_tree<-ggtree(tre1e, branch.length = "none" )
p_tree

pdf(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/tree_train.pdf"), height = 6, width = 3)
print(p_tree)
dev.off()

d = fortify(tre1e)
d = subset(d, isTip)
sp_selected_order<-with(d, label[order(y, decreasing=T)])

col_annot_row<-data.frame( sapply(sp_selected_order, function(x) as.character(all_class[x])))
col_annot_col<-data.frame( sapply(sp_full_order, function(x) as.character(all_class[x])))

colnames(col_annot_col)<-c("class")
colnames(col_annot_row)<-c("class")

ha = columnAnnotation(df = col_annot_col, col = list("class" = class_colors), show_legend=c(F))
ra = rowAnnotation(df = col_annot_row, col = list("class" = class_colors), show_legend=c(F))

col_fun = colorRamp2(c(0, 0.4,0.5,0.6, 1), c("#2b83ba","#abd9e9", "#ffffbf", "#fdae61","#d7191c"))

pdf(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/aucs_hm_nocl_cm.pdf"), height=2, width=4)
Heatmap(full_auc_mat[sp_selected_order, sp_full_order ],
        cluster_columns = F, cluster_rows = F,
        column_order = sp_full_order, row_order = sp_selected_order, 
        show_column_names = F, show_row_names = F,
        top_annotation = ha, name="AUC", col = col_fun, use_raster = T, raster_device = "png")+ra
dev.off()

##for actinopteri - median auc on test data
median_auc <- all_aucs %>%
  filter(class == "Actinopteri" & ifRand=="noRandTest" & train_class == "Actinopteri") %>%
  group_by(train_species) %>%
  summarise(median_auc = median(auc))

median_auc <- as.data.frame(median_auc)  
median_auc <- inner_join(median_auc, all_aucs[all_aucs$ifRand == "noRandTrain",c("train_species", "auc")], by = "train_species")

median_full_auc <- all_aucs %>%
  filter(ifRand=="noRandTest" & train_class == "Actinopteri") %>%
  group_by(train_species) %>%
  summarise(median_full_auc = median(auc))

median_auc <- inner_join(median_auc,median_full_auc, by = "train_species")

ggplot(median_auc, aes( x = median_auc, y=median_full_auc)) + geom_point() + geom_text_repel(aes(label = train_species))

write.csv(median_auc, paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/median_auc_Actinopteri.csv"))
all_aucs$class<- factor(all_aucs$class, levels = names(class_colors))
all_aucs$train_class<- factor(all_aucs$train_class, levels = names(class_colors))

all_pred_aucs <- all_aucs[all_aucs$ifRand == "noRandTest", c("type", "auc", "class", "train_species", "train_class")]
all_self_aucs <- all_aucs[all_aucs$ifRand == "noRandTrain", c( "auc", "train_species", "train_class")]
med_self_auc <- as.data.frame(setDT(all_self_aucs)[, j = list(median_self_auc = median(auc)), by = train_class])

colnames(all_self_aucs)[[1]] <- "self_auc"

left_join(all_pred_aucs, all_self_aucs, by = c("train_species", "train_class"))

ggplot(all_aucs[all_aucs$ifRand == "noRandTest", ], aes(x = class, y = auc )) + 
  geom_boxplot(aes(fill = class), outlier.shape = 21, varwidth = T)  +rotate_labels() + 
  scale_fill_manual(values = class_colors) + 
  geom_hline(data = med_self_auc, aes(yintercept = median_self_auc), linetype = "dashed")+
  facet_wrap(~train_class, ncol = 2)

ggplot(all_pred_aucs[all_pred_aucs$train_class == "Mammalia", ], aes(x = class, y = auc )) + 
  geom_boxplot(aes(fill = class), outlier.shape = 21, varwidth = T)  +rotate_labels() + 
  scale_fill_manual(values = class_colors)+geom_hline( yintercept  = 0.7662280)

##summarizing all in a bubble heatmap that reprsents the median value and the std of the AUCs in the data
all_aucs_classes<- all_aucs %>%
  filter(ifRand=="noRandTest") %>%
  group_by(.dots=c("train_class", "class")) %>%
  summarize(auc_median=median(auc), auc_std=sd(auc), auc_mean=mean(auc))
all_aucs_classes$train_class<-factor(all_aucs_classes$train_class, levels=names(class_colors))
all_aucs_classes$class<-factor(all_aucs_classes$class, levels=names(class_colors))

p_bubble<-ggplot(all_aucs_classes, aes(x=train_class, y=class, size=auc_median, color=auc_std))+
  geom_point()+rotate_labels()+scale_colour_gradient(low = "#fc8d62", high ="grey" ) + 
  scale_size_continuous(name="median\n ROC-AUC", range = c(1,8), 
          breaks=c(0.6, 0.9), limits=c(0.5, 1))+
  labs(colour="std of the\n ROC-AUC", x="group, trained on", y="group, tested on") +
  theme(axis.text =element_text(size = 15), axis.title =element_text(size = 15) )
p_bubble


pdf(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/aucs_bubble.pdf"), height=6, width=8)
print(p_bubble)
dev.off()




## mean methylation score
MSS_by_class<-function(auc){
  self_auc<-auc[auc$ifRand=="noRandTrain","auc"]
  med_auc<-auc %>% 
    filter(ifRand=="noRandTest") %>% 
    group_by(class) %>%
    summarize(med_auc=median(auc), species=first(train_species), train_class=first(train_class))
  med_auc$mss<-self_auc/med_auc$med_auc
  return(med_auc)
}

mss_by_class_list<-rbindlist(lapply(aucs, MSS_by_class))
mss_by_class_list$isTrainedOn<-mss_by_class_list$class==mss_by_class_list$train_class
mss_by_class_list$class<-factor(mss_by_class_list$class, levels=names(class_colors))

p_mss<-ggplot(mss_by_class_list, aes(x=class, y=mss, color=isTrainedOn, fill = class))+
  geom_boxplot(outlier.shape = 21)+ylab("Median Specificity Score")+rotate_labels()+
  geom_hline(yintercept=1, linetype="dashed")+facet_wrap(~train_class, nrow = 2)+
  scale_color_manual(values = c("TRUE"="red", "FALSE"="black"))+scale_fill_manual(values = class_colors)

pdf(paste0(analysis_dir, "/02_predict_meth/02.2_test_on_other_species/summary/MSS.pdf"), height=10, width=15)
print(p_mss)
dev.off()


##inverded - identify
s<-colMeans(full_auc_mat[,sp_selected_order])
inv<-names(s[s<0.45])


s<-rowMeans(full_auc_mat[sp_full_order,])
full_inv<-names(s[s<0.5])
