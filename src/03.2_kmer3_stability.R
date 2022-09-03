source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/03.0_predict_meth_init.R"))
library(phylolm)
library(dummies)
library(MASS)
library(boot)
library(purrr)
library(magrittr)
library(corrplot)
library(caret)
library(LiblineaR)


wd = file.path(analysis_dir, "03_prediction_lm", "03.2_kmer_stability")
dir.create(wd)
setwd(wd)

summary_df <- fread(file.path(analysis_dir, "01_basicStats/01.10_all_seq_features/feature_summary_filtered.tsv"), sep = "\t")
summary_df[color_class=="Jawless_vertebrate", color_class:="Invertebrata",]
summary_df$color_class <- factor(summary_df$color_class, levels=names(class_colors))
N_boot <- 100 #bootstrapping values
#N_boot <- 10
## identifying sets of features, we are interested at
target="mean_meth"

#groups of features
seq_features <- c("coveredCpGs", "CpG_O", "CG_freq", "CpG_OE")
cpg_island_features <- c("Gardiner_Garden", "Takai_Jones")
kmers1 <- c("A", "C", "G") # freq of T depends on others
kmers2 <- c("AA", "AC", "AG", "AT", "CA", "CC", "CpG_O", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG"  )
kmers3 <- colnames(summary_df)[c(28:90)]
  

## transferring the color_class into the dummy variable:
d <- dummy(summary_df$color_class, sep="_")[, -1]
colnames(d) <- paste0("is_", names(class_colors)[c(-1,-2)])
df <- cbind(summary_df, d)

phyl_features <- colnames(d)

all_features <- list("meth"=seq_features, "cpgs"=cpg_island_features,
                  "kmers1"=kmers1, "kmers2"=kmers2, "kmers3"=kmers3)

##stepwise selection in both directions, adding and deleting
both_ranks<-function(data, indices, target="mean_meth",features, use_phylo=F){
  d=data[indices,]
  d_scaled<-as.data.frame(scale(d[,features], center = TRUE, scale = TRUE))
  is.na(d_scaled)<-sapply(d_scaled, is.infinite)
  d_scaled[is.na(d_scaled)]<-0
  
  if(use_phylo){
    lm0<-lm(d[, target]~is_Actinopteri+is_Chondrichthyes+is_Amphibia+is_Reptilia+
              is_Aves+is_Marsupialia+is_Mammalia,  data=d_scaled)
  }else{
    lm0<-lm(d[, target]~1,  data=d_scaled)
  }
  
  lm_full<-formula(lm(d[, target]~., data=d_scaled))
  
  model<-stepAIC(lm0, scope=lm_full, direction = "both", k=log(nrow(d)), trace=0)
  
  #anova shows steps taken, that we transform into a matrix of features (on which step what happened)
  table<-model$anova
  table<-table[-1, ]
  table["rank"]<-as.integer(row.names(table))
  table["name"]<-unlist(lapply(table$Step, function(x) paste0(strsplit(as.character(x), " ")[[1]], collapse="")))
  table<-table[order(table$rank, decreasing = TRUE),]
  table<-table[!duplicated(table$name), ]
  v<-extract2(table, "rank") %>% set_names(table$name)
  
  
  features_s<-c(paste0("+", features),paste0("-",features)) 
  
  missed_names<-setdiff(features_s, names(v))
  v<-c(v, setNames(rep(0, length(missed_names)), missed_names))
  v<-v[order(names(v))]
  ##final weights of the features in the model
  w<-model$coefficients
  missed_names_w<-setdiff(features, names(w))
  w<-c(w, setNames(rep(0, length(missed_names_w)), missed_names_w))
  answer<-c(v[features_s],w)
  return(c(answer))
  #return(w)
}

df <- as.data.frame(df)
boot.3mers<-boot(df, statistic = both_ranks,R=N_boot,target = target, 
                features = kmers3)

#boot.3mers.groups<-boot(df, statistic = both_ranks,R=N_boot,target = target, 
 #                features = c(phyl_features, kmers3), use_phylo=T)

model_stability<-function(model, R, features){
  all_f<-model$t
  colnames(all_f)<-names(model$t0)
  scores<-all_f[, !colnames(all_f) %in% c("(Intercept)",phyl_features, unlist(all_features))]
  models<-as.data.frame(all_f[,colnames(all_f) %in% features])
  scores <- as.data.frame(scores)
  scores["boot_ids"]<-row.names(scores)
  r_m<-melt(scores, id.vars=c("boot_ids"))
  r_m<-as.data.frame(r_m)
  #reforming the output
  r_m["dir"]<-unlist(lapply(r_m$variable, function(x) { y=substr(as.character(x), 1,1); if(y=="+") {return(1)}                                                                      else{return(-1)}}))
  r_m["feature"]<-unlist(lapply(r_m$variable, function(x) substr(as.character(x), 2, nchar(as.character(x)))))
  
  rm_pos<-r_m %>% 
    filter((value!=0)) %>%
    dplyr::group_by(feature, boot_ids) %>%
    dplyr::mutate(value_m=max(value)) %>%
    filter(value==value_m & dir ==1) #%>%# we removed all the events, except the last with the feature
  ##now for all the cases adding the actual feature weight
  models$boot_ids <- row.names(models)
  models_m <- melt(models, id.vars = c("boot_ids"))
  colnames(models_m) <- c("boot_ids", "feature", "feature_weight")
  rm_pos <- left_join(rm_pos, models_m)
  
  rm_pos <- rm_pos %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(pos=sum(feature_weight>0)/R, neg =sum(feature_weight<0)/R, n = n())
  
  rm_neg<-r_m %>% 
    filter((value!=0)) %>%
    group_by(.dots=c("feature", "boot_ids")) %>%
    mutate(value_m=max(value)) %>%
    filter(value==value_m) %>%# we removed all the events, except the last with the feature
    filter(dir==-1) %>%
    group_by(feature) %>%
    summarize(n=n()/R)
  
  return(rm_pos)
}

##extracting the table with feature weights
feat_w <- boot.3mers$t
colnames(feat_w) <- names(boot.3mers$t0)
colSums(feat_w != 0)
feat_w_count <- data.frame(colSums(feat_w > 0), colSums(feat_w < 0), colSums(feat_w != 0))
colnames(feat_w_count) <- c("positive", "negative", "total")
feat_w_count$kmer <- row.names(feat_w_count)
feat_w_count$kmer <- factor(feat_w_count$kmer, levels = feat_w_count$kmer[order(feat_w_count$total, decreasing = T)])
feat_w_count_long <- melt(feat_w_count, measure.vars = c("positive", "negative"))

ggplot(feat_w_count_long[feat_w_count_long$total > 0,], aes(x = kmer, y = value, fill = variable)) + 
  geom_bar(position = "stack", stat = "identity") + rotate_labels()
ggsave("selection_stepwise_NEW.pdf", width = 6, height = 3)

res<-model_stability(boot.3mers, N_boot, kmers3)
res<-res[order(res$pos, decreasing = T),]
res$feature<-factor(res$feature, levels = res$feature)


## all the features that where selected in more than 10 cases from 100
pdf(paste0("selection_stepwise_", target, "_", N_boot, "_filtered_c.pdf"), width = 8, height = 4)
ggplot(res[c(1:30),], aes(x=feature, y=pos))+geom_bar(stat="identity", fill="#66c2a5")+rotate_labels()+
  labs(x="feature", y="stability")
dev.off()


pdf(paste0("selection_stepwise_top_bottom", target, "_", N_boot, "_filtered_c.pdf"), width = 4, height = 4)
ggplot(res[c(1:10),], aes(x=feature, y=pos))+geom_bar(stat="identity", fill="#66c2a5")+rotate_labels()+
  labs(x="", y="stability")+ylim(0,1)
ggplot(res[c(53:63),], aes(x=feature, y=pos))+geom_bar(stat="identity", fill="#66c2a5")+rotate_labels()+
  labs(x="", y="")+ylim(0,1)
dev.off()


                                
                                
### DEPRICATED #######    
res.groups<-model_stability(boot.3mers.groups, N_boot)
res.groups<-res.groups[order(res.groups$n, decreasing = T),]
res.groups$feature<-factor(res.groups$feature, levels = res.groups$feature)

common<-intersect(res[c(1:30),]$feature, res.groups[c(1:30), ]$feature)

col.ticks<-sapply(res[c(1:30), ]$feature, function(x) {if(x %in% common) {return ("red")} else {return ("black")}})                               
col.ticks.groups<-sapply(res.groups[c(1:30), ]$feature, function(x) {if(x %in% common) {return ("red")} else {return ("black")}})
## all the features that where selected in more than 10 cases from 100
pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/selection_stepwise_compare", target, "_", N_boot, "_filtered_c.pdf"), width = 8, height = 4)
ggplot(res[c(1:30),], aes(x=feature, y=n))+geom_bar(stat="identity", fill="#66c2a5")+rotate_labels()+
  labs(x="", y="stability")+ theme(axis.text.x = element_text(colour = col.ticks))+ggtitle("No prior")
ggplot(res.groups[c(1:30),], aes(x=feature, y=n))+geom_bar(stat="identity", fill="#66c2a5")+rotate_labels()+
  labs(x="feature", y="stability")+ theme(axis.text.x = element_text(colour = col.ticks.groups))+ggtitle("groups as prior")
dev.off()

                                
                                
                                
                                
## model with the tree inside!
tree<-read.tree(paste0(analysis_dir, "/99.2_ITOL/tree_species_as_tips.phy"))


tree.stepwise<-function(df, idx, features, target, tree){
  new_df<-df[idx,]
  y<-new_df[, target]
  row.names(new_df)<-new_df$species
  d=data.frame(scale(new_df[, features], scale = T, center = T))
  tips.to.drop<-setdiff(tree$tip.label, new_df$species)
  tree_sub<-drop.tip(tree, tips.to.drop)
  
  fit<-phylostep(y~., data = d[, features], phy = tree_sub, model = "BM", direction="both", k=log(NROW(d)))
  f<-fit$formula
  print(f)
  selected<-strsplit(f, " ")[[1]][sapply(strsplit(f, " ")[[1]], nchar)>1]
  ans<-c(setNames(rep(1,length(selected)), selected), 
         setNames(rep(0, length(features)-length(selected)), setdiff(features, selected)))
  return(ans)
}
f<-c(seq_features, cpg_island_features, kmers3)
b<-lapply(seq(N_boot), function(x) {ind<-sample(row.names(df), 0.9*NROW(df)); 
            a<-tree.stepwise(df, ind, kmers3, target, tree); return(a)})


f<-colSums(Reduce(rbind, b)/N_boot)
tax_res<-data.frame(names(f), f)
colnames(tax_res)<-c("feature", "n")

pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/selection_stepwise_taxtree_", target, 
           "_", N_boot, ".pdf"), width = 8, height = 4)
ggplot(tax_res[c(1:30),], aes(x=feature, y=n))+geom_bar(stat="identity", fill="#8da0cb")+rotate_labels()
dev.off()




