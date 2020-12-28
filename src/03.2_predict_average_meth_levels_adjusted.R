source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/03.0_predict_meth_init.R"))

library(phylolm)
library(tidytree)

#N_boot <- 10 #dev

# bootstrapping parameters
N_boot<-100 
group_colors<-c("no prior" ="#66c2a5", "phyl. groups"="#fc8d62", "taxonomy tree"="#8da0cb")

### loading the data (generated at 03.1_prediction_feautures)
summary_df<-read.csv(file.path(analysis_dir, "01_basicStats/feature_summary_filtered.tsv"), sep = "\t")

summary_df$color_class<-factor(summary_df$color_class, levels=names(class_colors))

##target to predict 
target="mean_meth"

#groups of features
seq_features<-c("coveredCpGs", "CpG_O", "CG_freq", "CpG_OE")
cpg_island_features<-c("Gardiner_Garden", "Takai_Jones")

kmers1<-c("A", "C", "G") # freq of T depends on others
kmers2<-c("AA", "AC", "AG", "AT", "CA", "CC", "CpG_O", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG"  )
kmers3<-colnames(summary_df)[c(28:89)]

##list of features for further use
all_features<-list("sequence features"=seq_features, "CpG criteria"=cpg_island_features,
                   "1mers"=kmers1, "2kmers"=kmers2, "3kmers"=kmers3)
                   

## transferring the color_class into the dummy variable:
d <- dummy(summary_df$color_class, sep="_")[, -1]
colnames(d) <- paste0("is_", names(class_colors)[-1])
df <- cbind(summary_df, d)
phyl_features <- colnames(d)
all_f <- unique(unlist(all_features))


## using the Yin and Fan (2001) label Wherry Formula-1 as the standart lm() function does
adjust_r2 <- function(df, r){
    
    k = ncol(df)
    N = NROW(df)
    
    return( 1 - ((1- r)*(N - 1))/(N - k - 1))
    
}


##baseline linear model
get_reg_baseline<-function(df, idx, target = "mean_meth", use_phylo = F){
  if(use_phylo){
    fit <- lm(df[idx,target]~., data = data.frame(scale(df[idx,phyl_features], scale=T, center=T)), subset=idx)
  }else{
    fit <- lm(df[idx,target]~1, data = data.frame(scale(df[idx,phyl_features], scale=T, center=T)), subset=idx)
  }
  r <- rr2::R2.pred(fit)
  r_adj <- adjust_r2(df[idx,], r)
    
  return(c(AIC(fit), r, r))
}


##function for normal linear models bootstrap

get_reg<-function(df, idx, target = "mean_meth",features, use_phylo = F){
    
  if(use_phylo){
      
    fit<-lm(df[idx,target]~., data = data.frame(scale(df[idx,c(features, phyl_features)], scale=T, center=T)),
            subset=idx)
      
    r <- rr2::R2.pred(fit)
    r_adj <- adjust_r2(df[idx,c(features, phyl_features)], r)
  
  }else{
      
    fit<-lm(df[idx,target]~., data = data.frame(scale(df[idx,features], scale=T, center=T)), subset=idx)
  
    r <- rr2::R2.pred(fit)
    r_adj <- adjust_r2(df[idx,features], r)
  }
  
    
  return(c(AIC(fit), r_adj, r))
}


##bootstrapping the baseline linear model (no features)
b <- boot(df, statistic = get_reg_baseline, R = N_boot, target = target)

simple_lm_booted_lm0 <- c(colMeans(b$t), apply(b$t, 2, sd))
print(simple_lm_booted_lm0)

##bootstrapping the linear model (no features)
simple_lm_booted<-lapply(all_features, function(x) {b<-boot(df, statistic = get_reg, R = N_boot, target = target, features = x); return(c(colMeans(b$t), apply(b$t, 2, sd)))})

## joining them together
simple_lm <- cbind(as.data.frame(simple_lm_booted), as.data.frame(simple_lm_booted_lm0))

print(head(simple_lm))

## transfrerring the data frame into the readble format
colnames(simple_lm)[[6]]<-"baseline"
simple_lm <- as.data.frame(t(simple_lm))

colnames(simple_lm) <- c("mean_AIC",  "mean_R2.adjusted", "mean_R2", "sd_AIC", "sd_R2.adjusted", "sd_R2")
simple_lm["type"] <- "no prior"

row.names(simple_lm) <- c("sequence features", "CpG criteria", "1mers", "2mers", "3mers", "no features")
simple_lm["feature_set"] <- row.names(simple_lm)

print("simple done")
save(simple_lm, file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/model_simple.RData")

## bootstrappping phylo-group based lm - baseline
b <- boot(df, statistic = get_reg_baseline, R = N_boot, target =target, use_phylo=T)
phylo_groups_lm_booted_lm0 <- c(colMeans(b$t), apply(b$t, 2, sd))


## bootstrapping phylo-based - different features
phylo_groups_lm_booted<-lapply(all_features, function(x) { b <- boot(df, statistic = get_reg, R = N_boot,
                                                                     target = target, features = x, use_phylo=T); 
                                return(c(colMeans(b$t), apply(b$t, 2, sd)))})

phylo_groups <- cbind(as.data.frame(phylo_groups_lm_booted), as.data.frame(phylo_groups_lm_booted_lm0))

## transferring into a readble format
colnames(phylo_groups)[[6]]<-"baseline"
phylo_groups <- as.data.frame(t(phylo_groups))

colnames(phylo_groups) <-  c("mean_AIC",  "mean_R2.adjusted", "mean_R2", 
                             "sd_AIC", "sd_R2.adjusted", "sd_R2")
phylo_groups["type"] <- "phyl. groups"

row.names(phylo_groups) <- c("sequence features", "CpG criteria", "1mers", 
                             "2mers", "3mers", "no features")

phylo_groups["feature_set"] <- row.names(phylo_groups)

print("phylo_groups done")
save(phylo_groups, file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/model_phylo_groups.RData"))

#######tree
tree <- read.tree(paste0(analysis_dir, "/99.2_ITOL/tree_species_as_tips.phy"))

## now the proper phylogenetic regression using the phylolm package:
## here we are implementing a bit different bootstrapping approach, without resamplinng, 
#due to the fact we can not perform resampling on the taxonomy tree

     
phyl_reg_baseline <- function(df, idx, target, tree){
  
  tips.to.drop <- setdiff(tree$tip.label, df[idx, ]$species)
  tree_sub <- drop.tip(tree, tips.to.drop)
  df_sub <- df[idx, ]
  row.names(df_sub) <- df_sub$species
  
    d <- data.frame(df_sub)
  phyl_m <- phylolm(df_sub[, target]~1, data=d, phy = tree_sub, boot=0)
  
  r <- rr2::R2.pred(mod = phyl_m, phy=tree_sub)
  r_adj <- adjust_r2(d[,c()], r)
  
  return(c(phyl_m$aic, r_adj, r))

  }
     
phyl_reg<-function(df, idx, target, features, tree){
  
  tips.to.drop<-setdiff(tree$tip.label, df[idx, ]$species)
  tree_sub<-drop.tip(tree, tips.to.drop)
  df_sub<-df[idx, ]
  row.names(df_sub)<-df_sub$species
  d<-data.frame(scale(df_sub[, features], center=T, scale=T))
  
  phyl_m<-phylolm(df_sub[, target]~., data=d, phy = tree_sub, boot=0)

  r <- rr2::R2.pred(mod = phyl_m, phy=tree_sub)
  #print(r)
  r_adj <- adjust_r2(d, r)
  
  return(c(phyl_m$aic, r_adj, r))
  
}


phyl_reg_boot<-function(df, target, features, tree, R){
  if (length(features)==0){
    ans <- mclapply(seq(R), function(x) {ind<-sample(row.names(df), 0.9*NROW(df)); print(x);
      r <- phyl_reg_baseline(df, ind, target, tree); return(r)})
      
  }else{
      
    ans<-mclapply(seq(R), function(x) {ind<-sample(row.names(df), 0.9*NROW(df));print(x);
                        r <- phyl_reg(df, ind, target, features, tree); return(r)})
  }
  print(ans)
  ans_df<-Reduce(rbind, ans)
  return(c(colMeans(ans_df), apply(ans_df, 2, sd)))
}
     
## running the baseline     
ans <- lapply(seq(N_boot), function(x) { ind<-sample(row.names(df), 0.9*NROW(df)); 
      r <- phyl_reg_baseline(df, ind, target, tree); 
                                        return(r)})
ans_df<-Reduce(rbind, ans)
     
phyl_boot_bl <- c(colMeans(ans_df), apply(ans_df, 2, sd))
save(phyl_boot_bl, 
     file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/model_phyl_boot_bl.RData"))

## running the bootsrtap on the tree
     
phyl_boot <- lapply(all_features, function(x) phyl_reg_boot(df, target, x, tree, 2))
save(phyl_boot, 
     file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/model_phyl_boot.RData"))
                    
                    
## joining and transforming

tax_tree <- cbind(as.data.frame(phyl_boot), phyl_boot_bl)

colnames(tax_tree)[[6]]<-"baseline"
tax_tree <- as.data.frame(t(tax_tree))

colnames(tax_tree) <- c("mean_AIC",  "mean_R2.adjusted", "mean_R2", 
                             "sd_AIC", "sd_R2.adjusted", "sd_R2")

row.names(tax_tree) <- c("sequence features", "CpG criteria", "1mers", 
                             "2mers", "3mers", "no features")
tax_tree["feature_set"] <- row.names(tax_tree)
                    
tax_tree["type"] <- "taxonomy tree"
print("tree done")
                    
####### joining and visualising
                    
df <- rbind(simple_lm, rbind(phylo_groups, tax_tree))
                    
df[df$feature_set == "sequence features",]$feature_set <- "CG comp."
df[df$feature_set == "CpG criteria",]$feature_set <- "CGI freq."
df$feature_set<-factor(df$feature_set, 
                       levels=c("no features", "CG comp.", "CGI freq.",
                                "1mers", "2mers", "3mers"))
write.table(df, 
            paste0(analysis_dir, "/02_predict_meth/02.6_predictability/features_adjustedR2.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE )     

## baseline-based models visualized
                    
ggplot(df[df$type != "no prior", ], aes(x = feature_set, y=mean_R2.adjusted, fill = mean_AIC))+
  geom_bar(stat="identity", position = "dodge") + rotate_labels()+
  geom_errorbar(aes(ymin=mean_R2.adjusted-sd_R2.adjusted, ymax=mean_R2.adjusted+sd_R2.adjusted), 
                width=.2,position=position_dodge(.9)) + facet_wrap(~type) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits = c(3300, 5000)) + ylim(c(0,1)) + 
  geom_hline(yintercept = df[df$feature_set=="3mers" & df$type=="no prior","mean_R2.adjusted"], linetype="dashed")+
 labs(x="",  y="mean R2 (adjusted)", fill= "mean AIC")
                    
ggsave(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/", target,"_phylfeatures_adjustedR2_", N_boot, ".pdf"), width = 6, height = 4)
                    
                    
## no prior model
                    
ggplot(df[df$type == "no prior", ], aes(x = feature_set, y=mean_R2.adjusted, fill = mean_AIC))+
  geom_bar(stat="identity", position = "dodge") + rotate_labels()+
  geom_errorbar(aes(ymin=mean_R2.adjusted-sd_R2.adjusted, ymax=mean_R2.adjusted+sd_R2.adjusted), 
                width=.2,position=position_dodge(.9)) + facet_wrap(~type) +
  scale_fill_gradient(low="#56B1F7", high="#132B43", limits = c(3300, 5000)) + ylim(c(0,1)) + 
  geom_hline(yintercept = df[df$feature_set=="3mers" & df$type=="no prior","mean_R2.adjusted"], linetype="dashed")+
 labs(x="",  y="mean R2 (adjusted)", fill= "mean AIC")
                    
ggsave(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/", target,"_noprior_adjustedR2_", N_boot, ".pdf"), width = 6, height = 4)

## comparing AICs:
cat(paste0(c("prior", "first", "second", "rel.lik."), sep = "\t"), 
    file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/rel_likelihood.txt"))

e = exp((df[df$type == "no prior" & df$feature_set == "3mers",]$mean_AIC - df[df$type == "no prior" & df$feature_set == "2mers",]$mean_AIC)/2) 
                    
cat(paste0(c("\nno prior", "2mers", "3mers", e), sep = "\t"), 
    file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/rel_likelihood.txt"), append = TRUE)
                    
e = exp((df[df$type == "phyl. groups" & df$feature_set == "3mers",]$mean_AIC - df[df$type == "phyl. groups" & df$feature_set == "no features",]$mean_AIC)/2) 
                    
cat(paste0(c("\nphyl. groups", "no features", "3mers", e), sep = "\t"), 
    file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/rel_likelihood.txt"), append = TRUE)
                    
e = exp((df[df$type == "taxonomy tree" & df$feature_set == "3mers",]$mean_AIC - df[df$type == "taxonomy tree" & df$feature_set == "no features",]$mean_AIC)/2) 

cat(paste0(c("\ntaxonomy tree", "no features", "3mers", e), sep = "\t"), 
    file = paste0(analysis_dir, "/02_predict_meth/02.6_predictability/rel_likelihood.txt"), append = TRUE)
                    
                    