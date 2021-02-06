source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/03.0_predict_meth_init.R"))

library(rr2)
library(phylolm)
library(tidytree)
library(ggtree)


N_boot <- 10 #dev
#N_boot<-100 #bootstrapping
group_colors<-c("no prior" ="#66c2a5", "phyl. groups"="#fc8d62", "taxonomy tree"="#8da0cb")
### loading the data
summary_df<-read.csv(file.path(analysis_dir, "01_basicStats/feature_summary_filtered.tsv"), sep = "\t")

summary_df$color_class<-factor(summary_df$color_class, levels=names(class_colors))

##target
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

##getting a baseline model
get_reg_baseline<-function(df, idx, target = "mean_meth", use_phylo = F, return_metric = "AIC"){
  if(use_phylo){
    fit <- lm(df[idx,target]~., data = data.frame(scale(df[idx,phyl_features], scale=T, center=T)), subset=idx)
  }else{
    fit <- lm(df[idx,target]~1, data = data.frame(scale(df[idx,phyl_features], scale=T, center=T)), subset=idx)
  }
  if(return_metric == "AIC") ans = AIC(fit)
  else if(return_metric == "var") ans = summary(fit)$adj.r.squared
  else ans = 0
  return(c(AIC(fit), summary(fit)$adj.r.squared))
}

##bootstrapping normal linear models
get_reg<-function(df, idx, target = "mean_meth",features, use_phylo = F, return_metric = "AIC"){
  if(use_phylo){
    fit<-lm(df[idx,target]~., data = data.frame(scale(df[idx,c(features, phyl_features)], scale=T, center=T)),
            subset=idx)
  }else{
    fit<-lm(df[idx,target]~., data = data.frame(scale(df[idx,features], scale=T, center=T)), subset=idx)
  }
  if(return_metric == "AIC") ans = AIC(fit)
  else if(return_metric == "var") ans = summary(fit)$adj.r.squared
  else ans = 0
  return(c(AIC(fit), summary(fit)$adj.r.squared))
}

b <- boot(df, statistic = get_reg_baseline, R = N_boot, target = target, return_metric = "var")

#simple_lm_booted_lm0 <- c(colMeans(b$t)[[3]], apply(b$t, 2, sd)[[3]])
simple_lm_booted_lm0 <- c(colMeans(b$t), apply(b$t, 2, sd))
                          
#simple_lm_booted<-lapply(all_features, function(x) {b<-boot(df, statistic = get_reg, R = N_boot,
 #                                                           target = target, features = x); 
  #                                  return(c(colMeans(b$t)[[3]], apply(b$t, 2, sd)[[3]]))})
simple_lm_booted<-lapply(all_features, function(x) {b<-boot(df, statistic = get_reg, R = N_boot,
                                                            target = target, features = x, return_metric = "var"); 
                                            return(c(colMeans(b$t), apply(b$t, 2, sd)))})

simple_lm <- cbind(as.data.frame(simple_lm_booted), as.data.frame(simple_lm_booted_lm0))

colnames(simple_lm)[[6]]<-"baseline"
simple_lm <- as.data.frame(t(simple_lm))

colnames(simple_lm) <- c("mean_AIC",  "mean_R2.adjusted","sd_AIC", "sd_R2.adjusted")
simple_lm["type"] <- "no prior"
row.names(simple_lm) <- c("sequence features", "CpG criteria", "1mers", "2mers", "3mers", "no features")
simple_lm["feature_set"] <- row.names(simple_lm)
print("simple done")

## bootstrappping phylo-group based lm
b <- boot(df, statistic = get_reg_baseline, R = N_boot, target =target, use_phylo=T)

#phylo_groups_lm_booted_lm0<-c(colMeans(b$t)[[3]], apply(b$t, 2, sd)[[3]])

phylo_groups_lm_booted_lm0 <- c(colMeans(b$t), apply(b$t, 2, sd))

#phylo_groups_lm_booted<-lapply(all_features, function(x) {b<-boot(df, statistic = get_reg, R = N_boot,target = target, features = x, use_phylo=T); 
                                   # return(c(colMeans(b$t)[[3]], apply(b$t, 2, sd)[[3]]))})

phylo_groups_lm_booted<-lapply(all_features, function(x) { b <- boot(df, statistic = get_reg, R = N_boot,
                                                                     target = target, features = x, use_phylo=T); 
                                return(c(colMeans(b$t), apply(b$t, 2, sd)))})

phylo_groups <- cbind(as.data.frame(phylo_groups_lm_booted), as.data.frame(phylo_groups_lm_booted_lm0))

colnames(phylo_groups)[[6]]<-"baseline"
phylo_groups <- as.data.frame(t(phylo_groups))

colnames(phylo_groups) <- c("mean_AIC", "mean_R2.adjusted", "sd_AIC", "sd_R2.adjusted")
phylo_groups["type"] <- "phyl. groups"

row.names(phylo_groups) <- c("sequence features", "CpG criteria", "1mers", "2mers", "3mers", "no features")
phylo_groups["feature_set"] <- row.names(phylo_groups)
print("phylo_groups done")

tree<-read.tree(paste0(analysis_dir, "/99.2_ITOL/tree_species_as_tips.phy"))

## now the proper phylogenetic regression using the phylolm package:
## here we are implementing a bit different bootstrapping approach, without resamplinng, 
#due to the fact we can not perform resampling on the taxonomy tree

phyl_reg<-function(df, idx, target, features, tree){
  
  tips.to.drop<-setdiff(tree$tip.label, df[idx, ]$species)
  tree_sub<-drop.tip(tree, tips.to.drop)
  df_sub<-df[idx, ]
  row.names(df_sub)<-df_sub$species
  d<-data.frame(scale(df_sub[, features], center=T, scale=T))
  
  phyl_m<-phylolm(df_sub[, target]~., data=d, phy = tree_sub, boot=0)
  #return(R2(mod = phyl_m, phy=tree_sub))
  return(phyl_m$aic)
}
packageVersion("phylolm")

phyl_reg_baseline <- function(df, idx, target, tree){
  
  tips.to.drop <- setdiff(tree$tip.label, df[idx, ]$species)
  tree_sub <- drop.tip(tree, tips.to.drop)
  df_sub <- df[idx, ]
  row.names(df_sub) <- df_sub$species
  d <- data.frame(df_sub)
  phyl_m <- phylolm(df_sub[, target]~1, data=d, phy = tree_sub, boot=0)
  #return(R2(mod = phyl_m, phy=tree_sub))
  return(phyl_m$aic)

  }

phyl_reg_boot<-function(df, target, features, tree, R){
  if (length(features)==0){
    ans<-mclapply(seq(R), function(x) {ind<-sample(row.names(df), 0.9*NROW(df));
      r<-phyl_reg_baseline(df, ind, target, tree); return(r)})
  }else{
    ans<-mclapply(seq(R), function(x) {ind<-sample(row.names(df), 0.9*NROW(df));
                        r<-phyl_reg(df, ind, target, features, tree); return(r)})
  }
  print(ans)
  ans_df<-Reduce(rbind, ans)
  return(c(mean(ans_df), sd(ans_df)))
}

phyl_boot<-lapply(all_features, function(x) phyl_reg_boot(df, target, x, tree, N_boot))

phyl_boot_bl<-phyl_reg_boot(df, target, c(), tree, N_boot)

tax_tree<-cbind(as.data.frame(phyl_boot), phyl_boot_bl)
colnames(tax_tree)[[6]]<-"baseline"
tax_tree <- as.data.frame(t(tax_tree))

colnames(tax_tree) <- c("mean", "sd")
tax_tree["type"] <- "taxonomy tree"
print("tree done")


df_a <- rbind(simple_lm, phylo_groups)

df_a[df_a$feature_set == "sequence features",]$feature_set <- "CG comp."
df_a[df_a$feature_set == "CpG criteria",]$feature_set <- "CGI freq."
df_a$feature_set<-factor(df_a$feature_set, 
                       levels=c("no features", "CG comp.", "CGI freq.",
                                "1mers", "2mers", "3mers"))

ggplot(df_a, aes(x = feature_set, y=mean_R2.adjusted, fill = mean_AIC))+
  geom_bar(stat="identity", position = "dodge") + rotate_labels()+
  geom_errorbar(aes(ymin=mean_R2.adjusted-sd_R2.adjusted, ymax=mean_R2.adjusted+sd_R2.adjusted), 
                width=.2,position=position_dodge(.9)) + facet_wrap(~type) +
  scale_fill_gradient(low="#56B1F7", high="#132B43") + 
  geom_hline(yintercept = df_a[df_a$feature_set=="3mers" & df_a$type=="no prior","mean_R2.adjusted"], linetype="dashed")
#+
  scale_fill_manual(values = group_colors)+labs(x="",  y="Variance explained", fill= "prior information")

##merging the results and saving the data (terrible!)
r2_sq<-t(rbind(rbind(tax_tree[1,], phylo_groups[1,]), simple_lm[1,]))
r2_sq<-cbind(r2_sq, as.data.frame(c("sequence features", "CpG criteria", "1mers", "2mers", "3mers", "no features")))
colnames(r2_sq)<-c("taxonomy tree", "phyl. groups", "no prior", "feature_set")
r2_m<-melt(r2_sq, id.vars="feature_set")
colnames(r2_m)[[3]]<-"mean_r2"

r2_sd_sq<-t(rbind(rbind(tax_tree[2,], phylo_groups[2,]), simple_lm[2,]))
r2_sd_sq<-cbind(r2_sd_sq, as.data.frame(c("sequence features", "CpG criteria", "1mers", "2mers", "3mers", "no features")))
colnames(r2_sd_sq)<-c("taxonomy tree", "phyl. groups", "no prior", "feature_set")
r2_sd_m<-melt(r2_sd_sq, id.vars="feature_set")
colnames(r2_sd_m)[[3]]<-"r2_sd"

r2<-full_join(r2_sd_m, r2_m)

r2$variable<-factor(r2_m$variable, levels=names(group_colors))

r2$feature_set<-factor(r2$feature_set, 
              levels=c("no features", "sequence features", "CpG criteria",
                                  "1mers", "2mers", "3mers"))

my_wt(r2, paste0(analysis_dir, "/02_predict_meth/02.6_predictability/features_", target,"_", N_boot, "_filtered.csv"))
r2<-read.csv(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/features_", target,"_", 
                                                                N_boot, "_filtered.csv"), sep="\t")

r2$feature_set <- as.character(r2$feature_set)
r2[r2$feature_set == "sequence features",]$feature_set <- "CG comp."
r2[r2$feature_set == "CpG criteria",]$feature_set <- "CGI freq."
r2$feature_set<-factor(r2$feature_set, 
                       levels=c("no features", "CG comp.", "CGI freq.",
                                "1mers", "2mers", "3mers"))



pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/features_", target,"_", N_boot, "_filtered.pdf"), 
    height = 4, width = 6)
ggplot(r2, aes(x = feature_set, y=mean_r2, fill = variable))+
  geom_bar(stat="identity", position = "dodge") + rotate_labels()+
  geom_errorbar(aes(ymin=mean_r2-r2_sd, ymax=mean_r2+r2_sd), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values = group_colors)+labs(x="",  y="Variance explained", fill= "prior information")
dev.off()

r2$variable <- factor(r2$variable, levels = c("taxonomy tree", "phyl. groups", "no prior"))

pdf(paste0(analysis_dir, "/02_predict_meth/02.6_predictability/features_grid_", target,"_", 
           N_boot, "filtered.pdf"), height = 3.5, width = 6)
ggplot(r2, aes(x = feature_set, y=mean_r2, fill = variable))+ facet_wrap(~variable, nrow = 1)+
  geom_bar(stat="identity", position = "dodge") + rotate_labels()+
  geom_errorbar(aes(ymin=mean_r2-r2_sd, ymax=mean_r2+r2_sd), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values = group_colors)+labs(x="",  y="Variance explained", fill= "prior information")+
  theme(legend.position = "None")+geom_hline(yintercept = r2[r2$feature_set=="3mers" & r2$variable=="no prior","mean_r2"], linetype="dashed")
dev.off()


##
