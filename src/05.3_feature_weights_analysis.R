source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(pheatmap)
library(ggseqlogo)
library(ComplexHeatmap)

wd=file.path(analysis_dir,"05_predict_meth", "05.3_feature_weight_analysis")
dir.create(wd)
setwd(wd)

#load the three-mer weights
if(!file.exists("all_kmerWeights_3.csv")){
    data_files = list.files("kmer3", pattern = "methPred_noRand_uc.RData", recursive = TRUE, full.names = TRUE)
    all_kmerWeights=data.table()
    for (file in data_files){
        species=strsplit( file,"/")[[1]][2]
        print(species)
        load(file)
        weights=ret$feature_weights
        setnames(weights,c("variable","value"),c("kmer",species))
        if (nrow(all_kmerWeights)==0){all_kmerWeights=weights}else{
    all_kmerWeights=merge(all_kmerWeights,weights,by="kmer",all=TRUE)
        }
    }
    head(all_kmerWeights)
    my_wt(all_kmerWeights, "all_kmerWeights_3.csv")
}

all_kmer3<-read.csv("all_kmerWeights_3.csv", row.names = 1, sep = "\t")
#all_kmer3<-as.data.frame(t(all_kmer3))
NROW(all_kmer3)

sp_order <- sp_df[species %in% colnames(all_kmer3),]$species
all_kmer3 <- all_kmer3[,sp_order]

sp_df <- as.data.frame(sp_df)
row.names(sp_df) <- sp_df$species
ha <- columnAnnotation(class=sp_df[colnames(all_kmer3),]$color_class, 
                       col=list(class=class_colors), show_legend = F)
##normalizing per column (a.k.a. species)
pdf("full_weights_pretty_3mers.pdf", height = 10, width = 10)
Heatmap(scale(all_kmer3), cluster_columns = F, show_column_names = F, top_annotation = ha, name = "scaled f.w.")
dev.off()


kmer_mean <- rowMeans(all_kmer3)
kmer_sd<-apply(all_kmer3, 1, sd)
kmer_over_all<-as.data.frame(t(rbind(kmer_sd, kmer_mean)))
kmer_over_all$kmer<-row.names(kmer_over_all)

kmer_over_all$kmer<-factor(kmer_over_all$kmer, levels = kmer_over_all$kmer)

kmer_over_all<-kmer_over_all[order(abs(kmer_over_all$kmer_mean), decreasing = T),]
ggplot(kmer_over_all[c(1:10),], aes(x=kmer, y = kmer_mean))+geom_bar(stat="identity")+ylim(c(-6,6))+
  geom_errorbar(aes(ymin=kmer_mean-kmer_sd, ymax=kmer_mean+kmer_sd),
                width=.2,position=position_dodge(.9)) + rotate_labels()
 

kmer_over_all<-kmer_over_all[order(abs(kmer_over_all$kmer_sd), decreasing = T),]
ggplot(kmer_over_all[c(1:10),], aes(x=kmer, y = kmer_mean))+geom_bar(stat="identity")+ylim(c(-6,6))+
  geom_errorbar(aes(ymin=kmer_mean-kmer_sd, ymax=kmer_mean+kmer_sd), 
                width=.2,position=position_dodge(.9))+rotate_labels()

all_kmer3 <- as.data.frame(t(all_kmer3))
all_kmer3$species <- row.names(all_kmer3)
all_kmer3 <- left_join(all_kmer3, sp_df[, c("species", "color_class")])
setDT(all_kmer3)
all_kmer3[color_class == "Jawless_vertebrate", color_class := "Invertebrata"]

group_mean <- all_kmer3 %>% select(-species) %>% group_by(color_class) %>%
  summarise_all(mean) %>% pivot_longer(-color_class, names_to = "kmer", values_to = "mean")

group_sd <- all_kmer3 %>% select(-species) %>% group_by(color_class) %>%
  summarise_all(sd) %>% pivot_longer(-color_class, names_to = "kmer", values_to = "sd")



group_cv<-merge(group_mean, group_sd)
group_cv$cv<-group_cv$sd/group_cv$mean
group_cv$letter<-sapply(group_cv$kmer, function(x) substr(x, 1, 1))

group_cv$group<-factor(unlist(lapply(group_cv$color_class, function(x) class_short[x])), levels=class_short)

ggplot(group_cv, aes(x=kmer, y=mean, fill = color_class, alpha = 0.7, color = color_class))+ 
  geom_point(shape = 21, aes(size=sd))+
  scale_color_manual(values = class_colors)+scale_fill_manual(values = class_colors)

ggplot(group_cv, aes(x=mean, y=sd, fill = color_class, alpha = 0.7, color = color_class))+ 
  geom_point(shape = 21)+
  scale_color_manual(values = class_colors)+scale_fill_manual(values = class_colors)

ggplot(group_cv)+geom_histogram(aes(x = cv))+facet_wrap(~color_class, ncol = 2, scale="free")

ggplot(group_cv, aes(x = kmer, y = color_class))+geom_point(aes(size=abs(cv)))

pdf("kmer_weights_by_class.pdf",width=8, height=8)
ggplot(group_cv, aes(x = kmer, y = group))+geom_point(aes(size=abs(mean), color = sd))+
  facet_wrap(~letter, scales = "free")+rotate_labels()+scale_color_gradient(high = "grey", low = "#56B1F7")
dev.off()

ggplot(all_kmer3)+geom_density(aes(x=GGG, fill = color_class, alpha = 0.5))+
  scale_fill_manual(values=class_colors)


## Now let draw a logo for all the positive  weighted kmers
get_pwm_vector<-function(df, n){
  df$letter_x<-sapply(df$kmer, function(x) substr(x, n, n))
  df_sum<-df %>% 
    group_by(letter_x) %>%
    summarise(s_p = sum(mean[mean>0]), s_n=sum(mean[mean<0]))
  colnames(df_sum)[[2]]<-paste0("s_p_", n)
  colnames(df_sum)[[3]]<-paste0("s_n_", n)
  return(df_sum)
}
                      
dir.create("logos_3mers")
                      
draw_logo<-function(df, class_name){
  print(class_name)
  l<-lapply(seq(1, 3), get_pwm_vector, df = df[df$color_class==class_name, ])
  full_l<-Reduce(full_join, l)
  setDT(full_l)
  
  logo_pos<-as.matrix(full_l[, c(2,4,6)])
  row.names(logo_pos)<-full_l$letter
  p1<-ggseqlogo(logo_pos)
  
  logo_neg<-as.matrix(full_l[, c(3,5,7)])
  row.names(logo_neg)<-full_l$letter
  weights<-list("neg. weights"=abs(logo_neg), "pos. weights" = logo_pos)
  p<-ggseqlogo(weights, method = "custom", ncol = 2,scales="free")
  
  pdf(paste0("logos_3mers", class_name, ".pdf" ), width = 4, height = 2)
  print(p)
  dev.off()
  
  return(weights)
  }

w<-lapply(unique(group_cv$color_class), draw_logo, df=group_cv)

all_kmer3 <- as.data.frame(all_kmer3)
##running stats test:
w_results<-lapply(colnames(all_kmer3)[-c(65,66)], function(k){ sapply(unique(all_kmer3$color_class),
      function(x) {p<-wilcox.test(all_kmer3[all_kmer3$color_class==x, k], 
          all_kmer3[all_kmer3$color_class!=x, k], paired=F)$p.value; return(setNames(p,x))})})
                      
w_results<-setNames(w_results, colnames(all_kmer3)[c(-65,-66)])
w_results_df<-as.data.frame(w_results)
w_results_df$color_class<-row.names(w_results_df)

write.csv(w_results_df, "significance_feature_weights.csv", quote = F)

w_results_df_m<-melt(w_results_df, id.vars=c("color_class"))
w_results_df_m$color_class<-factor(w_results_df_m$color_class, levels = names(class_colors))
w_results_df_m$val<-log10(w_results_df_m$value)*(-1)

pdf("p_value_kmers.pdf", width = 16, height = 16)
ggplot(w_results_df_m, aes(x=color_class, y = val))+
  geom_bar(stat = "identity")+coord_flip()+
  facet_wrap(~variable, ncol = 8)+labs(y="-log10(p.value)", x="phylogenetic group")
dev.off()

all_kmers_m<-melt(all_kmer3, id.vars=c("species", "color_class"))

kmer_sign_count<-all_kmers_m %>%
  group_by(.dots=c("color_class", "variable")) %>%
  summarize(n_pos=sum(value>0), n_neg=sum(value<0), n=n())
colnames(kmer_sign_count)[[2]]<- "kmer"
##and now by class letshow the biggest and the smallest kmers:
ggplot(w_results_df_m, aes(x=variable, y = val))+
  geom_bar(stat = "identity")+coord_flip()+
  facet_wrap(~color_class, ncol = 4)+labs(y="-log10(p.value)", x="phylogenetic group")

## we need to filter to show top 10 at each
colnames(w_results_df_m)[[2]]<-"kmer"
colnames(w_results_df_m)[[3]]<-"pval"
w_top<-merge(w_results_df_m, group_cv[, c("color_class", "kmer", "mean", "sd")], all.y=FALSE)

w_top<-w_top %>%
  group_by(color_class) %>%
  #selecting top 10 in each category 
  top_n(n = 10, wt = val) %>%
  ungroup() %>%
  arrange(color_class, abs(mean)) %>%
  mutate(order = row_number())

#w_top$neg<-w_top$val*w_top$n_neg/w_top$n
#w_top$pos<-w_top$val*w_top$n_pos/w_top$n
#w_top$zero<-w_top$val*(w_top$n - w_top$n_neg - w_top$n_pos)/w_top$n
#w_top<-melt(w_top, measure.vars=c("neg", "pos", "zero"))


pdf("weights_and_pval_kmers.pdf", width = 16, height = 6)
ggplot(w_top, aes(x=order, y = mean, fill = val))+
  geom_bar(stat = "identity")+
  facet_wrap(~color_class, ncol = 4, scales = "free_y")+
  scale_x_continuous(breaks = w_top$order, labels = w_top$kmer,
    expand = c(0,0)
  )+coord_flip()+
  labs(fill="-log10\n(p.value)", x="kmer", y="mean feature weight")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                width=.2,position=position_dodge(.9))+
  scale_fill_continuous(low = "#9ecae1", high = "#08519c")+
  geom_hline(yintercept = 0, alpha = 0.5)
dev.off()

ggplot(w_results_df_m)+geom_density(aes(x = val, fill = color_class), alpha = 0.5)+
  scale_fill_manual(values = class_colors)


ggplot(w_results_df_m[w_results_df_m$val>2, ])+geom_histogram(aes(x=color_class, fill=color_class), stat="count")+
  rotate_labels()+scale_fill_manual(values = class_colors)

w_top<-w_results_df_m %>%
  group_by(color_class) %>%
  top_n(n = 10, wt = val)

w_results_df_m<-w_results_df_m[order(w_results_df_m$kmer),]

pdf("p_value_top_kmers.pdf", width = 16, height = 16)
ggplot(w_results_df_m[w_results_df_m$kmer %in% w_top$kmer, ], aes(x = color_class, y = kmer))+
  geom_tile(aes(fill = val))+rotate_labels()+labs(x="phylogenetic group", y="", fill="-log10(p)")
dev.off()

