source(file.path("../../src/00.0_init.R"))

library(kebabs)
library(caret)
library(ROCR)

annot_full <- fread(file.path(analysis_dir, "validation", "03_WGBS", "WGBS_prediction_selection.tsv"))

head(annot_full)

path_to_results <- file.path(analysis_dir, "validation", "03_WGBS", "03.4_prediction")

setwd(path_to_results)

dir.create("summary")

models = list.files(pattern = "methPred", recursive = T) 

roc_res_full = data.table()
for(model_path in models){
    
model_name <- strsplit(model_path, "/")[[1]][4]
    
tryCatch({simpleCache(cacheName=gsub(".RData", "", model_name), instruction={ train_test(x_train=split_ds$x_train,
            x_test = split_ds$x_test, y_train = split_ds$y_train, y_test = split_ds$y_test,
                                                        ifRand='noRand', k=kmer, runid = 0)},
            cacheDir=gsub(model_name, "", model_path ), assignToVariable="res", recreate=FALSE)
roc_res=res$roc_dt   
roc_res$species <- strsplit(model_path, "/")[[1]][1]
roc_res$model <-model_name
roc_res_full <- rbind(roc_res_full, roc_res, fill = TRUE)}, error = function(cond){
    print(cond)
    print(model_name)
})
}

roc_res_full[, seedid:=gsub(".RData", "",
                            strsplit(model, "_")[[1]][length(strsplit(model, "_")[[1]])]),
                            by = row.names(roc_res_full)]
    
roc_res_full[, unique_run := paste0(run, "_", seedid), row.names(roc_res_full)]

setDT(roc_res_full)
my_wt(roc_res_full, file.path("summary", "model_results_combined.tsv"))

roc_res_full <- fread(file.path("summary", "model_results_combined.tsv"))

head(roc_res_full)

roc_res_full <-roc_res_full[seedid != "502"]

unique(roc_res_full$species)

auc_res <- unique(roc_res_full[ifRand == "noRand", c("auc", "ifRand", "run","seedid", "species")]) %>% 
group_by(species, ifRand) %>% summarise(mean_auc = mean(auc), auc_sd = sd(auc))

head(auc_res)

tail(annot_full)

setDT(auc_res)

auc_res[species == "Danio_rerio_GSE134055",species:="Danio_rerio",]

auc_res$species<- gsub("_", " ", auc_res$species)

annot_full[Species == "Danio_rerio_GSE134055",Species:="Danio_rerio",]

annot_full$Species<- gsub("_", " ", annot_full$Species)

auc_res <- left_join(auc_res, annot_full, by = c("species" = "Species"))

head(auc_res)

setDT(auc_res)

auc_rrbs <- fread(file.path(analysis_dir, "05_predict_meth", "05.1_within_species", "summary", "all_aucs.csv"))

head(auc_rrbs)

frogs <- unique(stats_annot[grep("frog", stats_annot$English), c("species", "English", "ncbi_name")])

auc_frogs <- auc_rrbs %>% filter(species %in% frogs$species) %>% group_by(color_class) %>% 
                summarize(AUC = mean(AUC))
auc_frogs$species <- "ALL FROGS"

auc_rrbs <- rbind(auc_rrbs, auc_frogs, fill = TRUE)

auc_res <- left_join(auc_res, auc_rrbs[, c("species", "AUC")], by = c("Species_match" ="species"))

head(auc_res)

k_res <- unique(roc_res_full[ifRand == "noRand", c("species", "k", "seedid")]) %>%
group_by(species) %>% summarise(K = paste(k, collapse = ", "))

k_res

k_frogs <- auc_rrbs %>% filter(species %in% frogs$species) %>% 
                group_by(color_class) %>% 
                summarize(k = median(k))
k_frogs

tail(auc_rrbs)

auc_rrbs[species == "ALL FROGS"]

auc_rrbs[species == "ALL FROGS"]$k <- k_frogs$k

k_res

k_res<- left_join(k_res, annot_full, by = c("species" = "Species"))

k_res <- left_join(k_res, auc_rrbs[, c("species", "k")], 
                   by = c("Species_match" ="species"))


setDT(k_res)

k_res[,x:=0.6,]
k_res[,y:=0.05,]

setDT(auc_res)

auc_res[,x:=0.6,]
auc_res[,y:=0.1,]

auc_res

k_res

roc_res_full[species == "Danio_rerio_GSE134055",species:="Danio_rerio",]

roc_res_full$species<- gsub("_", " ", roc_res_full$species)

unique(roc_res_full$species)

WGBS_species_order <- c('Branchiostoma lanceolatum','Danio rerio','Xenopus laevis','Chelydra serpentina', 'Gallus gallus', 'Phascolarctos cinereus', 'Mus musculus','Bos taurus')

roc_res_full$species <- factor(roc_res_full$species, levels = WGBS_species_order)

head(k_res)

options(repr.plot.width = 20, repr.plot.height = 8)
ggplot(roc_res_full, aes(x = fdr, y = tpr, color = ifRand)) + 
        geom_line(aes(group=unique_run, alpha=ifRand)) +  facet_wrap(.~species, ncol = 7) + 
    geom_text(data=auc_res,aes(x=x,y=y,
                               label=paste0("ROC-AUC = ",round(mean_auc,2), "±", round(auc_sd, 2), 
                                            " (",round(AUC,2), ")")), color = "black", size = 3)+
geom_text(data=k_res,aes(x=x,y=y,
                               label=paste0("k-mer length = {",K, "} (", k,")")), color = "black", size = 3)+
 scale_color_manual(values=c("rand"="grey","noRand"="blue"))+
    scale_alpha_manual(values=c("rand"=0.2,"noRand"=1)) + 
    labs(x = "False discovery rate", y = "True positive rate") +
    theme(text = element_text(size = 15)) + 
facet_wrap(.~factor(species, levels = WGBS_species_order), ncol = 4)
ggsave("summary/WGBS_auc.pdf", width = 10, height = 8)

options(repr.plot.width = 8, repr.plot.height = 10)
ggplot(roc_res_full, aes(x = fdr, y = tpr, color = ifRand)) + 
        geom_line(aes(group=unique_run, alpha=ifRand)) +  facet_wrap(.~species, ncol = 7) + 
    geom_text(data=auc_res,aes(x=x,y=y,
                               label=paste0("ROC-AUC = ",round(mean_auc,2), "±", round(auc_sd, 2), 
                                            " (",round(AUC,2), ")")), color = "black", size = 3)+
geom_text(data=k_res,aes(x=x,y=y,
                               label=paste0("k-mer length = {",K, "} (", k,")")), color = "black", size = 3)+
 scale_color_manual(values=c("rand"="grey","noRand"="blue"))+
    scale_alpha_manual(values=c("rand"=0.2,"noRand"=1)) + 
    labs(x = "False discovery rate", y = "True positive rate") +
    theme(text = element_text(size = 15)) + 
facet_wrap(.~factor(species, levels = WGBS_species_order), ncol = 2)
ggsave("summary/WGBS_auc_vert.pdf", width = 8, height = 10)


