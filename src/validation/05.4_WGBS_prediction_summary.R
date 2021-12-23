source(file.path("../../src/00.0_init.R"))

annot_full <- fread("../meta/WGBS_prediction_selection.tsv")

head(annot_full)

path_to_results <- "../validation/WGBS_public/prediction/"

kebabs_dirs <- list.dirs(path = "../validation/WGBS_public/prediction/Gallus_gallus/galGal5/", recursive = F)

kebabs_dirs

full_full_df <- data.table()
for(i in seq(1:NROW(annot_full))){
    species <- annot_full$Species[[i]]
    assembly <- annot_full$assembly[[i]]
    full_df <- data.table()
    kebabs_dirs <- list.dirs(path = paste0(path_to_results, "/", species, "/", assembly), recursive = F)
    
    for(i in seq_along(kebabs_dirs)){
        print(kebabs_dirs[i])
    df <- fread(paste0(kebabs_dirs[i], "/", assembly, "roc_res.csv"))
        
    df$seed_id <- strsplit(strsplit(kebabs_dirs[i], "/")[[1]][8], "_")[[1]][3]
    print(df$seed_id[[1]])
    full_df <- rbind(full_df, df)
   
    }
    full_df$species <- gsub("_", " ", species)
    full_df$assembly <- assembly
    full_full_df <- rbind(full_full_df, full_df)
}

                  

head(full_full_df)

full_full_df[is.na(seed_id)]$seed_id <- "1234"

full_full_df[, unique_run:=paste0(run, "_", seed_id),]

head(full_full_df)

full_full_df[species == "Danio rerio GSE134055",species:="Danio rerio",]

rrbs_auc <- fread("/nobackup/lab_bock/projects/DNAmeth500species/results_analysis_v1/02_predict_meth/02.1_within_species/summary/all_aucs.csv")
head(rrbs_auc)

head(annot_full[,-c("Species")])

full_full_df <- left_join(full_full_df, annot_full[,-c("Species")])

auc_res=full_full_df[,list(auc=mean(auc),auc_sd = sd(auc)),by=list(ifRand, species, Species_match)]
auc_res[,x:=0.7,]
auc_res[,y:=ifelse(ifRand=="rand", 0.08, 0.13),]

auc_res

auc_res <- left_join(auc_res, rrbs_auc[, c("species", "AUC")], by = c("Species_match" = "species"))

head(auc_res)

unique(full_full_df$species)

full_full_df$species <- factor(full_full_df$species, levels = c('Branchiostoma lanceolatum','Danio rerio', 'Gallus gallus', 'Phascolarctos cinereus','Mus musculus','Bos taurus'  ))

head(full_full_df$species)

auc_res$species <- factor(auc_res$species, levels = c('Branchiostoma lanceolatum','Danio rerio', 'Gallus gallus', 'Phascolarctos cinereus','Mus musculus','Bos taurus'  ))

options(repr.plot.width = 15, repr.plot.height = 4)
ggplot(full_full_df, aes(x = fdr, y = tpr, color = ifRand)) + geom_line(aes(group=unique_run, alpha=ifRand)) + 
    geom_text(data=auc_res,aes(x=x,y=y,label=paste0("auc=",round(auc,2), "±", round(auc_sd, 2)))) +
     geom_text(data=auc_res[ifRand == "noRand"],aes(x=x,y=y+0.05,label=paste0("auc RRBS =",round(AUC,2))), color = "red") + 
    scale_color_manual(values=c("rand"="grey","noRand"="blue", "RRBS" = "red"))+
    scale_alpha_manual(values=c("rand"=0.2,"noRand"=1)) + 
    labs(x = "False discovery rate", y = "True positive rate") +
    theme(text = element_text(size = 15)) + facet_wrap(.~species, ncol = 6)
ggsave("../validation/WGBS_public/prediction/WGBS_auc.pdf", width = 20, height = 4)

write.table(unique(full_full_df[ifRand == "noRand", c("species", "seed_id", "k", "auc")]), "../validation/WGBS_public/prediction/WGBS_auc_and_k.csv", sep = ";", quote = F, row.names = F)
