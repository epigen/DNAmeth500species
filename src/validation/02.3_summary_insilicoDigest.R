## ---------------------------
##
## Insilico digestion analysis
##
## Summary of the simulation of the RRBS prototcol insilicos and evaluation of the percentage of CpGs covered in regulatory elements of interest
##
## Authors: Daria Romanovskaia
##
## Date Created: 2021-10-28
##
##
## ---------------------------

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

library(Biostrings)

read_df <- function(file_path){
print(file_path)
df <- read.table(file_path, sep = "\t")
df <- df[, c(-3)]
colnames(df) <- c('CpGs','CpGs_in_RRBS')
df["type"] <- row.names(df)
df["genome_id"] <- strsplit(file_path, "/")[[1]][10]
return(df)
}

output_dir <- file.path(analysis_dir, "validation", "02_insilico_digest")

annot_df <- fread(file.path(output_dir, "genomes_to_run.tsv"))
colnames(annot_df) <- c("ucsc_db", "ucsc_species", "species", "scientific_name", "class")

annot_df[class== "Ascidiacea",class := "Invertebrata",]
annot_df[class== "Malacostraca",class := "Invertebrata",]
annot_df[class== "Echinoidea",class := "Invertebrata",]
annot_df[class== "Hyperoartia",class := "Jawless_vertebrate",]
head(annot_df)

annot_df[ucsc_db=="petMar3"]

if(!file.exists(file.path(output_dir, "summary_df.tsv"))){
    file_list = list.files(output_dir, pattern="*/*_overlap_stats_final.tsv", recursive = T, full.names = T)
    df_list <- sapply(file_list, read_df, simplify = FALSE)

    summary_df <- rbindlist(df_list)
    setDT(summary_df)

    summary_df[, ratio:=summary_df$CpGs_in_RRBS/summary_df$CpGs,]
    
    summary_df <- inner_join(summary_df, annot_df[, c("ucsc_db", "ucsc_species", "class")], by = c("genome_id" = "ucsc_db"))
    summary_df <- filter(summary_df, genome_id!="bosTau9_noconcat")
    
    ##cleaning up the dataframe
    summary_df[summary_df$class == "Marsupiala", class:= "Marsupialia",]
    summary_df[summary_df$type == "rmsk", type:= "repeats",]
    summary_df[summary_df$type == "simpleRepeat", type:= "repeats",]
    summary_df[summary_df$type == "cpgIslandExt", type:= "CpG Islands",]
    summary_df[summary_df$type == "promoters1000_500", type:= "promoters",]
    
    ##sorting the dataframe
    summary_df$class <- factor(summary_df$class, levels = names(class_colors))
    summary_df <- summary_df[order(summary_df$class),]
    
    summary_df$type <- factor(summary_df$type, levels = c("total", "transcripts", "promoters",
                                                      "repeats", "CpG Islands"))
    
    summary_df[, group:=class_short[class],]
    summary_df$group <- factor(summary_df$group, levels = class_short)
    my_wt(summary_df, file.path(output_dir, "summary_df.tsv")) 
    
    
}else{
    summary_df <- fread(file.path(output_dir, "summary_df.tsv"))
}

summary_df[summary_df$ucsc_species == "Lamprey",group:="Inv."]
summary_df[summary_df$ucsc_species == "Lamprey",class:="Invertebrata"]

summary_df$class <- factor(summary_df$class, levels = names(class_colors))
summary_df <- summary_df[order(summary_df$class),]
    
summary_df$type <- factor(summary_df$type, levels = c("total", "transcripts", "promoters",  "repeats", "CpG Islands"))
    
summary_df[, group:=class_short[class],]

summary_df$group <- factor(summary_df$group, levels = class_short)


ggplot(summary_df, aes(x = group, y = ratio, color = class)) + geom_boxplot(outlier.shape = NA) + 
       geom_jitter( alpha = 0.5, shape = 16) +
        geom_point(data = summary_df[summary_df$ucsc_species == "Lamprey",], aes(x = group, y = ratio), color = class_colors[[2]]) + 
        facet_wrap(~type, nrow = 1) + 
        theme_bw() + ylim(c(0,1)) + 
        theme( legend.position = "None") +
        scale_color_manual(values = class_colors) + 
        geom_hline(yintercept = 0.25, linetype = 2, color = "#bdbdbd") + 
        geom_hline(yintercept = 0.5, linetype = 2, color = "#bdbdbd") +
        geom_hline(yintercept = 0.75, linetype = 2, color = "#bdbdbd")+stat_summary(fun.data = give.n,fun.args = c(y=0.1),geom = "text", angle = 90, hjust =1)
ggsave(file.path(output_dir, "summary_boxplots.pdf"), width = 12, height = 3 )

ggplot(summary_df[class!="Jawless_vertebrate"], aes(x = group, y = ratio, color = class)) + geom_boxplot(outlier.shape = NA) + 
        geom_jitter(data = summary_df, aes(x = group, y = ratio, color = class), alpha = 0.5) + 
        facet_wrap(~type, nrow = 1) + 
        theme_bw() + ylim(c(0,1)) + 
        theme( legend.position = "None") +
        scale_color_manual(values = class_colors) + 
        geom_hline(yintercept = 0.25, linetype = 2, color = "#bdbdbd") + 
        geom_hline(yintercept = 0.5, linetype = 2, color = "#bdbdbd") +
        geom_hline(yintercept = 0.75, linetype = 2, color = "#bdbdbd")  + 
        geom_label_repel(data = summary_df[ucsc_species %in% c("Human", "Mouse", "Zebrafish")], 
                        aes(x = group, y = ratio, label = substr(ucsc_species, 1, 1)), color = "black", alpha = 0.8)
ggsave(file.path(output_dir, "summary_boxplots_names.pdf"), width = 12, height = 3 )


