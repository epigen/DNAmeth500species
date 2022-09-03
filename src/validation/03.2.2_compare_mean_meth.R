## ---------------------------
##
## WGBS-based validation
##
## comparing the WGBS and RRBS data
##
## Authors: Daria Romanovskaia
##
## Date Created:  2021-05-20
##
##
## ---------------------------

source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

COV=5

wd=file.path(analysis_dir, "validation", "03_WGBS", "03.2_mean_meth")
setwd(wd)

##collect WGBS
files_list <-list.files( pattern = "*.csv")
files_list <- files_list[!grepl("summary",files_list)]

wgbs_mean = data.table()
for(path in files_list){
    print(path)
    df <- read.csv(path, sep = ";")
    df$source <- gsub(".csv", "",path)
    wgbs_mean <- rbind(df, wgbs_mean, fill = TRUE)
    }


wgbs_mean <- wgbs_mean[, c("path", "mean_ratio", "tissue","replica", "source")]

## clean-up of the tissue types
wgbs_mean[, latin_name:= sapply(wgbs_mean$path, function(x) strsplit(as.character(x), "/")[[1]][5]),]
wgbs_mean[latin_name == "Branchiostoma_lanceolatum_1", latin_name:= "Branchiostoma_lanceolatum",]
wgbs_mean[,tissue:=tolower(tissue),]
wgbs_mean[,tissue:=as.character(sapply(wgbs_mean$tissue, function(x) ifelse(length(grep("brain", x) > 0), "brain", x))),]
wgbs_mean[tissue == "frontal_cortex",tissue:="brain"]
wgbs_mean[tissue == "cerebellum",tissue:="brain"]
wgbs_mean[tissue == "right_hemisphere_cerebral_lobe",tissue:="brain"]
wgbs_mean[tissue == "olfactory_bulb",tissue:="brain"]  
#wgbs_mean[latin_name == "DNAmeth500species", latin_name:="Chelydra_serpentina" ]

wgbs_mean[source=="Chelydra_serpentina", latin_name := "Chelydra_serpentina"]



wgbs_mean[source=="Parus_major_PRJNA574487", 
          latin_name := "Parus_major"]

wgbs_mean[source== "Danio_rerio_GSE134055", replica:=3,]

wgbs_mean[latin_name== "Branchiostoma_lanceolatum"]

wgbs_mean[source == "Chelydra_serpentina"]

wgbs_mean[latin_name=="Danio_rerio",replica:=4]

wgbs_mean[latin_name=="Danio_rerio_gemBS", latin_name:="Danio_rerio",]

wgbs_mean[latin_name=="Branchiostoma_lanceolatum" & source == "Combined_study_GSE141609", replica:=2,]

head(wgbs_mean)

unique(wgbs_mean$latin_name)

## reading in the matches in the RRBS data

match <- fread(file.path(analysis_dir,  "validation", "03_WGBS","WGBS_RRBS_match.tsv"), 
                 fill = TRUE, header = TRUE, sep = "\t")


match[,latin_name := sapply(`scientific name`, function(x) paste(strsplit(x, " ")[[1]], collapse = "_")), ]

match[,latin_name :=sapply(latin_name, function(x) paste0(toupper(substr(x, 1,1)), substr(x, 2,nchar(x)))),]

match[latin_name== "Xenopus_laevis", c("latin_name", "Species")]

match[latin_name== "Parus_major", c("latin_name", "Species")]

match[latin_name== "Xenopus_laevis", Species:="ALL FROGS",]

# adding frogs for Xenopus

frogs <- unique(stats_annot[grep("frog", stats_annot$English), c("species", "English", "ncbi_name")])

frogs$latin_name <- "Xenopus_laevis"

colnames(frogs)[[1]] <- "Species"

wgbs_mean <- inner_join(wgbs_mean, unique(match[, c("latin_name", "Species")]))

wgbs_mean[latin_name== "Xenopus_laevis"]$Species

wgbs_mean[latin_name== "Parus_major"]

my_wt(wgbs_mean, "WGBS_mean_per_sample_summary.tsv")

stats_annot_mean_meth  <- stats_annot %>% filter(conversion_type == "converted") %>% group_by(color_class, species, English) %>% 
                summarize(m_rrbs = mean(CpG_meth), min_rrbs =min(CpG_meth), max_rrbs=max(CpG_meth) ) %>% ungroup()

stats_annot_mean_meth_frogs <- stats_annot %>% filter(species %in% frogs$Species & conversion_type == "converted") %>% group_by(color_class) %>% 
                summarize(m_rrbs = mean(CpG_meth), min_rrbs =min(CpG_meth), max_rrbs=max(CpG_meth) ) %>% ungroup()

head(stats_annot_mean_meth)

head(stats_annot_mean_meth_frogs)

stats_annot_mean_meth_frogs$English <- "Clawed frog" ##the English name of the WGBS match
stats_annot_mean_meth_frogs$species <- "ALL FROGS"

stats_annot_mean_meth <- rbind(stats_annot_mean_meth, stats_annot_mean_meth_frogs)

setDT(stats_annot_mean_meth)

setDT(wgbs_mean)

wgbs_mean[Species == "ALL FROGS"]

wgbs_mean <- inner_join(wgbs_mean, stats_annot_mean_meth, by = c("Species" = "species"))
#ordering by evolutionary classes
wgbs_mean$color_class <- factor(wgbs_mean$color_class, levels = names(class_colors))
wgbs_mean <- wgbs_mean[order(color_class), ]
wgbs_mean$Species <- factor(wgbs_mean$Species, levels = unique(wgbs_mean$Species))
wgbs_mean$English <- factor(wgbs_mean$English, levels = unique(wgbs_mean$English))

wgbs_mean[mean_ratio > 1, mean_ratio:=mean_ratio/100]

ggplot(wgbs_mean, aes(x = English, y = mean_ratio, fill = color_class)) + geom_boxplot(outlier.shape =  21) + theme_bw() + 
theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) + 
scale_fill_manual(values = class_colors) + labs(x = "", y = "mean WGBS meth. level", fill = "class")
ggsave("WGBS_boxplot.pdf", height = 6, width = 8)

ggplot(wgbs_mean, aes(x = English, y = mean_ratio)) + geom_boxplot(outlier.shape = NA, aes(color = color_class)) + theme_bw() + 
geom_jitter(shape = 21, size = 2,aes(fill = tissue)) + scale_color_manual(values = class_colors) + 
theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "", y = "mean WGBS meth. level", color = "class")
ggsave("WGBS_boxplot_tissue.pdf", height = 8, width = 15)

wgbs_mean_per_tissue <- wgbs_mean %>% group_by(Species, latin_name, tissue, English) %>% summarize(m = mean(mean_ratio))

options(repr.plot.width = 20, repr.plot.height = 10)
ggplot(wgbs_mean_per_tissue, aes(x = English, y = tissue,fill = m, size = m)) + geom_point(shape = 21) + 
    theme_bw()+
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(low = "#3288bd", mid = "#ffffbf", high = "#9e0142", midpoint = 0.5) + 
    labs(x = "", y = "", fill = "mean WGBS meth.", size = "mean WGBS meth." )
ggsave("WGBS_mean_per_tissue.pdf", height = 10, width = 10)  

wgbs_mean_per_species <- wgbs_mean %>% filter(tissue!= "placenta") %>% group_by(latin_name, Species) %>% summarise(m_wgbs = mean(mean_ratio), min_wgbs = min(mean_ratio),
                                                                                                                   max_wgbs=max(mean_ratio), n = n())

wgbs_mean_per_species <- inner_join(wgbs_mean_per_species, stats_annot_mean_meth, by = c("Species" = "species"))



options(repr.plot.width = 10, repr.plot.height = 10)
ggplot(wgbs_mean_per_species, aes(x = m_wgbs*100, y = m_rrbs, color = color_class)) + 
geom_abline(slope = 1, alpha = 0.5, linetype = "dashed") +
geom_point(aes(size = n))  +
theme_bw()+
theme(text = element_text(size = 8)) + 
scale_color_manual(values = class_colors) + 
xlim(c(0,100)) + ylim(c(0,100)) + coord_equal() + 
labs(x = "mean WGBS meth.", y = "mean RRBS meth.", color = "class", size = "number of\nsamples(WGBS)")  + 
ggtitle(paste0("corr = ", round(cor(wgbs_mean_per_species$m_rrbs,wgbs_mean_per_species$m_wgbs),2))) + 
 geom_errorbar(aes(xmin=100*min_wgbs, xmax=100*max_wgbs), width=.5,
                 position=position_dodge(.9)) +
geom_errorbar(aes(ymin=min_rrbs, ymax=max_rrbs), width=.5,
                 position=position_dodge(.9))+ geom_text_repel(aes(label = English))
ggsave("WGBS_RRBS.pdf", height = 4, width = 5) 


options(repr.plot.width = 10, repr.plot.height = 10)
ggplot(wgbs_mean_per_species, aes(x = m_wgbs*100, y = m_rrbs, color = color_class)) + 
geom_abline(slope = 1, alpha = 0.5, linetype = "dashed") +
geom_point(aes(size = n))  +
theme_bw()+
theme(text = element_text(size = 8)) + 
scale_color_manual(values = class_colors) + 
xlim(c(0,100)) + ylim(c(0,100)) + coord_equal() + 
labs(x = "mean WGBS meth.", y = "mean RRBS meth.", color = "class", size = "number of\nsamples(WGBS)")  + 
ggtitle(paste0("corr = ", round(cor(wgbs_mean_per_species$m_rrbs,wgbs_mean_per_species$m_wgbs),2))) + 
 geom_errorbar(aes(xmin=100*min_wgbs, xmax=100*max_wgbs), width=.5,
                 position=position_dodge(.9)) +
geom_errorbar(aes(ymin=min_rrbs, ymax=max_rrbs), width=.5,
                 position=position_dodge(.9))#+ geom_text_repel(aes(label = English))
ggsave("WGBS_RRBS_nolabels.pdf", height = 6, width = 6) 

my_wt(wgbs_mean_per_tissue, "WGBS_summary_tissues.csv")

my_wt(wgbs_mean_per_species, "WGBS_RRBS_summary.csv")

getwd()


