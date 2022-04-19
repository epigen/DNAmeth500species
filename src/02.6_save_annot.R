source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

column_annotation <- read.csv(file.path(Sys.getenv("CODEBASE"), "DNAmeth500species/meta","suppl_annotatons_where_to.csv"), na.strings = "", row.names = 1, sep=";")
column_annotation <- column_annotation[!is.na(column_annotation$New.name),]

dir.create(file.path(analysis_dir, "supplementary_tables"))
##making stats and checking for unumbigouicity:
##sanity checks between scientific/english names and abbreviatons
stats_annot %>% 
  group_by(scientific_name) %>%
  summarize(n = length(unique(species))) %>%
  filter(n > 1)

stats_annot %>% 
  group_by(species) %>%
  summarize(n = length(unique(scientific_name))) %>%
  filter(n > 1)

stats_annot %>% 
  group_by(species) %>%
  summarize(n = length(unique(English))) %>%
  filter(n > 1)

##fixing the english names for the future:
english_fix_typos = c("Blaubauchracke" = "Blue-bellied roller", 
                      "Common Octopus" = "Common octopus", 
                      "Ostrich"  = "Common ostrich",
                      "Wolf" = "Polarwolf", 
                      "Spot-Fin porcupinefish" = "Spot-fin porcupinefish",
                      "stone marten" = "Stone marten", 
                      "Monterey stalked tunicate" = "Stalked tunicate",
                      "Yellowfin snook" = "YellowFin snook")

for (i in names(english_fix_typos)){
  stats_annot[stats_annot$English == i,]$English <- english_fix_typos[i]
}

unique(stats_annot$English[grepl("[^a-zA-Z0-9 -]", stats_annot$English)])## check for weird symbols

stats_annot[species == "CHD"]$English <- "Charadriidae"
stats_annot[species == "KAN"]$English <- "Kangaroo"
stats_annot[species == "WEA"]$English <- "Weasel"

##some stats per abbreviation:
count_samp <- stats_annot %>% 
  group_by(species, ncbi_name) %>% 
  summarize(N_samples = n(), N_individuals = length(unique(replicate)), tissue_types = paste(sort(unique(Tissue)), collapse = ", "))

count_samp_repl <- stats_annot %>% 
  group_by(species, ncbi_name, replicate) %>% 
  summarize(N_samples = n(), tissue_types = paste(sort(unique(Tissue)), collapse = ", "))

##ordering phylogenetically (pre-done order in 00_init)
stats_annot$species <- factor(stats_annot$species, levels = sp_df$species)

stats_annot <- stats_annot[order(stats_annot$species), ]

##recalculating conversion rate
#stats_annot$conversionRate <- stats_annot$conversionRate/100
stats_annot$k1_unmeth <- 100*(1-stats_annot$k1_unmeth) ##ideal: 100%
stats_annot$k3_meth <- 100*(1-stats_annot$k3_meth) ##ideal: 0

stats_annot$mean_overlap_perc <- 100*(stats_annot$mean_overlap/stats_annot$coveredCpGs)

##now saving the tables:

#####on species level

names_species <- as.character(column_annotation[column_annotation$level== "species",]$column_list)

names_species_new <- as.character(column_annotation[column_annotation$level== "species",]$New.name)
df_species <- unique(stats_annot[,..names_species])

##adding counts
df_species <- left_join(df_species, count_samp)
colnames(df_species) <- c(names_species_new, "Number of samples", "Number of individuals", "Tissue types")


#order: matching key, NCBI annotation, per species statistics
colnames_order_species <- c("Common name", "NCBI scientific name","Taxonomic group", "NCBI ID", "NCBI order", "NCBI class", "NCBI group", "Species abbreviation",  "Number of individuals","Number of samples", "Consensus reference fragments", "Tissue types")

##add tissue types
df_species <- df_species[,..colnames_order_species]

write.table(df_species, file.path(analysis_dir, "supplementary_tables","ST2_species.tsv"),
            quote = F, sep = "\t", row.names = F)

###locations: adding the information about the origin of each sample
locations <- read.csv(file.path(meta_dir, "locations.csv"))

locations$species <- sapply(locations$Abbreviation, function(x) strsplit(as.character(x), "_")[[1]][1] )
locations$replicate <- sapply(locations$Abbreviation, function(x) 
  as.integer(strsplit(as.character(x), "_")[[1]][2] ))
for (i in names(english_fix_typos)){
  locations[locations$English == i,]$English <- english_fix_typos[i]
}
                           
locations$Institution <- as.character(locations$Institution)  
setDT(locations)
locations[Institution == "fiwi", Institution:="FIWI",]
locations[Institution == "MUW", Institution:="MedUni",]
unique(locations$Institution)
                              
#####on individual( = replica) level
names_individuals <- c("species", "English", "ncbi_name", "color_class", as.character(column_annotation[column_annotation$level== "individual",]$column_list))

names_individuals_new <- c("Species abbreviation", "Common name", "NCBI scientific name", "Taxonomic group",as.character(column_annotation[column_annotation$level== "individual",]$New.name))

df_individuals <- unique(stats_annot[,..names_individuals])

df_individuals<- left_join(df_individuals, unique(locations[, c("species","replicate",  "Institution", "Location")]),
          by = c("species","replicate"))

##adding counts
df_individuals <- left_join(df_individuals, count_samp_repl)
                              
colnames(df_individuals) <- c(names_individuals_new, "Institution", "Sample source", "Number of samples", "Tissue types")


#order: matching key, replicas and count per, info, info collection
colnames_order_individuals <- c("Common name", "NCBI scientific name", "Taxonomic group","Individual", "Species abbreviation","Number of samples", "Sex","Age", "Pathological diagnosis", "Preservation", "Tissue types", "Dissection date",   "Institution", "Sample source") # Check institution spelling
                            
df_individuals <- df_individuals[,..colnames_order_individuals]

head(df_individuals)
write.table(df_individuals, file.path(analysis_dir, "supplementary_tables","ST3_individuals.tsv"),
            quote = F, sep = "\t", row.names = F, na = "N/A")

#####on sample level
names_samples <- c("species", "English", "ncbi_name", "color_class","replicate", "mean_overlap_perc",             as.character(column_annotation[column_annotation$level== "sample",]$column_list))
                              
names_samples_new <- c("Species abbreviation", "Common name", "NCBI scientific name", "Taxonomic group", "Individual", "Mean CpG overlap (%)", as.character(column_annotation[column_annotation$level== "sample",]$New.name))

df_samples <- unique(stats_annot[,..names_samples])
                              
df_samples$lowest_common_blastS1_num <- as.numeric(df_samples$lowest_common_blastS1) 

df_samples$lowest_common_blastS2_num <- as.numeric(df_samples$lowest_common_blastS2)  

df_samples[lowest_common_blastS2_num == max_lowest_common, blast_genus :=lowest_common_blastS2, by = row.names(df_samples[lowest_common_blastS2_num == max_lowest_common])]
                              
df_samples[lowest_common_blastS2_num == max_lowest_common, blast_hit :=blast_species2, by = row.names(df_samples[lowest_common_blastS2_num == max_lowest_common])]

df_samples[lowest_common_blastS1_num == max_lowest_common, blast_genus :=lowest_common_blastS1, by = row.names(df_samples[lowest_common_blastS1_num == max_lowest_common])]
                              
df_samples[lowest_common_blastS1_num == max_lowest_common, blast_hit :=blast_species1, by = row.names(df_samples[lowest_common_blastS1_num == max_lowest_common])]

df_samples$mean_overlap_perc <- round(df_samples$mean_overlap_perc, 2)
df_samples$mapping_efficiency <- round(df_samples$mapping_efficiency, 2)   
                              
df_samples$CpG_meth <- round(df_samples$CpG_meth, 2)
df_samples$CpG_meth_cphph<- round(df_samples$CpG_meth_cphph, 2)
                              
df_samples$conversionRate <- round(df_samples$conversionRate, 2) 
df_samples$k1_unmeth <- round(df_samples$k1_unmeth,2)
df_samples$k3_meth <- round(df_samples$k3_meth,2)

df_samples$PDR_mean <- round(df_samples$PDR_mean,2)
df_samples$cont_rat <- round(df_samples$cont_rat,2)                              
colnames(df_samples) <- c(names_samples_new, 'lowest_common_blastS1_num', 'lowest_common_blastS2_num', "BLAST hit", "lowest common genus with BLAST hit")


#order: matching key/id, blast, dedRef stats, methylation stats, cont stats,experimental stats
colnames_order_samples_converted <- c("Common name", "NCBI scientific name",  "Taxonomic group", "Tissue", "Individual", "Sample ID", "Species abbreviation", "Total reads", "Mapped reads", "Alignment rate (%)", "Average CpG methylation (%)", "Covered CpGs", "Mean CpG overlap","Mean CpG overlap (%)",'Covered CpHs', 'Average CpH methylation (%)', 'Average PDR', 'PDR sites',"Conversion rate (%, non-CpGs)","Conversion rate (%, unmethylated spike-in)","Conversion rate (%, methylated spike-in)","CGG reads (%)",  "TGG reads (%)","CGA reads (%)", "TGA reads (%)", "Adenines (%)", "Thymines (%)","Cytosines (%)", "Guanines (%)", "Ns(%)", "Prefragmentation (%)","PCR enrichment cycles","Contamination rate","BLAST hit", "lowest common genus with BLAST hit","Species confirmation")
                        
                              
 df_samples_conv <-  df_samples[`Converion status` == "converted",..colnames_order_samples_converted]                            
                              

write.table(df_samples_conv, file.path(analysis_dir, "supplementary_tables","ST1_samples.tsv"),
            quote = F, sep = "\t", row.names = F, na = "N/A")
df_samples_conv[Tissue %in% core, Tissue_vis:= Tissue, by = `Sample ID`]
df_samples_conv[!Tissue %in% core, Tissue_vis:= "Other", by = `Sample ID`]
  write.table(df_samples_conv, file.path(analysis_dir, "supplementary_tables","ST1_samples_forwebsite.tsv"),
            quote = F, sep = "\t", row.names = F, na = "N/A")                            
### unconv
df_samples_unconv <-  df_samples[`Converion status` == "unconverted",..colnames_order_samples_converted]                            
                              

write.table(df_samples_unconv, file.path(analysis_dir, "supplementary_tables","ST4_unconv_samples.tsv"),
            quote = F, sep = "\t", row.names = F, na = "N/A") 
                              
                              

#####Saving stats, required for the paper:
######Introduction:

cat(paste0("DNA methylation profiles: ", NROW(stats_annot[stats_annot$conversion_type == "converted",])),file = file.path(analysis_dir, "01_basicStats", "count_summary.txt"), sep = "\n")
                              
cat(paste0("taxonomic groups (“species”):  ", length(unique(stats_annot$species))), append = TRUE, file = file.path(analysis_dir, "01_basicStats", "count_summary.txt"), sep = "\n")

cat(paste0("vertebrates:  ", length(unique(stats_annot[stats_annot$ncbi_group == "Vertebrata",]$species))), append = TRUE, file = file.path(analysis_dir, "01_basicStats", "count_summary.txt"), sep = "\n")
                              
cat(paste0("Invertebrates:  ", length(unique(stats_annot[stats_annot$ncbi_group == "Invertebrata",]$species))), append = TRUE, file = file.path(analysis_dir, "01_basicStats", "count_summary.txt"), sep = "\n")


##methods, sample collection part:
##for each group print out N of samples, species and orders
                              
for (class_id in levels(stats_annot$color_class)){
  cat(paste0(class_id, ": ", NROW(stats_annot[stats_annot$color_class==class_id, ]), " samples, ", 
               NROW(stats_annot[conversion_type == "converted" & color_class==class_id, ]), " conv. samples, ",
              length(unique(stats_annot[stats_annot$color_class==class_id, ]$species)), " species, ",
              length(unique(stats_annot[stats_annot$color_class==class_id, ]$ncbi_order)), " orders" ), append = TRUE, file = file.path(analysis_dir, "01_basicStats", "count_summary.txt"), sep = "\n")
}
length(unique(stats_annot$Tissue))
##stats which locations things come from

location_stats <- locations %>% 
  group_by(Institution, Location) %>%
  summarise(N_species = length(unique(species)), N_samples = n())

stats_annot[Tissue %in% core, Tissue_label := Tissue, by = Sample_Name]
stats_annot[!Tissue %in% core, Tissue_label := "other", by = Sample_Name]
                              
tissue_stat <- stats_annot %>% filter(conversion_type == "converted") %>%
                group_by(color_class, Tissue_label) %>% 
                summarise(n = n()) %>% ungroup() %>% pivot_wider(names_from = Tissue_label, values_from = n)
                    

write.csv(location_stats, file.path(analysis_dir, "01_basicStats","location_stats.csv") )
