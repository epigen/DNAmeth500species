source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

column_annotation <- read.csv(file.path(Sys.getenv("CODEBASE"),
                                        "DNAmeth500species/meta","suppl_annotatons_where_to.csv"), 
                              na.strings = "", row.names = 1, sep=";")
column_annotation <-column_annotation[!is.na(column_annotation$New.name),]

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
                      "YellowFin snook" = "YellowFin snook"
                      )
for (i in names(english_fix_typos)){
  stats_annot[stats_annot$English == i,]$English <- english_fix_typos[i]
}

##some stats per abbreviation:
count_samp <- stats_annot %>% 
  group_by(species, ncbi_name) %>% 
  summarize(N_samples = n(), N_individuals = length(unique(replicate)))

count_samp_repl <- stats_annot %>% 
  group_by(species, ncbi_name, replicate) %>% 
  summarize(N_samples = n())

##ordering phylogenetically (pre-done order in 00_init)
stats_annot$species <- factor(stats_annot$species, levels = sp_df$species)
stats_annot <- stats_annot[order(stats_annot$species), ]

##recalculating conversion rate
stats_annot$conversionRate <- stats_annot$conversionRate/100
stats_annot$k1_unmeth <- 1-stats_annot$k1_unmeth
stats_annot$k3_meth <- 1-stats_annot$k3_meth

##now saving the tables:
#####on species level

names_species <- as.character(column_annotation[column_annotation$level== "species",]$column_list)
names_species_new <- as.character(column_annotation[column_annotation$level== "species",]$New.name)
df_species <- unique(stats_annot[,..names_species])

##adding counts
df_species <- left_join(df_species, count_samp)
colnames(df_species) <- c(names_species_new, "Number of samples", "Number of individuals")


#order: matching key, NCBI annotation, per species statistics
colnames_order_species <- c("Species abbreviation", "Common name", "NCBI scientific name", "NCBI ID",
                    "NCBI order", "NCBI class", "NCBI group", "Practical class", "Scientific name",
                    "Number of individuals","Number of samples", "Deduced reference fragments" )
df_species <- df_species[,colnames_order_species]
write.table(df_species, file.path(analysis_dir, "supplementary_tables","ST2_species.tsv"),
            quote = F, sep = "\t", row.names = F)

###locations: adding the information about the origin of each sample
locations <- read.csv(file.path(Sys.getenv("CODEBASE"),"compEpi/meta/locations.csv"))

locations$species <- sapply(locations$Abbreviation, function(x) strsplit(as.character(x), "_")[[1]][1] )
locations$replicate <- sapply(locations$Abbreviation, function(x) 
  as.integer(strsplit(as.character(x), "_")[[1]][2] ))
for (i in names(english_fix_typos)){
  locations[locations$English == i,]$English <- english_fix_typos[i]
}
#####on individual( = replica) level
names_individuals <- c("species", "English", "ncbi_name", "color_class",
                   as.character(column_annotation[column_annotation$level== "individual",]$column_list))
names_individuals_new <- c("Species abbreviation", "Common name", "NCBI scientific name", "Practical class",
            as.character(column_annotation[column_annotation$level== "individual",]$New.name))

df_individuals <- unique(stats_annot[,..names_individuals])

df_individuals<- left_join(df_individuals, unique(locations[, c("species","replicate",  "Institution", "Location")]),
          by = c("species","replicate"))

##adding counts
df_individuals <- left_join(df_individuals, count_samp_repl)
colnames(df_individuals) <- c(names_individuals_new, "Institution", "Location", "Number of samples")


#order: matching key, replicas and count per, info, info collection
colnames_order_individuals <- c("Species abbreviation", "Replicate", "Common name", "NCBI scientific name", 
                                "Practical class","Number of samples", 
                            "Sex","Age",  "Disease","Preservation",  
                            "Dissection date", "Internal ID",  "Institution", "Location")
                            
df_individuals <- df_individuals[,colnames_order_individuals]


write.table(df_individuals, file.path(analysis_dir, "supplementary_tables","ST3_individuals.tsv"),
            quote = F, sep = "\t", row.names = F, na = "N/A")

#####on sample level
names_samples <- c("species", "English", "ncbi_name", "color_class","replicate",
                       as.character(column_annotation[column_annotation$level== "sample",]$column_list))
names_samples_new <- c("Species abbreviation", "Common name", "NCBI scientific name", "Practical class", "Replicate",
                           as.character(column_annotation[column_annotation$level== "sample",]$New.name))

df_samples <- unique(stats_annot[,..names_samples])
colnames(df_samples) <- names_samples_new


#order: matching key/id, blast, dedRef stats, methylation stats, cont stats,experimental stats
colnames_order_samples <- c("Sample ID", "Species abbreviation", "Replicate", "Tissue","Converion status",
                            "Common name", "NCBI scientific name", "Practical class",
                             
                            
                            "BLAST hit 1",  "BLAST hit 2", "BLAST hit 1 (lowest common rank)", 
                            "BLAST hit 2 (lowest common rank)",  "Best lowest common rank",
                            
                            "Total reads", "Mapped reads", "Mapping efficiency",
                            "Average CpG methylation", "Covered CpGs", "Mean CpG overlap", 
                            "Minimum CpG overlap", "Maximum CpG overlap", 
                            "Minimum CpG overlap (%)", "Maximum CpG overlap (%)",
                            
                            "Conversion rate (non-CpG Cs)","Conversion rate (unmethylated spike-in)",
                            "Conversion rate (methylated spike-in)",
                            
                            "CGG reads (%)",  "TGG reads (%)","CGA reads (%)", "TGA reads (%)", "Prefragmentation (%)",
                            "Adenines (%)", "Thymines (%)","Cytosines (%)", "Guanines (%)", "Ns(%)",
                            
                            "Top contaminating species", "Contamination rate",
                            
                            "DNA isolation date", "PCR enrichment cycles", "PCR enrichment cycles (unconverted)",
                            "Adapter", "Adapter unconverted",                             
                            "Flowcell","Flowcell unconverted",                                
                            "Lane", "Lane unconverted")



df_samples <- df_samples[,..colnames_order_samples]
write.table(df_samples, file.path(analysis_dir, "supplementary_tables","ST1_samples.tsv"),
            quote = F, sep = "\t", row.names = F, na = "N/A")


#####Saving stats, required for the paper:
######Introduction:

print(paste0("DNA methylation profiles: ", NROW(stats_annot[stats_annot$conversion_type == "converted",])))
print(paste0("taxonomic groups (“species”):  ", length(unique(stats_annot$species))))
print(paste0("vertebrates:  ", length(unique(stats_annot[stats_annot$ncbi_group == "Vertebrata",]$species))))
print(paste0("Invertebrates:  ", length(unique(stats_annot[stats_annot$ncbi_group == "Invertebrata",]$species))))


##methods, sample collection part:
##for each group print out N of samples, species and orders
for (class_id in levels(stats_annot$color_class)){
  print(paste0(class_id, ": ", NROW(stats_annot[stats_annot$color_class==class_id, ]), " samples, ", 
              length(unique(stats_annot[stats_annot$color_class==class_id, ]$species)), " species, ",
              length(unique(stats_annot[stats_annot$color_class==class_id, ]$ncbi_order)), " orders" ))
}

##stats which locations things come from

location_stats <- locations %>% 
  group_by(Institution, Location) %>%
  summarise(N_species = length(unique(species)), N_samples = n())

write.csv(location_stats, file.path(analysis_dir, "01_basicStats","location_stats.csv") )
