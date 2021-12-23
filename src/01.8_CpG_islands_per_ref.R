source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

args = commandArgs(trailingOnly=TRUE)
species = args[1]


wd = file.path(processed_dir,species)
setwd(wd)
ded_ref = readDNAStringSet("reduced/consensus/toSelf_filtered_0.08mm_final.fa")


C <- letterFrequency(ded_ref, "C")/width(ded_ref)
G <- letterFrequency(ded_ref, "G")/width(ded_ref)

freq <- dinucleotideFrequency(ded_ref)/(width(ded_ref)-1)
CG_freq <- freq[, "CG"]
  
GC_content <- C+G
CpG_OE <- CG_freq/(C*G)

criteria <- data.frame(GC_content, CpG_OE)
colnames(criteria) <- c("CG_cont", "CpG_OE")
criteria["GG"] <- (criteria$CG_cont>0.5)&(criteria$CpG_OE>0.6)
criteria["TJ"] <- (criteria$CG_cont>0.55)&(criteria$CpG_OE>0.65)

ans <- data.frame(NROW(criteria[criteria$GG,])/NROW(criteria),NROW(criteria[criteria$TJ,])/NROW(criteria))
colnames(ans) <- c("GardG", "TJ")
row.names(ans) <- species

write.csv(ans, paste0(analysis_dir, "/01_basicStats/01.8_CpG_islands_per_ref/", species,".csv" ), quote = F)
