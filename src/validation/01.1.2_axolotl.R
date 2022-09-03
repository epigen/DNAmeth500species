source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

#The latest version of the genome and annotation are downloaded from https://www.axolotl-omics.org/assemblies

wd = file.path(data_dir, "resources", "reference_genomes","AmexG_v6.0-DD")
setwd(wd)

full_fasta = readDNAStringSet("AmexG_v6.0-DD_full.fa")

df <- fread("GCA_002915635.3_AmbMex60DD_assembly_report.txt", fill = TRUE, skip = 30)


df[`Assigned-Molecule-Location/Type`!= "na",sequence_name:=paste0("chromosome ", `Assigned-Molecule`), by = row.names(df[`Assigned-Molecule-Location/Type`!= "na"]) ]

df[`Assigned-Molecule-Location/Type`== "na",sequence_name:=`# Sequence-Name`, by = row.names(df[`Assigned-Molecule-Location/Type`== "na"]) ]


df$chr = paste0(df$`GenBank-Accn`, " Ambystoma mexicanum strain DD151 ", df$sequence_name, ", whole genome shotgun sequence")

write.table(df[, c("chr", "Sequence-Length")], "AmbMex60DD_full.chrom.sizes", quote = F, row.names = F, col.names = F, sep = "\t")

df <- df[`Assigned-Molecule-Location/Type`!= "na"]


for(nameid in names(full_fasta)){
    if(length(grep("chr", nameid)) == 1){
        writeXStringSet(full_fasta[nameid], "AmexG_v6.0-DD.fa", append=TRUE)
        print(nameid)
    }
}

write.table(df[, c("chr", "Sequence-Length")], "AmbMex60DD.chrom.sizes", quote = F, row.names = F, col.names = F, sep = "\t")


