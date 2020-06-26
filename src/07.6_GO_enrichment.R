#BiocManager::install("topGO")
source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
wd = file.path(analysis_dir, "03_motifAnalysis/03.7_AME")
setwd(wd)
outdir = file.path(wd, "GO")
dir.create(outdir)

library(topGO)

tf_db <- read.csv("../JASPAR.csv", sep = ";")

gene_clusters <- read.table("summary/clusters_and_selex.csv", sep = " ", header = 1)


##reading annotations:



##applying it to our TFs (all on the human level)
all_tf <- as.character(tf_db$Name)

all_tf <- as.character(sapply(all_tf, function(x) strsplit(x, "[(]")[[1]][1]))


require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
musGenes <- all_tf[sapply(all_tf, function(x) x!=toupper(x))]
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", musGenes, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- setNames(genesV2[, 2], genesV2[, 1])

all_tfs_id <- geneName2ID[c(all_tf[sapply(all_tf, function(x) x==toupper(x))], humanx)]
all_tfs_id <- unique(all_tfs_id[!is.na(all_tfs_id)])

gene_clusters$geneID <- geneName2ID[as.character(gene_clusters$TF_upd)]

###uploading gene2go
if(!file.exists("/scratch/lab_bock/dromanovskaia/resources/NCBI/gene2go_right_format_HS.csv")){
  gene2go <- read.table("/scratch/lab_bock/dromanovskaia/resources/NCBI/geneID2go_homo_sapiens", sep="\t", fill = F, quote="\"")
  colnames(gene2go) <- c("tax_id", "GeneID",	"GO_ID",	"Evidence", "Qualifier",	"GO_term",	"PubMed",	"Category")
  gene2gobygene <- gene2go %>%
    group_by(GeneID) %>%
    summarize(GOS = paste(GO_ID, collapse = ", "))
  write.table(gene2gobygene, "/scratch/lab_bock/dromanovskaia/resources/NCBI/gene2go_right_format_HS.csv", sep = "\t", row.names = F, quote = F)
}

gene2go_list <- readMappings("/scratch/lab_bock/dromanovskaia/resources/NCBI/gene2go_right_format_HS.csv")

transform <- function(s){
  l <- as.integer(length(s)/2)
  new_s <- paste(c(paste(s[c(1:l)], collapse = " "), paste(s[c((l+1):length(s))], collapse = " ")), collapse = "\n")
  return(new_s)
}

##desining the gene enricment for each cluster
runGOEA <- function(gene_set, NAME){
  #gene_set <- geneName2ID[toupper(strsplit(as.character(genes), ", ")[[1]])]
  ##creating gene universe
  geneList <- factor(as.integer(all_tfs_id %in% gene_set))
  names(geneList) <- all_tfs_id
  
  ##creating GOdata object
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, 
                annot = annFUN.gene2GO, gene2GO = gene2go_list)
  ##runiing fisher´s test
  resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  df <- as.data.frame(GenTable(GOdata, Fisher = resultFisher.elim))
  dTerm_full <- select(GO.db, keys = df$GO.ID, columns = c("TERM"))
  dTerm_full$TERM_nice <-sapply(dTerm_full$TERM, function(x) {s<-strsplit(x, " ")[[1]]; 
  if(length(s) < 6) {return (x)} else {return (transform(s))}} )
  df <- left_join(df, dTerm_full, by = c("GO.ID" = "GOID"))
  df$Fisher <- as.numeric(df$Fisher)
  df$logP <- -log10(df$Fisher)
  df <- df[order(df$logP, decreasing = F),]
  df$TERM_nice <- factor(df$TERM_nice, levels = df$TERM_nice)
  ggplot(df, aes(x = as.factor(TERM_nice), y = logP)) + geom_bar(stat = "identity") + 
    coord_flip() + labs(y = "-log10(p.Value)", x = "") 
  ggsave(file.path(outdir, paste0(NAME,"_bars.pdf")), width = 8, height = 5)
  write.table(GenTable(GOdata,  elim = resultFisher.elim), 
              file.path(outdir, paste0(NAME,".csv")), quote = F, row.names = F) 
  pdf(file.path(outdir, paste0(NAME,".pdf")), width = 15, height = 15)
  showSigOfNodes(GOdata, score(resultFisher.elim), firstSigNodes = 5, useInfo = 'all')
  dev.off()
  return(1)
}



runGOEA(gene_clusters[gene_clusters$cl == 1 & gene_clusters$Call_upd!="MethylPlus",]$geneID, "liver_act")
runGOEA(gene_clusters[gene_clusters$cl == 2 & gene_clusters$Call_upd!="MethylPlus",]$geneID, "heart_act")

runGOEA(gene_clusters[gene_clusters$cl == 1 & gene_clusters$Call_upd=="MethylPlus",]$geneID, "liver_methylPlus")

sapply(seq(NROW(gene_clusters)), 
       function(x) runGOEA(gene_clusters[x, "clusters"], as.character(x)))

runGOEA(paste(as.character(gene_clusters[2,]), as.character(gene_clusters[3,]), sep=", "), "all_orange")
runGOEA(paste(c(as.character(gene_clusters[1,]), as.character(gene_clusters[5,]), as.character(gene_clusters[4,])), collapse=", "), "all_blue")

