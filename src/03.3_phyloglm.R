source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(ape)
library(gee)
library(patchwork)

wd = file.path(analysis_dir, "03_prediction_lm", "03.4_phyloglm")
dir.create(wd, recursive = TRUE)
setwd(wd)
# devtools::install_github("lamho86/phylolm")
## feature table
data <- as.data.frame(fread(file.path(analysis_dir, "01_basicStats", "01.10_all_seq_features", "feature_summary_filtered.tsv")))
row.names(data) <- data$species
## phylo tree
tree=ape::read.tree("../../01_basicStats/01.6_ITOL/tree_species_as_tips_2021_2.phy")

str(tree)

## functions 

### parces the comap.gee output
calculate_coefficions <- function(fit){
    nas <- is.na(fit$coef)

coef <- fit$coef[!nas]

cnames <- names(coef)
coef <- matrix(rep(coef, 4), ncol = 4)
dimnames(coef) <- list(cnames,
                       c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
df <- fit$dfP - dim(coef)[1]
coef[, 2] <- sqrt(diag(fit$W))
coef[, 3] <- coef[, 1]/coef[, 2]
if (df < 0) {
        warning("not enough degrees of freedom to compute P-values.")
        coef[, 4] <- NA
    } else coef[, 4] <- 2 * (1 -  pt(abs(coef[, 3]), df))
    
    return(coef)
}

run_on_subset <- function(subdf){
    colnames(subdf)[2] <- "kmer"
    sub_fit <- compar.gee(mean_meth ~ kmer,data = subdf, family = "gaussian", phy = tree)
    coef_sub <- as.data.frame(calculate_coefficions(sub_fit))
    coef_sub$name <- row.names(coef_sub)
    return(coef_sub)    
}

run_on_subset_nophy <- function(subdf){
    colnames(subdf)[2] <- "kmer"
    sub_fit_nophy <- glm(mean_meth ~ kmer,data = subdf, family = "gaussian")
    ans <- as.data.frame(coef(summary(sub_fit_nophy)))
    ans$name <- row.names(ans)
    return(ans)    
}


pval_per_kmer <- data.table()

## iterating over all 3-mers
for(kmer_name in colnames(data)[28:90]){
    
t1 <- run_on_subset(data[,c("mean_meth", kmer_name)])
t2 <-run_on_subset_nophy(data[,c("mean_meth", kmer_name)])

df <- left_join(t1[,c(1,4,5)],t2[,c(4,5)])[2,]
    
df$kmer <- kmer_name
colnames(df) <- c("Esitmate_phy", "pVal_phy", "name", "pVal_nophy", "kmer_name")

pval_per_kmer <- rbind(pval_per_kmer, df)
    
}

pval_per_kmer$pVal_phy_adj <- p.adjust(pval_per_kmer$pVal_phy, method = "bonferroni")
pval_per_kmer$pVal_nophy_adj <- p.adjust(pval_per_kmer$pVal_nophy, method = "bonferroni")

my_wt(pval_per_kmer, "pvalues_phy_vs_nonphy_per_kmer.tsv")

p1 <- ggplot(pval_per_kmer, aes(x = -log(pVal_phy_adj))) + geom_histogram(fill = "lightblue", color = "black") + geom_vline(xintercept = -log(0.05), linetype = "dashed", color = "red")
p1

p2 <- ggplot(pval_per_kmer, aes(x = -log(pVal_nophy_adj))) + geom_histogram(fill = "lightblue", color = "black") + geom_vline(xintercept = -log(0.05), linetype = "dashed", color = "red")
p2

g <- p1 | p2

ggsave("pval_adjusted_distribution.pdf", width = 8, height = 4)


ggplot(pval_per_kmer, aes( x =  -log(pVal_phy_adj), y =  -log(pVal_nophy_adj), color = Esitmate_phy)) + 
        geom_point() + 
        geom_text_repel(data = pval_per_kmer[pVal_phy_adj < 0.05 | pVal_nophy_adj<0.05], aes( x = -log(pVal_phy_adj), y =  -log(pVal_nophy_adj), label = kmer_name), max.overlaps = 25) + 
        scale_color_gradient2(low = "#313695",mid = "grey", high = "#a50026", midpoint = 0) + 
         geom_vline(xintercept = -log(0.05), linetype = "dashed", color = "red", alpha = 0.5) + 
        geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "red", alpha=0.5) + 
theme(text = element_text(size = 15))
ggsave("pval_adjusted_comparison.pdf", width = 6, height = 5)



## combined model (not used)
fit_nophy= glm(mean_meth~AAA + AAC + AAG + AAT + ACA + ACC + ACG + ACT + AGA + AGC + AGG + AGT + ATA + ATC + ATG + ATT + CAA + CAC + CAG + CAT + CCA + CCC + CCG + CCT + CGA + CGC + CGG + CGT + CTA + CTC + CTG + CTT + GAA + GAC + GAG + GAT + GCA + GCC + GCG + GCT + GGA + GGC + GGG + GGT + GTA + GTC + GTG + GTT + TAA + TAC + TAG + TAT + TCA + TCC + TCG + TCT + TGA + TGC + TGG + TGT + TTA + TTC + TTG,data=data,family="gaussian")
##




fit <- compar.gee(mean_meth ~ AAA + AAC + AAG + AAT + ACA + ACC + ACG + ACT + AGA + AGC + AGG + AGT + ATA + ATC + ATG + ATT + CAA + CAC + CAG + CAT + CCA + CCC + CCG + CCT + CGA + CGC + CGG + CGT + CTA + CTC + CTG + CTT + GAA + GAC + GAG + GAT + GCA + GCC + GCG + GCT + GGA + GGC + GGG + GGT + GTA + GTC + GTG + GTT + TAA + TAC + TAG + TAT + TCA + TCC + TCG + TCT + TGA + TGC + TGG + TGT + TTA + TTC + TTG,data = data, family = "gaussian", phy = tree)

summary(fit)




coef <- calculate_coefficions(fit)

coef_nophy <- coef(summary(fit_nophy))

df_pval <- data.frame(coef[,4], coef_nophy[,4])
colnames(df_pval) <- c("phy", "nophy")

ggplot(df_pval, aes(x = -log10(phy), y = -log10(nophy))) + geom_point() + geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed")  + geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") + ggtitle(paste0("corr = ", cor(df_pval$nophy, df_pval$phy)))




pval_per_kmer <- list()
for(kmer_name in colnames(data)[28:90]){
    print(kmer_name)
    pval_per_kmer[[kmer_name]] <- run_on_subset(data[,c("mean_meth", kmer_name)] )
    
}

df_kmers <- as.data.frame(t(data.frame(pval_per_kmer)))
df_kmers$kmer <- row.names(df_kmers)
df_kmers$pvalue_adj <- p.adjust(df_kmers$V1, method = 'bonferroni')

ggplot(df_kmers, aes(x = -log(V1))) + geom_histogram(fill = "lightblue", color = "black") + geom_vline(xintercept = -log(0.05), linetype = "dashed", color = "red")



row.names(df_kmers[-log(df_kmers$pvalue_adj) > 20,])
head(df_kmers)


## color the slope of ohylogenetic model?
