source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

library(dummies)
library(MASS)
library(boot)
library(magrittr)
library(corrplot)
library(caret)
library(LiblineaR)


#custom pvalue matrix for pairwise correlation matrix
#needed for the corrplot function!

cor.mtest.cust <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


##draws and saves the pairwise correlation matrix
draw_feature_correlations<-function(feature_list, name, k){
  print(name)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  m<-scale(summary[,..feature_list], center=TRUE, scale=TRUE)
  
  res <- cor(m)
  p.mat<-cor.mtest.cust(m)
  
  pdf(paste0(analysis_dir, "/03_prediction_lm/03.1_features/", name, ".pdf"), width = k, height = k)
  corrplot(res, method="color", col=col(200), type="upper",order="hclust",
           addCoef.col = "black", addCoefasPercent = TRUE,# Add coefficient of correlation
           tl.col="black", tl.srt=45, #label parametrs
           p.mat = p.mat, sig.level = 0.05, 
           diag=FALSE)
  dev.off()
}


#function, that draws a scatter plot for two vectors, correlation line
#and the correlation and pvalues 
#used in 03.1

draw_correlation<-function(v1, v2, x_lab, y_lab){
  c<-cor.test(v1, v2)$estimate[["cor"]]
  p.val<-cor.test(v1, v2)$p.value
  
  print(c)

  p<-ggplot(summary, aes(x=v1, y=v2)) + 
    geom_point(aes( fill = color_class, alpha = 0.7, color = color_class)) + 
    scale_color_manual(values=class_colors)+scale_fill_manual(values=class_colors)+
    geom_smooth(method = "lm", size=0.7)+theme(legend.position = "None")+
    labs(x=x_lab, y=y_lab)
  
  return(p)
}


