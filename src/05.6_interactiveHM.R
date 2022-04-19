source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(treeio)
library(circilize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)


set.seed(123)
m = matrix(rnorm(100*100), 100)
rownames(m) = paste0("R", 1:100)
colnames(m) = paste0("C", 1:100)

ht1 = Heatmap(m, 
	top_annotation = HeatmapAnnotation(foo1 = runif(100)),
	left_annotation = rowAnnotation(bar1 = anno_points(1:100)),
	show_row_names = FALSE, show_column_names = FALSE)

m = matrix(rnorm(100*100), 100)
rownames(m) = paste0("R", 1:100)
colnames(m) = paste0("C", 1:100)

ht2 = Heatmap(m, 
	bottom_annotation = HeatmapAnnotation(foo2 = runif(100)),
	right_annotation = rowAnnotation(bar2 = anno_points(1:100)),
	show_row_names = FALSE, show_column_names = FALSE)

htShiny(ht1 + ht2, width1 = 600)
