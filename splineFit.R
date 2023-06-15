############################################
##   Creates gene expression plots using  ##
##    predicted values from pseudotime    ##
##     gene expression using splines      ##
############################################

library(Seurat)
library(dplyr)
library(pheatmap)
library(dtwclust)

S.O.int.full <- readRDS("./S.O.int.full.rds")
troph.prog <- subset(S.O.int.full, idents = c('early_troph', 'progenitor'))
DefaultAssay(troph.prog) <- 'RNA'

# create expression matrix of DEG
markers_troph_prog <- FindAllMarkers(troph.prog, only.pos = TRUE)
markers_troph_prog_filter <- markers_troph_prog %>% filter(avg_log2FC > 1 & pct.1 > 0.5)
exp_mat <- as.matrix(troph.prog@assays$RNA@data)
genes <- which(rownames(exp_mat) %in% rownames(markers_troph_prog_filter))
exp_mat_filter <- exp_mat[genes, ]
dim(exp_mat_filter)

sex_pseudotime <- readRDS("../sex_pseudotime_expression.rds")

cells <- colnames(sex_pseudotime)
exp_mat_filter <- exp_mat_filter[, cells]

stages <- as.data.frame(t(sex_pseudotime['stage', ]))
newCols  <- colorRampPalette(grDevices::rainbow(length(unique(stages$stage))))
annoCol  <- newCols(length(unique(stages$stage)))
names(annoCol)  <- unique(stages$stage)
annoCol  <- list(category = annoCol)
heat_map <- pheatmap(exp_mat_filter, annotation_col = stages, annotation_colors = annoCol, scale = "column", 
                   cluster_cols = FALSE, show_rownames = F, show_colnames = F)
genes_order_heatmap <- rownames(exp_mat_filter[heat_map$tree_row[["order"]], ])

# predict gene expression using splines
max_time <- as.numeric(sex_pseudotime[1, ncol(sex_pseudotime)])
pseudotime <- as.numeric(sex_pseudotime[1, ])
scaled_pseudo <- pseudotime/max_time

first_gene <- exp_mat_filter[genes_order_heatmap[1], ]
first_gene_spline <- data.frame(scaled_pseudo, first_gene)
spline <- smooth.spline(first_gene_spline$scaled_pseudo, first_gene_spline$first_gene)
prediction  <- predict(spline, seq(0, 1, by = 0.01))
exp_mat_filter_pred <- as.data.frame(prediction$y)

for (i in 2:nrow(exp_mat_filter)){
  pred_gene <- exp_mat_filter[genes_order_heatmap[i], ]
  pred_gene_spline <- data.frame(scaled_pseudo, pred_gene)
  spline <- smooth.spline(pred_gene_spline$scaled_pseudo, pred_gene_spline$pred_gene)
  prediction  <- predict(spline, seq(0, 1, by = 0.01))
  exp_mat_filter_pred <- cbind(exp_mat_filter_pred, prediction$y) 
}

colnames(exp_mat_filter_pred) <- genes_order_heatmap
exp_mat_filter_pred <- t(exp_mat_filter_pred)
colnames(exp_mat_filter_pred) <- seq(0, 1, 0.01)
write.csv(exp_mat_filter_pred, "../predicted_exp_mat_troph_prog.csv")

# plot predicted gene expression
exp_mat_filter_pred <- read.csv("../predicted_exp_mat_troph_prog.csv")
rownames(exp_mat_filter_pred) <- exp_mat_filter_pred$X
exp_mat_filter_pred$X <- NULL
gene_exp <- tsclust(exp_mat_filter_pred, type = 'h', k = 13L, distance = 'dtw', centroid = shape_extraction, preproc = zscore, 
                  trace = T)
plot(gene_exp, type = 'sc')
