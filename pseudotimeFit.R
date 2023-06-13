##########################################
##   Separates cells of interest, fits  ##
##        a pseudotime curve using      ##
##   slingshot and creates a pseudotime ##
##       matrix of gene expression      ##
##########################################


library(Seurat)
library(slingshot)

source('C:/Users/Tais.Vieira/Desktop/research/scripts/util_funcs.R')

S.O.int.full <- readRDS("C:/Users/Tais.Vieira/Desktop/research/git/S.O.int.full.rds")

troph.prog <- subset(S.O.int.full, idents = c('early_troph', 'progenitor'))
sling_data <- slingshot(Embeddings(troph.prog, "pca"), clusterLabels = troph.prog$stage, start.clus = "early_troph")
sling_data <- SlingshotDataSet(sling_data)
plot3d.SlingshotDataSet(sling_data)
cell_id <- as.data.frame(slingBranchID(sling_data))

cells_troph <- rownames(sling_data@clusterLabels)[which(sling_data@clusterLabels[, 1] =  = 1)]
cells_prog <- rownames(sling_data@clusterLabels)[which(sling_data@clusterLabels[, 2] =  = 1)]
pseudotime <- as.data.frame(slingPseudotime(sling_data))
pseudotime$cell <- rownames(pseudotime)
pseudotime$stage <- NA
pseudotime[which(pseudotime$cell %in% cells_prog), 'stage'] <- 'progenitor'
pseudotime[which(pseudotime$cell %in% cells_troph), 'stage'] <- 'trophozoite'

ind <- colnames(S.O.int.full) %in% pseudotime$cell
sex_comm <- pseudotime$Lineage1
S.O.int.full@meta.data$Pseudotime_Sexual_Commit <- ifelse(colnames(S.O.int.full) %in% pseudotime$cell, sex_comm, NA)

my_cols <- c("#FFA500", "#440154FF")
FeaturePlot(S.O.int.full, features = "Pseudotime_Sexual_Commit", reduction = "pca", cols = my_cols, pt.size = 1.5)

gene_exp <- as.data.frame(S.O.int.full@assays$RNA@counts)
ind <- which(pseudotime$cell %in% rownames(cell_id))
sex_pseudotime <- pseudotime[ind, c("Lineage1", "cell", 'stage')]
temp <- sex_pseudotime[order(sex_pseudotime$Lineage1), ]
ordered_cells <- rownames(temp)
sex_pseudotime <- as.data.frame(t(sex_pseudotime))
ind2 <- which(colnames(gene_exp) %in% rownames(cell_id))
sex_pseudotime <- rbind(sex_pseudotime, gene_exp[, ind2])
sex_pseudotime <- sex_pseudotime[, ordered_cells]
saveRDS(sex_pseudotime, "../sex_pseudotime_expression.rds")
