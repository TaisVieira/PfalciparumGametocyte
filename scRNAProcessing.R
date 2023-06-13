#############################
## Creates a Seurat object ##
##     with scRNA data     ##
#############################

library(Seurat)
library(readxl)
library(Matrix)

source('../util_funcs.R')

# experiment data
input.dir <-  "../scRNAseq_Pf/outs/"

barcode.path <- paste(input.dir, "filtered_feature_bc_matrix/barcodes.tsv.gz", sep = "")
features.path <- paste(input.dir, "filtered_feature_bc_matrix/features.tsv.gz", sep = "")
matrix.path <- paste(input.dir, "filtered_feature_bc_matrix/matrix.mtx.gz", sep = "")

mat <- readMM(file = matrix.path)

feature.names = read.delim(features.path, 
                           header = FALSE, 
                           stringsAsFactors = FALSE)

barcode.names = read.delim(barcode.path, 
                           header = FALSE, 
                           stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

expr <- as.data.frame(as.matrix(mat))

S.O.RNA <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
VlnPlot(S.O.RNA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.RNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
quant1 <- quantile(S.O.RNA$nFeature_RNA, probs = c(0.1, 0.98))
quant2 <- quantile(S.O.RNA$nCount_RNA, probs = c(0.1, 0.98))
S.O.RNA <- subset(S.O.RNA, subset = nFeature_RNA > quant1[1] & nFeature_RNA < quant1[2] & 
                  nCount_RNA > quant2[1] & nCount_RNA < quant2[2])
S.O.RNA <- prep_S.O(S.O.RNA)
DimPlot(S.O.RNA, reduction = 'pca')
DimPlot(object = S.O.RNA, label = TRUE, reduction = 'umap') + NoLegend()

saveRDS(S.O.RNA, '../S.O.exp.rds')

# published data
gam_file <- read.csv("../scRNAseq_Pf/Fig1-10X-data.csv", sep = ";")

gam_file <- t(gam_file)
colnames(gam_file) <- gam_file[1, ]
gam_file <- as.data.frame(gam_file[-1, ])

counts <- gam_file[which(gam_file$Genotype  =  =  "WT"), 5:ncol(gam_file)]
meta <- gam_file[, 1:4]

expr <- t(counts)

S.O.gam <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
VlnPlot(S.O.gam, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.gam, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
quant1 <- quantile(S.O.gam$nFeature_RNA, probs = c(0.1, 0.98))
quant2 <- quantile(S.O.gam$nCount_RNA, probs = c(0.1, 0.98))
S.O.gam <- subset(S.O.gam, subset = nFeature_RNA > quant1[1] & nFeature_RNA < quant1[2] & 
                  nCount_RNA > quant2[1] & nCount_RNA < quant2[2])
S.O.gam <- prep_S.O(S.O.gam)
DimPlot(S.O.gam, reduction = 'pca')
DimPlot(object = S.O.gam, label = TRUE, reduction = 'umap') + NoLegend()

S.O.gam <- AddMetaData(S.O.gam, meta)
Idents(S.O.gam) <- 'Cell_Type'

saveRDS(S.O.gam, '../S.O.pub.rds')

# MCA data
file.counts <- read.csv("../scRNAseq_Pf/MalariaCellAtlas-master/Expression_Matrices/10X/pf10xIDC/pf10xIDC_counts.csv")
genes <- file.counts$X
expr <- file.counts[, -1]
rownames(expr) <- genes

S.O.MCA <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
VlnPlot(S.O.MCA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.MCA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
quant1 <- quantile(S.O.MCA$nFeature_RNA, probs = c(0.1, 0.98))
quant2 <- quantile(S.O.MCA$nCount_RNA, probs = c(0.1, 0.98))
S.O.MCA <- subset(S.O.MCA, subset = nFeature_RNA > quant1[1] & nFeature_RNA < quant1[2] & 
                  nCount_RNA > quant2[1] & nCount_RNA < quant2[2])
S.O.MCA <- prep_S.O(S.O.MCA)
DimPlot(S.O.MCA, reduction = 'pca')
DimPlot(object = S.O.MCA, label = TRUE, reduction = 'umap') + NoLegend()

meta_data <- read.csv("../scRNAseq_Pf/MalariaCellAtlas-master/Expression_Matrices/10X/pf10xIDC/pf10xIDC_pheno.csv")
cell_id <- meta_data$X
meta_data <- meta_data[, -1]
rownames(meta_data) <- cell_id
meta_data <- meta_data[, c('sample_id', 'clock_pseudotime', 'bulk')]
colnames(meta_data) <- c('sample_id', 'clock_pseudotime', 'stage')
S.O.MCA <- AddMetaData(S.O.MCA, meta_data)
Idents(S.O.MCA) <- 'stage'
S.O.MCA@meta.data$group <- 'idc'
DimPlot(S.O.MCA, reduction = 'pca')

saveRDS(S.O.gam, '../S.O.mca.rds')
