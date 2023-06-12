#############################
## Creates a Seurat object ##
##     with scRNA data     ##
#############################

library(Seurat)
library(readxl)
library(Matrix)

source('../util_funcs.R')

input.dir <-  "C:/Users/Tais.Vieira/Desktop/research/Pf/scRNAseq_Pf/outs/"

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

write.csv(expr, paste(input.dir, "pf_scrna_count_table.csv", sep=""))

S.O.RNA<-CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
VlnPlot(S.O.RNA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.RNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
quant1<-quantile(S.O.RNA$nFeature_RNA, probs = c(0.1,0.98))
quant2<-quantile(S.O.RNA$nCount_RNA, probs = c(0.1,.98))
S.O.RNA<-subset(S.O.RNA, subset = nFeature_RNA > quant1[1] & nFeature_RNA < quant1[2] & 
                  nCount_RNA > quant2[1] & nCount_RNA < quant2[2])
S.O.RNA <- prep_S.O(S.O.RNA)
DimPlot(S.O.RNA, reduction = 'pca')
DimPlot(object = S.O.RNA, label = TRUE, reduction = 'umap') + NoLegend()

saveRDS(S.O.RNA, '../S.O.proc')
