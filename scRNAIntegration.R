################################################
##   Integrates gametocyte data from scRNA    ##
## experiment with published gametocyte data  ##
################################################

library(Seurat)
library(ggplot2)

source('../util_funcs.R')

gam.exp<-readRDS("../S.O.exp")
gam.pub<-readRDS("../S.O.pub")

# integrate both objects
S.O.list <- list(gam.exp = gam.exp, gam.pub = gam.pub)
S.O.list <- lapply(S.O.list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = S.O.list)
spps <- names(S.O.list)
reference_dataset <- 2
anchors <- FindIntegrationAnchors(object.list = S.O.list, anchor.features = features, dims = 1:30)
gam_integrated <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 65)

# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(gam_integrated) <- "integrated"
gam_integrated <- ScaleData(object = gam_integrated, verbose = FALSE)
gam_integrated <- RunPCA(gam_integrated, features = VariableFeatures(object = gam_integrated))
gam_integrated <- FindNeighbors(gam_integrated, dims = 1:10, reduction = 'pca')
gam_integrated <- FindClusters(gam_integrated, resolution = 0.5)
gam_integrated <- RunUMAP(gam_integrated, dims = 1:13)
Idents(gam_integrated)<-'Cell_Type'
DimPlot(gam_integrated, reduction = 'pca')
DimPlot(object = gam_integrated, label = TRUE, reduction = 'umap') + NoLegend()

# transfer labels to experimental data
gam_integrated@meta.data$orig.ident<-ifelse(gam_integrated@meta.data$orig.ident == "NF54", "pub", 'exp')
Idents(gam_integrated) <- 'orig.ident'
gam1_sub <- subset(gam_integrated, ident = 'exp')
gam2_sub <- subset(gam_integrated, ident = 'pub')

anchors <- FindTransferAnchors(reference = gam2_sub, query = gam1_sub, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = gam2_sub@meta.data$Cell_Type,dims = 1:30)
gam1_sub <- AddMetaData(object = gam1_sub, metadata = predictions)
gam1_sub@meta.data$Cell_Type <- gam1_sub@meta.data$predicted.id
Idents(gam1_sub) <- 'Cell_Type'
DimPlot(gam1_sub, reduction = 'pca')

ind1 <- gam_integrated@meta.data$orig.ident == 'exp'
ind2 <- match(rownames(gam_integrated@meta.data)[ind1], rownames(gam1_sub@meta.data))
gam_integrated@meta.data$Cell_Type[ind1] <- gam1_sub@meta.data$Cell_Type[ind2]

Idents(gam_integrated) <- 'Cell_Type'
gam_integrated$Cell_Type <- factor(gam_integrated$Cell_Type, levels = c('progenitor', 'female', 'male'))
gam_integrated@meta.data$Cell_Type <- factor(gam_integrated@meta.data$Cell_Type, levels = c('progenitor', 'female', 'male'))
p <- DimPlot(gam_integrated, reduction = "pca", 
             pt.size = 1,
             label = T, label.size = 5) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

saveRDS(gam_integrated, "../S.O.integrated_gam.rds")
