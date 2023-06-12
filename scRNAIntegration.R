##########################################
##  Integrates scRNA experimental data  ##
##    with published gametocyte data    ##
##     and Malaria Cell Atlas data      ##
##########################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(plotly)

source('../util_funcs.R')

S.O.exp<-readRDS("../S.O.exp.rds")
S.O.pub<-readRDS("../S.O.pub.rds")
S.O.mca<-readRDS("../S.O.mca.rds")

S.O.exp.idc<-subset(S.O.exp, ident=c('0','1','2','4','5'))
S.O.exp.gam<-subset(S.O.exp, ident='3')

# integrate experimental data with MCA data
S.O.mca$orig.ident<-'mca'
S.O.exp.idc$orig.ident<-'exp'
S.O.list <- list(S.O.exp.idc = S.O.exp.idc, S.O.mca = S.O.mca)
idc_integrated<-integrate_S.O(S.O.list, 2)
Idents(idc_integrated)<-'stage'
DimPlot(idc_integrated, reduction = 'pca')
DimPlot(object = idc_integrated, label = TRUE, reduction = 'umap') + NoLegend()

# transfer labels from MCA to experimental data
Idents(idc_integrated) <- 'orig.ident'
exp_sub <- subset(idc_integrated, ident = 'exp')
mca_sub <- subset(idc_integrated, ident = 'mca')

anchors <- FindTransferAnchors(reference = mca_sub, query = exp_sub, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = mca_sub@meta.data$stage, dims = 1:30)
exp_sub <- AddMetaData(object = exp_sub, metadata = predictions)
exp_sub@meta.data$stage <- exp_sub@meta.data$predicted.id
Idents(exp_sub) <- 'stage'
DimPlot(exp_sub, reduction = 'pca')

ind1 <- idc_integrated@meta.data$orig.ident == 'exp'
ind2 <- match(rownames(idc_integrated@meta.data)[ind1], rownames(exp_sub@meta.data))
idc_integrated@meta.data$stage[ind1] <- exp_sub@meta.data$stage[ind2]

Idents(idc_integrated) <- 'stage'
idc_integrated$stage <- factor(idc_integrated$stage, levels = c('late_troph', 'early_troph', 'schizont', 'ring'))
idc_integrated@meta.data$stage <- factor(idc_integrated@meta.data$stage, 
                                         levels = c('late_troph', 'early_troph', 'schizont', 'ring'))
DimPlot(idc_integrated, reduction = 'pca')
DimPlot(object = idc_integrated, label = TRUE, reduction = 'umap') + NoLegend()

saveRDS(idc_integrated, "../S.O.integrated_idc.rds")

##############################

# integrate published data with experimental gametocyte data
S.O.pub$orig.ident<-'pub'
S.O.exp.gam$orig.ident<-'exp'
S.O.list <- list(S.O.exp.gam = S.O.exp.gam, S.O.pub = S.O.pub)
gam_integrated<-integrate_S.O(S.O.list, 2)
Idents(gam_integrated)<-'Cell_Type'
DimPlot(gam_integrated, reduction = 'pca')
DimPlot(object = gam_integrated, label = TRUE, reduction = 'umap') + NoLegend()

# transfer labels to experimental data
Idents(gam_integrated) <- 'orig.ident'
exp_sub <- subset(gam_integrated, ident = 'exp')
pub_sub <- subset(gam_integrated, ident = 'pub')

anchors <- FindTransferAnchors(reference = pub_sub, query = exp_sub, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = pub_sub@meta.data$Cell_Type, dims = 1:30)
exp_sub <- AddMetaData(object = exp_sub, metadata = predictions)
exp_sub@meta.data$Cell_Type <- exp_sub@meta.data$predicted.id
Idents(exp_sub) <- 'Cell_Type'
DimPlot(exp_sub, reduction = 'pca')

ind1 <- gam_integrated@meta.data$orig.ident == 'exp'
ind2 <- match(rownames(gam_integrated@meta.data)[ind1], rownames(exp_sub@meta.data))
gam_integrated@meta.data$Cell_Type[ind1] <- exp_sub@meta.data$Cell_Type[ind2]

gam_integrated$stage<-gam_integrated$Cell_Type
Idents(gam_integrated) <- 'stage'
gam_integrated$stage <- factor(gam_integrated$stage, levels = c('progenitor', 'female', 'male'))
gam_integrated@meta.data$stage <- factor(gam_integrated@meta.data$stage, 
                                         levels = c('progenitor', 'female', 'male'))
DimPlot(gam_integrated, reduction = 'pca')
DimPlot(object = gam_integrated, label = TRUE, reduction = 'umap') + NoLegend()

saveRDS(gam_integrated, "../S.O.integrated_gam.rds")

##############################

# integrate IDC object with gametocytes object
idc_integrated$orig.ident<-'idc'
gam_integrated$orig.ident<-'gam'
S.O.list <- list(idc_integrated = idc_integrated, gam_integrated = gam_integrated)
full_integrated<-integrate_S.O(S.O.list, 1)
Idents(full_integrated)<-'stage'
DimPlot(full_integrated, reduction = 'pca')

S.O.intra <- RunUMAP(gam_integrated, dims = 1:13, n.components = 3)
pca.full <- getPcaUmapMetaData(S.O.intra)

fig <- plot_ly(pca.full, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~stage, 
               colors = c('#f8766d', '#a2a400', '#00bf7d', '#04b0f6', '#e76bf3')) 
fig <- fig %>% add_markers(size = 2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))
fig
