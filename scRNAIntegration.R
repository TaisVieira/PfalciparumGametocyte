##########################################
##  Integrates scRNA experimental data  ##
##    with published gametocyte data    ##
##########################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(plotly)

source('../util_funcs.R')

S.O.exp <- readRDS("../S.O.exp.rds")
S.O.pub <- readRDS("../S.O.pub.rds")
S.O.mca <- readRDS("../S.O.mca.rds")

# divide experimental data into IDC and gametocyte
S.O.exp.idc <- subset(S.O.exp, ident = c('0', '1', '2', '4', '5'))
S.O.exp.gam <- subset(S.O.exp, ident = '3')

# get MCA labels for IDC
anchors <- FindTransferAnchors(reference = S.O.mca, query = S.O.exp.idc, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.mca@meta.data$stage, dims = 1:30)
S.O.exp.idc <- AddMetaData(object = S.O.exp.idc, metadata = predictions)
S.O.exp.idc@meta.data$stage <- S.O.exp.idc@meta.data$predicted.id

# integrate experimental gametocyte with published gametocyte
S.O.pub$orig.ident <- 'pub'
S.O.exp.gam$orig.ident <- 'exp'

anchors <- FindTransferAnchors(reference = S.O.pub, query = S.O.exp.gam, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = S.O.pub@meta.data$Cell_Type, dims = 1:30)
S.O.exp.gam <- AddMetaData(object = S.O.exp.gam, metadata = predictions)
S.O.exp.gam@meta.data$Cell_Type <- S.O.exp.gam@meta.data$predicted.id

S.O.list <- list(S.O.exp.gam = S.O.exp.gam, S.O.pub = S.O.pub)
S.O.int.gam <- integrate_S.O(S.O.list, 2)
S.O.int.gam$stage <- S.O.int.gam$Cell_Type
Idents(S.O.int.gam) <- 'stage'
DimPlot(S.O.int.gam, reduction = 'pca')
DimPlot(object = S.O.int.gam, label = TRUE, reduction = 'umap') + NoLegend()

# integrate IDC with gametocyte
S.O.list <- list(S.O.exp.idc = S.O.exp.idc, S.O.int.gam = S.O.int.gam)
S.O.int.full <- integrate_S.O(S.O.list, 1)
S.O.int.full$stage <- factor(S.O.int.full$stage, levels = c('ring', 'early_troph', 'late_troph', "schizont", 'progenitor', 'female', 'male'))
S.O.int.full@meta.data$stage <- factor(S.O.int.full@meta.data$stage, levels = c('ring', 'early_troph', 'late_troph', "schizont", 'progenitor', 'female', 'male'))
Idents(S.O.int.full) <- 'stage'
DimPlot(S.O.int.full, reduction = 'pca')
DimPlot(object = S.O.int.full, label = TRUE, reduction = 'umap') + NoLegend()

S.O.intra <- RunUMAP(S.O.int.full, dims = 1:13, n.components = 3)
pca.pf <- getPcaUmapMetaData(S.O.intra)

fig <- plot_ly(pca.pf, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~stage, 
               colors = c('#f8766d', '#a2a400', '#00bf7d', '#04b0f6', '#e76bf3')) 
fig <- fig %>% add_markers(size = 2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'), 
                                   yaxis = list(title = 'PC_2'), 
                                   zaxis = list(title = 'PC_3')))
fig

saveRDS(S.O.int.full, "../S.O.int.full.rds")
