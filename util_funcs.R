prep_S.O <- function(S.O, res = 0.1, var.features = F, down.sample = F){
  set.seed(100)
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 6000)
  if(var.features){
    ## Work on variable features only
    S.O <- subset(S.O, features = VariableFeatures(S.O))
  }
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:10, reduction = 'pca')
  S.O <- FindClusters(S.O, resolution = res)
  if(down.sample){
    S.O <- subset(x = S.O, downsample = 800)
    S.O <- FindNeighbors(S.O, dims = 1:13)
    S.O <- FindClusters(S.O, resolution = res)
  }
  S.O <- RunUMAP(S.O, dims = 1:13, n.components = 3L)
  return(S.O)
}

integrate_S.O<-function(S.O.list, ref.dataset){
  set.seed(100)
  S.O.list <- lapply(S.O.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA'
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  
  features <- SelectIntegrationFeatures(object.list = S.O.list)
  spps <- names(S.O.list)
  reference_dataset <- ref.dataset
  anchors <- FindIntegrationAnchors(object.list = S.O.list, anchor.features = features, dims = 1:30)
  S.O.integrated <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 65)
  
  # switch to integrated assay. Make sure to set to RNA for Differential Expression
  DefaultAssay(S.O.integrated) <- "integrated"
  S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
  S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
  S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
  S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.5)
  S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)
  
  return(S.O.integrated)
}

getPcaUmapMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2, UMAP_3 = UMAP_3)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), clusters = S.O@meta.data$seurat_clusters,
                          stage = S.O@meta.data$stage, orig.ident = S.O@meta.data$orig.ident)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}
