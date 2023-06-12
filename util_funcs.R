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
