# PfalciparumGametocyte
Exploration of Plasmodium falciparum gametocyte's genetic expression profile.

* scRNAProcessing: creates a Seurat object with scRNA data. Experimental data came from a 10X Genomic scRNAseq experiment provided by Dr. David Serre, published data came from the published paper "A transcriptional switch controls sex determination in Plasmodium falciparum" (https://doi.org/10.1038/s41586-022-05509-z) and MCA data came from the Malaria Cell Atlas (https://github.com/vhowick/MalariaCellAtlas).
* scRNAIntegration: integrates experimental data with published gametocyte data using MCA stage labels for IDC. 
* pseudotimeFit: separates cells of interest from the integrated object, fits a pseudotime curve using slingshot and creates a pseudotime matrix of gene expression.
* splineFit: creates gene expression plots using predicted values from pseudotime gene expression using splines.
* ap2gPromoter: checks if the promoter of genes of interest have AP2G and/or AP2G-2 motifs. Promoters are considered as the 1kbp upstream of the gene.
* util_funcs: functions prep_S.O, integrate_S.O, getPcaUmapMetaData .
 
