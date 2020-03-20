library(Seurat) # Seurat 3
library(dplyr)
library(Matrix)
library(reshape2)
#library(outliers)
library(gmodels)
library(RANN) # for nn2 function
library(Biobase)
library(gridExtra)
library(ggplot2)
library(cowplot)

objlist = list(GSM3489182,GSM3489187,
           GSM3489189,GSM3489191,GSM3489193,
           GSM3489195,Donor1)

anchors <- FindIntegrationAnchors(object.list = objlist, 
                                    dims = 1:30,
                                    anchor.features = 10000)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(object = combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- ScaleData(combined)

DefaultAssay(object = combined) <- "integrated"

combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
ElbowPlot(object = combined)

pc.n <- 15

combined <- RunTSNE(object = combined, reduction = "pca", dims = 1:pc.n)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:pc.n)

combined <- FindNeighbors(combined)
combined <- FindClusters(combined)

saveRDS(combined, file="LungControl.rds")
