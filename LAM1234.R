library(Seurat)
library(ggplot2)

db = "LAM1234"

# LAM1,...,LAM4 are Seurat3 objects of individual LAM lung data
objlist = list(LAM1,LAM2,LAM3,LAM4)

anchors <- FindIntegrationAnchors(object.list = objlist, 
                                  dims = 1:30, 
                                  anchor.features = 2000)
obj <- IntegrateData(anchorset = anchors, dims = 1:30)

obj = ScaleData(obj, vars.to.regress = "nCount_RNA")
obj = RunPCA(obj, npcs=200)

obj = RunUMAP(obj, reduction="pca", dims=1:30, umap.method="uwot", 
              min.dist=0.2, n.neighbors = 30, metric="cosine", n.epochs =500)

obj = FindNeighbors(obj, reduction = "pca", dims=1:100) 
obj = FindClusters(obj, resolution = 1.2)

saveRDS(obj, file=paste0(db, ".rds"))
