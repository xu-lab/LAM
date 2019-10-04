library(Seurat)
library(dplyr)
library(Matrix)
library(reshape2)
library(outliers)
library(gmodels)
library(RANN)
library(Biobase)
library(gridExtra)
library(sSeq)


# Build a nearest neighbor graph with or without edge weights, and return an adjacency matrix
# code from Shekhar et al., 2016
get_edges=function(X,nn=30,do.jaccard=TRUE) {
  nearest=RANN::nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx = nearest$nn.idx[,-1]
  nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
  
  
  edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
  edges$B = edges$C; edges$C=1
  
  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  
  if (do.jaccard){
    
    NN = nearest$nn.idx
    jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    
    edges$C = jaccard_dist
    edges = subset(edges, C != 0)
    edges$C = edges$C/max(edges$C)
  }
  
  
  
  
  Adj = matrix(0, nrow=nrow(X), ncol=nrow(X))
  rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
  Adj[cbind(edges$A,edges$B)] = edges$C
  Adj[cbind(edges$B,edges$A)] = edges$C
  return(Adj)
  
}

# code adapted from Shekhar et al., 2016
graph.clustering.0 <- function(data.use, num.nn=30, do.jaccard=FALSE, method="Louvain") {
  
  
  if (do.jaccard){
    weights=TRUE;
    method_print = paste0(method,"-","Jaccard")
  } else {
    weights=NULL;
    method_print = method
  }
  
  print(paste0("Performing ", method_print, " clustering. Using ", num.nn, " nearest neighbors"))
  
  Adj = get_edges(data.use,nn=num.nn,do.jaccard=do.jaccard)
  
  
  g=igraph::graph.adjacency(Adj, mode = "undirected", weighted=weights)
  if (method=="Louvain") graph.out = igraph::cluster_louvain(g)
  if (method=="Infomap") graph.out = igraph::cluster_infomap(g)
  
  clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
  names(clust.assign) = graph.out$names
  k=order(table(clust.assign), decreasing = TRUE)
  new.levels = rep(1,length(unique(graph.out$membership)))
  new.levels[k] = 1:length(unique(graph.out$membership))
  levels(clust.assign) = new.levels
  clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
  print("Outputting clusters ..")
  
  return(clust.assign) 
}


min.gene.mean <- 0.1
min.gene.dispersion <- 0.25

genes.use <- c()
genes.commmon <- c()
objlist <- list()

cat("Processing dataset: LAM1\n")
obj <- readRDS(file=paste("LAM1.rds", sep=""))
obj@meta.data$group <- obj@meta.data$Group
obj@meta.data <- obj@meta.data[, c("nGene","nUMI","Cell","group", "celltype_1","celltype_2","celltype_3","celltype_4")]
obj <- FindVariableGenes(obj, do.plot=FALSE)
obj@var.genes <- rownames(subset(obj@hvg.info, gene.mean>min.gene.mean & gene.dispersion.scaled>min.gene.dispersion))
cat(length(obj@var.genes), "var genes\n")
genes.use <- c(genes.use, obj@var.genes)
genes.commmon <- rownames(obj@data)
objlist[["LAM1"]] <- obj


cat("Processing dataset: LAM3\n")
obj <- readRDS(file=paste("LAM3.rds", sep=""))
obj@meta.data$group <- obj@meta.data$Group
obj@meta.data <- obj@meta.data[, c("nGene","nUMI","Cell","group", "celltype_1","celltype_2","celltype_3","celltype_4")]
obj <- FindVariableGenes(obj, do.plot=FALSE)
obj@var.genes <- rownames(subset(obj@hvg.info, gene.mean>min.gene.mean & gene.dispersion.scaled>min.gene.dispersion))
cat(length(obj@var.genes), "var genes\n")
genes.use <- c(genes.use, obj@var.genes)
genes.commmon <- intersect(genes.commmon, rownames(obj@data))
objlist[["LAM3"]] <- obj


genes.use <- sort(unique(genes.use))
print(length(genes.use))
genes.use <- intersect(genes.use, genes.commmon)
print(length(genes.use))

num.cc <- 50
combined <- RunCCA(objlist[[1]], objlist[[2]], genes.use=genes.use, 
                   num.cc = num.cc)
combined@var.genes <- genes.use

acc.n <- 50
combined <- AlignSubspace(combined, reduction.type = "cca", 
                          grouping.var = "group", 
                          dims.align = 1:acc.n)

saveRDS(combined, file="LAM13.rds")

combined <- RunPCA(object = combined, pc.genes = combined@var.genes, 
                   do.print = FALSE, pcs.print = 1:5, pcs.compute = acc.n, 
                   genes.print = 5)
PCElbowPlot(object = combined)

combined <- RunUMAP(combined, reduction.use = "cca.aligned",
                    n_neighbors = 500,
                    min_dist = 0.3,
                    dims.use = 1:acc.n)

data.use <- combined@dr$cca.aligned@cell.embeddings[, 1:acc.n]
ident <- graph.clustering.0(data.use, num.nn=25, do.jaccard=T, method="Louvain")

combined <- SetIdent(combined, ident.use=ident)
combined@meta.data$louvain <- ident


# visualize the clustering results in umap plot
DimPlot(combined, reduction.use = "umap", do.return = T, pt.size = 0.5, 
              do.label = T, label.size = 5)


