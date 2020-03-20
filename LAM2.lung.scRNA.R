library(Seurat) # Seurat 2
library(dplyr)
library(Matrix)
library(reshape2)
library(outliers)
library(gmodels)
library(RANN) # for nn2 function
library(Biobase)
library(gridExtra)

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



sdp = readRDS(file="LAM2.rds")
sdp <- FindVariableGenes(sdp, do.plot = F)
genes.use <- rownames(subset(sdp@hvg.info, gene.mean>0.25 & gene.dispersion.scaled>0.25))
sdp@var.genes <- genes.use
cat("Reducing dimensions using ", length(sdp@var.genes), " highly variable genes\n", sep="")

sdp <- RunPCA(object = sdp, pc.genes = sdp@var.genes, 
              do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
PCElbowPlot(object = sdp)

pc.n <- 15
sdp <- RunTSNE(sdp, reduction.use = "pca", 
               dims.use = 1:pc.n, 
               tsne.method = "Rtsne")

sdp <- RunUMAP(sdp, reduction.use = "pca", 
               dims.use = 1:pc.n)

data.use <- sdp@dr$pca@cell.embeddings[, 1:pc.n]
ident <- graph.clustering.0(data.use, num.nn=10, do.jaccard=T, method="Louvain")

sdp <- SetIdent(sdp, ident.use=ident)
sdp@meta.data$louvain <- ident

saveRDS(sdp, file="LAM2.rds")
