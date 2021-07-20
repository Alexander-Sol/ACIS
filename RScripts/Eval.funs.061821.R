Evaluate.clustering <- function(object){

  # library(aricode)
  # library(cidr)
  # library(RaceID)
  # library(SC3)
  # library(Seurat)

  if(class(object) == "scData.wLabs"){

    labels <- object@labels
    n.clusters <- object@nCluster
    clusters <- as.factor(object@clusters)

  }else if(class(object) == "SCseq"){

    labels <- object@cluster$labels
    n.clusters <- max(object@cpart)
    clusters <- as.factor(object@cpart)
    
  }else if(class(object) == "SingleCellExperiment"){

    labels <- colData(object)$labels
    n.clusters <- length(levels(colData(object)[[grep(pattern = "\\clusters", colnames(colData(object)))]]))
    clusters <- colData(object)[[grep(pattern = "\\clusters", colnames(colData(object)))]]
    
  }else if(class(object) == "Seurat"){

    labels <- object@meta.data$labels
    n.clusters <- length(levels(Idents(object)))
    clusters <- Idents(object)

  }
  
  adj.Rand <- ARI(c1 = clusters, c2 = labels)
  #index.list <- Comp.ICVI()

  return(value = list(n.clusters, adj.Rand))
}



Comp.ICVI <- function(object, clusters, dim){
  
  # library(aricode)
  # library(clusterSim)
  # library(clValid)
  # library(fpc)
  # library(SingleCellExperiment)
  
  index.list <- c(
    ARI(c1 = colData(object)$labels, c2 = clusters),
    calinhara(reducedDim(object, type = "pca")[, 1:dim], clustering = clusters),
    calinhara(scale(reducedDim(object, type = "pca")[, 1:dim]), clustering = clusters),
    index.DB(reducedDim(object, type = "pca")[, 1:dim], cl = clusters)$DB,
    index.DB(reducedDim(object, type = "pca")[, 1:dim], cl = clusters, p = 1)$DB,
    index.DB(scale(reducedDim(object, type = "pca")[, 1:dim]), cl = clusters)$DB,
    index.DB(scale(reducedDim(object, type = "pca")[, 1:dim]), cl = clusters, p = 1)$DB,
    dunn(clusters = clusters, Data = reducedDim(object, type = "pca")[, 1:dim]),
    dunn(clusters = clusters, Data = reducedDim(object, type = "pca")[, 1:dim], method = "manhattan"),
    dunn(clusters = clusters, Data = scale(reducedDim(object, type = "pca")[, 1:dim])), #scale: set mean to 0, set SD to 1
    dunn(clusters = clusters, Data = scale(reducedDim(object, type = "pca")[, 1:dim]), method = "manhattan"),
    index.S(d = dist(reducedDim(object, type = "pca")[, 1:dim]), cl = clusters),
    index.S(d = dist(reducedDim(object, type = "pca")[, 1:dim], method = "manhattan"), cl = clusters),
    index.S(d = dist(scale(reducedDim(object, type = "pca")[, 1:dim])), cl = clusters),
    index.S(d = dist(scale(reducedDim(object, type = "pca")[, 1:dim]), method = "manhattan"), cl = clusters)
  )
  
  return(index.list)
}