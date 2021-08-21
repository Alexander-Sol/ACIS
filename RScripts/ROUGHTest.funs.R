library(aricode)
library(cidr)
library(RaceID)
library(SC3)
library(Seurat)
library(tidyverse)

Test.method <- function(object, method, best.IKAP = NULL, expr.meas = "umi"){
  
  if(method == "cidr"){
    n.clusters <- object@nCluster
    labels <- object@labels
    clusters <- as.factor(object@clusters)
    
  }else if(method == "raceid"){
    n.clusters <- max(object@cpart)
    labels <- object@cluster$labels
    clusters <- as.factor(object@cpart)
    
  }else if(method == "sc3"){
    # cluster.column <- 4
    # if (expr.meas != "umi") {cluster.column <- 2}
    # n.clusters <- max(as.integer(object@colData[ , cluster.column]))
    labels <- object@colData[["labels"]]
    clusters <- as.factor(object@colData[ , cluster.column])
    
  }else if(method == "autoclustr" | method == "cellfindr" | method == "seurat"){
    n.clusters <- length(levels(object))
    labels <- object$labels
    clusters <- as.factor(Idents(object))
  
  }else if(method == "ikap"){
    Idents(object) <- object[[object[["best.param"]][1,1]]]
    n.clusters <- length(levels(object))
    labels <- object$labels
    clusters <- as.factor(Idents(object))
    
  }else if(method == "celltrails"){
    clusters <- states(object)
    n.clusters <- length(levels(clusters))
    labels <- object@colData[["labels"]]
  }
  
  ari <- ARI(c1 = clusters, c2 = labels)

  return(value = list(n.clusters = n.clusters, ARI = ari))
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