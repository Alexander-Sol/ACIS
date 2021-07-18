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

    labels <- object@metadata$labels
    n.clusters <- length(levels(Idents(object)))
    clusters <- Idents(object)

  }
  
  adj.Rand <- ARI(c1 = clusters, c2 = labels)

  return(value = list(n.clusters, adj.Rand))
}