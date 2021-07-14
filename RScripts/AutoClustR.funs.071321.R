Proc.data <- function(data, algorithm = "seurat", expr.meas = "umi", n.features = 2000, seed = NULL){

  library(cidr)
  library(Rtsne)
  library(Seurat)
  library(SingleCellExperiment)

  if(algorithm == "cidr"){

    if(expr.meas == "umi"){
      expr.meas <- "raw"
    }
    
    object <- scDataConstructor(tags = counts(data), tagType = expr.meas)
    scData.wLabs <- setClass("scData.wLabs", slots = c(labels = "vector"), contains = "scData")
    object <- as(object = object, Class = "scData.wLabs")
    object@labels <- colData(data)$labels
    
    object <- determineDropoutCandidates(object = object) %>%
      wThreshold() %>%
      scDissim() %>%
      scPCA(plotPC = FALSE)

  }else if(algorithm == "seurat"){

    object <- CreateSeuratObject(counts = counts(data))
    object@meta.data$labels <- colData(data)$labels
    object@meta.data$n.features <- colData(data)$n.features
    object@meta.data$pct.mito <- colData(data)$pct.mito
    
    if(expr.meas == "umi"){
      
      object <- NormalizeData(object = object, scale.factor = 1e6) %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress = c("pct.mito", "n.features")) 
      
    }else if(expr.meas == "tpm" | expr.meas == "cpm"){
      
      object@assays$RNA <- CreateAssayObject(data = log1p(counts(data)))
      object <- FindVariableFeatures(object = object, selection.method = "mvp") %>%
        ScaleData(vars.to.regress = c("pct.mito", "n.features"))
      
    }
    
    object <- RunPCA(object = object, seed.use = seed)

  }else if(algorithm == "tsne-kmeans"){

    object <- data
    
    if(expression.measure == "umi"){
      logcounts(object) <- log1p(1e6 * proportions(counts(object), margin = 2))
    }else if(expr.meas == "tpm" | expr.meas == "cpm"){
      logcounts(object) <- log1p(counts(object))
    }
    
    model.summary.table <- modelGeneVar(object)
    model.summary.table <- model.summary.table[order(model.summary.table[, 4], decreasing = TRUE), ]
    variable.features <- rownames(model.summary.table)[1:n.features]
    object <- object[variable.features, ]
    
    pca.results <- prcomp(t(logcounts(object)), scale. = TRUE)
    reducedDim(object, type = "pca") <- pca.results$x
    metadata(object)$roots <- pca.results$sdev

  }

  return(value = object)
}

ClustR <- function(object = object, n.pcs = 10, nCluster = NULL, k.param = 20, resolution = 0.8, centers = 10, perplexity = 30){

if(class(object) == "scData"){
  
  object <- scCluster(object = object, nCluster = nCluster, nPC = n.pcs)
  
  comp.scores <- object@PC[, 1:n.pcs]
  clusters <- object@clusters
  
}else if(class(object) == "Seurat"){
  
  # if(is.null(object@graphs[[paste0("RNA_snn_", k.param)]])){
     object <- FindNeighbors(object = object, dims = 1:n.pcs, k.param = k.param, prune.SNN = 1/15)
  #   object@graphs[[paste0("RNA_snn_", k.param)]] <- object@graphs$RNA_snn
  # }else{
  #   object@graphs$RNA_snn <- object@graphs[[paste0("RNA_snn_", k.param)]]
  # }
    
  object <- FindClusters(object = object, resolution = resolution)
  
  projections <- Embeddings(object = object)[, 1:n.pcs]
  clusters <- as.integer(Idents(object))
  
}else if(class(object) == "SingleCellExperiment"){
  
  reducedDim(object, type = "tsne") <- Rtsne(reducedDim(object, type = "pca")[, 1:n.pcs], perplexity = perplexity, pca = FALSE)$Y
  colData(object)$clusters <- kmeans(reducedDim(object, type = "tsne"), centers = centers)$cluster
  
  projections <- reducedDim(object)[, 1:n.pcs]
  clusters <- as.integer(colData(object)$clusters)
  
}

ICVI <- index.S(d = dist(projections, method = "manhattan"), cl = clusters)

if(!is.finite(ICVI)){
  ICVI <- -1
}

return(value = list(Score = ICVI))
}

AutoClustR <- function(object, file.path, method = "Bayesian", x.bounds = list(k.param = c(2, 160), resolution = c(0.0, 2.0)), n.priors = 15, n.starts = 15, n.pcs = NULL){
  
  library(clusterSim)
  library(lhs)
  library(lme4)
  library(ParBayesianOptimization)
  library(Seurat)
  library(stringi)
  
  if(is.null(n.pcs)){
    if(class(object) == "scData"){
      roots <- object@variation
    }else if(class(object) == "Seurat"){
      roots <- object@reductions[["pca"]]@stdev
    }else if(class(object) == "SingleCellExperiment"){
      roots <- metadata(object)$roots
    }
    n.pcs <- Select.nPC(roots, file.path = file.path)
  }
  
  n.pcs <- n.pcs[[2]]

  if(method == "Bayesian"){
    output <- bayesOpt(do.call(function(){ClustR(object = object, n.pcs = n.pcs, ...)}, list(k.param = k.param, resolution = resolution)),
      iters.n = n.starts - 1, 
      initPoints = n.priors, 
      bounds = x.bounds)
    
    x.best <- getBestPars(output)
    
  }else if(method == "Nelder-Mead"){
    
    x.mins <- sapply(x.bounds, "[", 1)
    x.maxes <- sapply(x.bounds, "[", 2)
    
    x0 <- maximinLHS(n = n.starts, k = length(x.bounds))
    
    for(i in 1:ncol(x0)){
      x0[, i] <- x0[, i] * (x.maxes[[i]] - x.mins[[i]]) + x.mins[[i]]
    }
    
    y.list <- list()
    x.list <- list()
    
    for(i in 1:n.starts){
      
      output <- Nelder_Mead(
        fn = function(x){
          eval(parse(text = paste0("do.call(ClustR, list(object = object, n.pcs = n.pcs, ", stri_join(names(x.bounds), stri_join("x[[", 1:length(x.bounds), "]]"), sep='=', collapse=','), "))")))},
        par = x0[i, ], lower = x.mins, upper = x.maxes)

      y.list[[i]] <- output$fval
      x.list[[i]] <- output$par
    }

    x.best <- x.list[[which.max(y.list)]]
        
  }

  
  return(value = list(object, method, n.pcs, x.best, x.bounds))
}
output <- AutoClustR(object = object, file.path = "C:/Users/15635/Documents/")