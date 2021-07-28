Proc.data <- function(data, algorithm = "seurat", expr.meas = "umi", n.features = 2000, seed = NULL){

  # library(cidr)
  # library(Rtsne)
  # library(Seurat)
  # library(SingleCellExperiment)

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

  }else if (algorithm %in% c("seurat", "IKAP") ){

    object <- CreateSeuratObject(counts = SingleCellExperiment::counts(data))
    object@meta.data$labels <- colData(data)$labels
    object@meta.data$n.features <- colData(data)$n.features
    object@meta.data$pct.mito <- object@meta.data$percent.mito <- colData(data)$pct.mito
    object@meta.data$nUMI <- object@meta.data$nCount_RNA
    
    if (algorithm == "seurat") {
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
      
    } else if (algorithm == "IKAP") {
      if(expr.meas == "umi"){
        
        object <- NormalizeData(object = object, scale.factor = 1e6)
        
      }else if(expr.meas == "tpm" | expr.meas == "cpm"){
        
        object@assays$RNA <- CreateAssayObject(data = log1p(counts(data)))
        
      }
    }
    
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

ClustR <- function(object = object,
                   algorithm = "seurat",
                   n.pcs = 10,
                   nCluster = NULL,
                   k.param = NULL,
                   resolution= NULL,
                   centers = 10,
                   perplexity = 30,
                   seed = 0){
  
  if(class(object) == "scData" | class(object) == "scData.wLabs"){
    
    object <- scCluster(object = object, nCluster = nCluster, nPC = n.pcs)
    
    comp.scores <- object@PC[, 1:n.pcs]
    clusters <- object@clusters
    
  } else if(class(object) == "Seurat"){
    
    object <- FindNeighbors(object = object, dims = 1:n.pcs, k.param = k.param, prune.SNN = 1/15)
    print(paste("K.Param:", k.param))
      
    object <- FindClusters(object = object, resolution = resolution, random.seed = 534555234)
    print(paste("Resolution:", resolution))
    
    projections <- Embeddings(object = object)[, 1:n.pcs]
    clusters <- as.integer(Idents(object))
    
  } else if(class(object) == "SingleCellExperiment"){
    
    reducedDim(object, type = "tsne") <- Rtsne(reducedDim(object, type = "pca")[, 1:n.pcs], perplexity = perplexity, pca = FALSE)$Y
    colData(object)$clusters <- kmeans(reducedDim(object, type = "tsne"), centers = centers)$cluster
    
    projections <- reducedDim(object)[, 1:n.pcs]
    clusters <- as.integer(colData(object)$clusters)
    
  }
  
  # TODO: User should be able to select from whichever ICVI they choose
  ICVI <- index.S(d = dist(projections, method = "manhattan"), cl = clusters)
  
  if(!is.finite(ICVI)){
    ICVI <- -1
  }
  print(paste("Silhouette Score:", ICVI))
  return(value = list(Score = ICVI))
}

Comp.ICVI <- function(object, clusters, dim){
  
  # library(aricode)
  # library(clusterSim)
  # library(clValid)
  # library(fpc)
  # library(SingleCellExperiment)
  
  if(class(object) == "Seurat"){
    object <- as.SingleCellExperiment(object)
  }
  
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

AutoClustR <- function(object,
                       file.path,
                       method = "Bayesian",
                       x.bounds = list(k.param = c(2L, 160L),
                                       resolution= c(0.0, 2.0)),
                       n.priors = 15,
                       n.starts = 15,
                       max.iterations = 50,
                       n.pcs = NULL){
  
  # library(clusterSim)
  # library(lhs)
  # library(lme4)
  # library(ParBayesianOptimization)
  # library(Seurat)
  # library(stringi)
  
  # TODO: Going to have to add some checks to make sure input (i.e., x.bounds)
  #       is properly formatted. This is probably the most delicate part
  
  if(is.null(n.pcs)){
    if(class(object) == "scData"){
      roots <- object@variation
    }else if(class(object) == "Seurat"){
      roots <- object@reductions[["pca"]]@stdev
    }else if(class(object) == "SingleCellExperiment"){
      roots <- metadata(object)$roots
    }
    n.pcs <- Select.nPC(roots, file.path = file.path, do.plot = T)
  }
  
  n.pcs <- n.pcs[[2]]

  if(method == "Bayesian"){
    
    output <- ParBayesianOptimization::bayesOpt(
      FUN = function(...){ 
        ClustR(object = object, n.pcs = n.pcs, ...)
        },
      iters.n = n.starts - 1, 
      initPoints = n.priors, 
      bounds = x.bounds,
      acq = "ei"
      )
    
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
          x <- as.list(x)
          names(x) <- names(x.bounds)
          do.call(function(...){ClustR(object = object, n.pcs = n.pcs, ...)}, x)
          },
        par = x0[i, ],
        lower = x.mins,
        upper = x.maxes,
        control = list(maxfun = max.iterations) #Set max number of function calls
        )

      y.list[[i]] <- output$fval
      x.list[[i]] <- output$par
    }

    x.best <- x.list[[which.max(y.list)]]
        
  }
  
  return(value = list(object, method, n.pcs, x.best, x.bounds))
}
# output <- AutoClustR(object = object, file.path = "C:/Users/15635/Documents/")

source("RScripts/nPC.funs.R")

#Should be sourcing stuff here, but
# TODO: Add verbosity option to AutoClustR
b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv")
b1 <- Prep.data(b1)
b1 <- Proc.data(b1)

kolodz <- Kolodz.to.SCexp("Data/Kolodz/counttable_es.csv")
kolodz <- Prep.data(kolodz)
kolodz <- Proc.data(kolodz, expr.meas = "cpm")
ACTest <- AutoClustR(kolodz,
                     file.path = "",
                     method = "Bayesian",
                     n.priors = 16,
                     n.starts = 16)
sc.runtime <- Sys.time()
sc.test <- subcluster(ACTest)
sc.runtime <- Sys.time() - sc.runtime

subcluster <- function(ACTest) {
  ac.obj <- FindNeighbors(ACTest[[1]], k.param = ACTest[[4]]$k.param, dims = 1:ACTest[[3]]) %>%
    FindClusters(resolution = ACTest[[4]]$resolution)
  n.pcs <- ACTest[[3]]
  initial.eval <- Test.method(ac.obj, method = "autoclustr")
  accept = T
  
  while(accept){
    clusterSilhouettes <- getSilScores(ac.obj, npc = 8)
    check.clusters <- TRUE
    cl.iter <- 1
    while(check.clusters) {
      sc.object <- subset(ac.obj, idents = clusterSilhouettes[cl.iter, 1])
      sc.k.param.range <- k.param.range <- c(2,140)
      if (k.param.range[2] > ncol(sc.object)/2) {sc.k.param.range[2] <- round(ncol(sc.object)/2) }
      if (sc.k.param.range[2] <= 4) {cl.iter <- cl.iter + 1} else {check.clusters <- FALSE}
    }
    x.bounds <- list(k.param = sc.k.param.range, resolution = c(0.0, 2.0))
    
    sc.results <- ParBayesianOptimization::bayesOpt(
      FUN = function(...){ 
        ClustR(object = sc.object, n.pcs = n.pcs, ...)
      },
      iters.n = 8 - 1, 
      initPoints = 8, 
      bounds = x.bounds,
      acq = "ei"
    )
    
    sc.object <- bestBayes(sc.object, sc.results, n.pcs = n.pcs)
    adjust.output <- adjustIdents(ac.obj, temp.object = sc.object, npc = n.pcs)
    ac.obj <- adjust.output$object
    accept <- adjust.output$accept
  }
  
  idents <- as.int.factor(ac.obj@active.ident) # as.int.factor should probably be a class function
  levels(ac.obj@active.ident) <- 0:sum(1:max(idents) %in% unique(idents))
  
  return(
    list(
      object = ac.obj,
      initial.eval = initial.eval,
      final.eval = Test.method(ac.obj, method = "autoclustr")
      )
    )
}

bestBayes <- function(object, bayes.out, n.pcs) {
  best.params <- getBestPars(bayes.out)
  object <- FindNeighbors(object, k.param = best.params$k.param, dims = 1:n.pcs) %>%
    FindClusters(resolution = best.params$resolution)
  return(object)
}

getSilScores <- function(object, npc) {
  embeds <- Embeddings(object)[ , 1:npc]
  embed.table <- cbind(embeds,
                       as.integer(levels(object@active.ident))[object@active.ident])
  sil.scores <- cluster::silhouette(as.integer(embed.table[ , ncol(embed.table)]),
                                    dist(embed.table[, 1:(ncol(embed.table-1))]))
  sil <- data.frame(cluster = sil.scores[,1], si_width = sil.scores[,3])
  scores.by.cluster <- aggregate(sil$si_width, list(sil$cluster), mean) %>%
    arrange(x)
  names(scores.by.cluster) <- c("Seurat_cluster", "avg.SI.score")
  return(scores.by.cluster)
}

getICVI <- function(object, nPC = NULL) {
  icvi <- suppressWarnings(
    index.S(
      d = dist(
        Embeddings(object, reduction = "pca")[, 1:nPC],
        method = "manhattan"),
      cl = as.integer(Idents(object))
      )
    )
  if(is.nan(icvi)) {icvi <- 0}
  return(icvi)
}

as.int.factor <- function(x) {as.integer(levels(x))[x]}

adjustIdents <- function(object, temp.object, npc) {
  embeds <- Embeddings(object)[ , 1:npc]
  new.idents <- as.int.factor(temp.object@active.ident) + 
    1 + 
    max(as.int.factor(object@active.ident))
  names(new.idents) <- Cells(temp.object)
  combined.new <- as.int.factor(object@active.ident)
  names(combined.new) <- Cells(object)
  combined.new[names(new.idents)] <- new.idents
  original.SI <- getICVI(object, npc)
  # original.SI <- intCriteria(embeds, as.int.factor(object@active.ident), crit = "sil")[[1]]
  new.SI <-  suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:npc],
                                               method = "manhattan"),
                                      cl = combined.new %>% as.factor %>% as.numeric))
  # new.SI <- intCriteria(embeds, as.integer(combined.new), "sil")[[1]]
  if(original.SI > new.SI) {
    print(paste(object@project.name, "REJECTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
    accept <- FALSE
  } else if (original.SI <= new.SI) {
    object@active.ident <- as.factor(combined.new)
    print(paste(object@project.name, "ACCEPTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
    accept <- TRUE
  }
  return(list(object=object, accept = accept))
}
