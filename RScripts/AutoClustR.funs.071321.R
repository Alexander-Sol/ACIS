# library(clusterSim) clusterSim is causing Xquartz to launch in dock when loaded into environment (on OSX Big Sur, R 4.1.0)
# It may be one of the dependencies, either cluster or MASS packages
library(clusterSim)
library(lhs)
library(lme4)
library(ParBayesianOptimization)
library(stringi)
library(cidr)
library(Rtsne)
library(Seurat)
library(SingleCellExperiment)


# acq.fun = acq argument for ParBayesOpt. Can be one of "ucb", "ei", "eips", "poi"
#Bayes op requires at least 4 points

AutoClustR <- function(object,
                       # file.path,
                       method = "Bayesian",
                       x.bounds = list(
                         k.param = c(2L, 160L),
                         resolution= c(0.0, 2.0)
                         ),
                       n.priors = 12,
                       n.starts = 12,
                       max.iterations = 50,
                       n.pcs = NULL,
                       acq.func = "ucb",
                       subcluster = T){
  
  # TODO: Going to have to add some checks to make sure input (i.e., x.bounds)
  #       is properly formatted. This is probably the most delicate part
  
  if(is.null(n.pcs)) {
    if(class(object) == "scData"){
      roots <- object@variation
    }else if(class(object) == "Seurat"){
      roots <- object@reductions[["pca"]]@stdev
    }else if(class(object) == "SingleCellExperiment"){
      roots <- metadata(object)$roots
    }
    n.pcs <- Select.nPC(roots,  do.plot = T) #TODO: Remove file.path argument, also maybe don't show the n.pc vs Max.PC plot unless they ask
  }

  if(method == "Bayesian") {
    
    output <- ParBayesianOptimization::bayesOpt(
      FUN = function(...){ 
        ClustR(object = object, n.pcs = n.pcs, ...)
        },
      iters.n = n.starts - 1, 
      initPoints = n.priors, 
      bounds = x.bounds,
      acq = acq.func
      )
    
    object <- bestBayes(object, output, n.pcs)
    
  }
  
  if(subcluster) {
    sc.results <- subClustR(object, 
                            x.bounds = x.bounds,
                            n.priors = n.priors,
                            n.starts = n.starts,
                            n.pcs = n.pcs, 
                            acq.func = acq.func)
    object <- sc.results$object
    pre.sc <- sc.results$initial.eval
    post.sc <- sc.results$final.eval
    sc.iterations <- sc.results$accept.count
  }
  
  return(
    list(
      object = object,
      method = method,
      n.pcs = n.pcs,
      pre.sc = pre.sc,
      post.sc = post.sc,
      sc.iterations = sc.iterations,
      x.best = getBestPars(output),
      x.bounds = x.bounds
    )
  )
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

#subClustR should also work as a stand alone function. Maybe we should just print the results
#i.e., "subClustR split two (2) clusters into five (5) sub clusters
#Currently, subClustR halts whenever a sub-clustered solution proves worse than the previous clustering solution
#We should give an option to tune that parameter in some way
subClustR <- function(object,
                      x.bounds = list(
                        k.param = c(2L, 160L),
                        resolution= c(0.0, 2.0)
                      ),
                      n.priors = 12,
                      n.starts = 12,
                      n.pcs = NULL,
                      acq.func = "ei") {

  # Test.method returns list of n.clusters and ARI
  # So this should remain internal/ be removed before package is published
  initial.eval.Test <- Test.method(object, method = "autoclustr") 
  
  projections <- Embeddings(object = object)[, 1:n.pcs] #projections won't change, so can be stored and re=used
  initial.clusters <- as.integer(Idents(object))
  initial.ICVI <- index.S(d = dist(projections, method = "manhattan"), cl = initial.clusters)
  accept = T
  accept.count = 0
  
  while(accept) {
    
    clusterSilhouettes <- getSilScores(object, npc = n.pcs) #Returns ordered list of seurat clusters from worst to best (as ranked by SI)
    print(clusterSilhouettes)
    # This code checks to ensure that only clusters > ~8 cells are submitted for 
    # sub clustering. 
    #TODO: Clean up this code
    check.clusters <- TRUE
    cl.iter <- 1
    while(check.clusters) {
      sc.object <- subset(object, idents = clusterSilhouettes[cl.iter, 1])
      sc.k.param.range <- k.param.range <- c(2,140)
      if (k.param.range[2] > ncol(sc.object)/2) {sc.k.param.range[2] <- round(ncol(sc.object)/2) }
      if (sc.k.param.range[2] <= 4) {cl.iter <- cl.iter + 1} else {check.clusters <- FALSE}
    }
    
    x.bounds <- list(k.param = sc.k.param.range, resolution = c(0.0, 2.0))
    
    sc.output <- ParBayesianOptimization::bayesOpt(
      FUN = function(...){ 
        ClustR(object = sc.object, n.pcs = n.pcs, ...)
      },
      iters.n = n.starts, 
      initPoints = n.priors, 
      bounds = x.bounds,
      acq = acq.func
    )
    
    sc.object <- bestBayes(sc.object, sc.output, n.pcs = n.pcs)
    
    adjust.output <- adjustIdents(object, temp.object = sc.object, n.pcs = n.pcs)
    #TODO: I hate storing variable like this. It's probably a sign that I need to create a class and go more object oriented. *shrug* 
    object <- adjust.output$object
    accept <- adjust.output$accept
    accept.count <- accept.count + 1
  }
  
  idents <- as.int.factor(object@active.ident) # as.int.factor should probably be a class function
  levels(object@active.ident) <- 0:sum(1:max(idents) %in% unique(idents)) # ensures ident levels are sequential integers
  
  return(
    list(
      object = object,
      initial.eval = initial.eval.Test,
      final.eval = Test.method(object, method = "autoclustr"),
      accept.count = accept.count
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

#TODO: Remove this or rewrite it
getICVI <- function(object, npc = NULL) {
  icvi <- suppressWarnings(
    index.S(
      d = dist(
        Embeddings(object, reduction = "pca")[, 1:npc],
        method = "manhattan"),
      cl = as.integer(Idents(object))
    )
  )
  if(is.nan(icvi)) {icvi <- 0}
  return(icvi)
}

as.int.factor <- function(x) {as.integer(levels(x))[x]}

adjustIdents <- function(object, temp.object, n.pcs) {
  embed.dist <- Embeddings(object)[ , 1:n.pcs] %>% dist(method = "manhattan")
  
  new.idents <- as.int.factor(temp.object@active.ident) + 1 + 
    max(as.int.factor(object@active.ident))
  names(new.idents) <- Cells(temp.object)
  
  combined.new <- as.int.factor(object@active.ident)
  names(combined.new) <- Cells(object)
  combined.new[names(new.idents)] <- as.integer(new.idents)
  
  original.SI <- suppressWarnings(
    index.S(
      d = embed.dist,
      cl = as.integer(Idents(object))
    )
  )
  new.SI <-  suppressWarnings(
    index.S(
      d = embed.dist,
      cl = test.cl#combined.new
      ) # It's concerning that index.S is so fragile and our whole package depends on it...
    )
  test.cl <- integer()
  for( i in seq_len( length(combined.new) ) ){
    test.cl[i] <- as.integer(unname(combined.new)[i])
  }
  #if(is.nan(new.SI)) {new.SI <- 0} ; if(is.nan(original.SI)) {original.SI <- 0}
  # This is, like, not a great form of error checking, but here we are
  if(is.nan(new.SI) || is.nan(original.SI)) {
    embed.pca <- Embeddings(object)[ , 1:n.pcs]
    original.SI <- intCriteria(embed.pca, as.int.factor(object@active.ident), crit = "sil")[[1]]
    new.SI <- intCriteria(embed.pca, combined.new, crit = "sil")[[1]]
    index.s.fails <<- index.s.fails + 1
  }
  if(original.SI > new.SI) {
    print(paste(object@project.name, "REJECTED sub clustering"))
    print(paste("Original:", round(original.SI, 4), "New:", round(new.SI, 4)))
    accept <- FALSE
  } else if (original.SI <= new.SI) {
    object@active.ident <- as.factor(combined.new)
    print(paste(object@project.name, "ACCEPTED sub clustering"))
    print(paste("Original:", round(original.SI, 4), "New:", round(new.SI, 4)))
    accept <- TRUE
  }
  return(list(object=object, accept = accept))
}

Proc.data <- function(data, algorithm = "seurat", expr.meas = "umi", n.features = 2000, seed = NULL){
  
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

source("RScripts/nPC.funs.R")

#Should be sourcing stuff here, but
# TODO: Add verbosity option to AutoClustR
b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv")
b1 <- Prep.data(b1)
b1 <- Proc.data(b1)

b2 <- Baron.to.SCexp("Data/Baron-2/GSM2230758_human2_umifm_counts.csv")
b2 <- Prep.data(b2)
b2 <- Proc.data(b2)

resTest <- AutoClustR.flow(data = b1, expr.meas = "umi")

kolodz <- Kolodz.to.SCexp("Data/Kolodz/counttable_es.csv")
kolodz <- Prep.data(kolodz)
kolodz <- Proc.data(kolodz, expr.meas = "cpm")

ac.runtime <- Sys.time()
acResults <- AutoClustR(b2,
                     # file.path = "",
                     method = "Bayesian",
                     n.priors = 16,
                     n.starts = 16)
ac.runtime <- Sys.time() - ac.runtime
