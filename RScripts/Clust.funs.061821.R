Run.CellFindR <- function(data, file.path, expr.meas = "umi"){

  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)

  # Note: Set RunUMAP(verbose = FALSE) in CellFindR source code
  # Note: For Seurat v4, set Find[All]Markers(base = exp(1)) in CellFindR source code
  object <- CreateSeuratObject(counts = counts(data))
  object@meta.data$labels <- colData(data)$labels
  object@meta.data$n.features <- colData(data)$n.features
  object@meta.data$pct.mito <- colData(data)$pct.mito

  if(expr.meas == "umi"){

    object <- NormalizeData(object = object, scale.factor = 1e6) %>%
      FindVariableFeatures() %>%
      ScaleData(vars.to.regress = c("pct.mito", "n.features")) 

  }else if(expr.meas == "tpm" | expr.meas == "cpm"){

    # Note: For expr.meas = "tpm" or "cpm", set FindVariableFeatures(selection.method = "mvp") in CellFindR source code 
    object@assays$RNA <- CreateAssayObject(data = log1p(counts(data)))
    object <- FindVariableFeatures(object = object, selection.method = "mvp") %>%
      ScaleData(vars.to.regress = c("pct.mito", "n.features"))

  }
  
  object <- RunPCA(object = object) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.1) %>%
    RunUMAP(dims = 1:20)

  resolution <- find_res(tenx = object)

  object <- FindClusters(object = object, resolution = resolution) %>%
    sub_clustering(output_folder = file.path)

  return(value = object)
}

Run.CIDR <- function(data, expr.meas = "umi"){

  library(cidr)
  library(SingleCellExperiment)
  library(tidyverse)

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
    scPCA(plotPC = FALSE) %>%
    nPC()

  object <- scCluster(object = object, nPC = object@nPC)

  return(value = object)
}

Run.IKAP <- function(data, file.path, expr.meas = "umi", seed = NULL){

  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)

  if(is.null(seed)){
    seed <- round(1e6 * runif(n = 1))
  }
  
  # Note: For Seurat v4, change "avg_logFC" to "avg_log2FC" in IKAP source code
  object <- CreateSeuratObject(counts = counts(data))
  object@meta.data$labels <- colData(data)$labels
  object@meta.data$n.features <- colData(data)$n.features
  object@meta.data$pct.mito <- colData(data)$pct.mito

  if(expr.meas == "umi"){

    object <- NormalizeData(object = object, scale.factor = 1e6)
    object <- IKAP(
      sobj = object, 
      out.dir = file.path, 
      confounders = c("pct.mito", "n.features"), 
      plot.decision.tree = FALSE, 
      random.seed = seed
      )

  }else if(expr.meas == "tpm" | expr.meas == "cpm"){

    object@assays$RNA <- CreateAssayObject(data = log1p(counts(data)))
    # Note: Set FindVariableFeatures(selection.method = "mvp") in IKAP source code
    object <- IKAP(
      sobj = object, 
      out.dir = file.path, 
      confounders = c("pct.mito", "n.features"), 
      plot.decision.tree = FALSE, 
      random.seed = seed
      )

  }

  return(value = list(object, seed))
}

Run.RaceID <- function(data, seed = NULL){

  library(RaceID)
  library(SingleCellExperiment)
  library(tidyverse)

  if(is.null(seed)){
    seed <- round(1e6 * runif(n = 1))
  }

  object <- SCseq(expdata = counts(data))

  object <- filterdata(object = object, mintotal = 1) %>%
    compdist() %>%
    clustexp(rseed = seed) %>%
    findoutliers()

  object@cluster$labels <- colData(data)$labels

  return(value = list(object, seed))
}

Run.SC3 <- function(data, expr.meas = "umi", seed = NULL){

  library(SC3)
  library(SingleCellExperiment)
  library(tidyverse)

  if(is.null(seed)){
    seed <- round(1e6 * runif(n = 1))
  }

  rowData(data)$feature_symbol <- rownames(data)
  
  if(expr.meas == "umi"){
    logcounts(data) <- log1p(1e6 * proportions(counts(data), margin = 2))
  }else if(expr.meas == "tpm" | expr.meas == "cpm"){
    logcounts(data) <- log1p(counts(data))
  }

  object <- sc3_prepare(object = data, rand_seed = seed) %>%
    sc3_estimate_k() %>%
    sc3_calc_dists() %>%
    sc3_calc_transfs()

  object <- sc3_kmeans(object = object, ks = metadata(object)$sc3$k_estimation) %>%
    sc3_calc_consens()

  if(dim(data)[[2]] > 5000){
    object <- sc3_run_svm(object = object, ks = metadata(object)$sc3$k_estimation)
  }

  return(value = list(object, seed))
}

Run.Seurat <- function(data, expr.meas = "umi", seed = NULL){

  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)

  if(is.null(seed)){
    seed <- round(1e6 * runif(n = 1))
  }
  
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

  object <- RunPCA(object = object, seed.use = seed) %>%
    FindNeighbors() %>%
    FindClusters(random.seed = seed)

  return(value = list(object, seed))
}