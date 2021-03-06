Run.CellTrails <- function(data, expr.meas = "umi") {
  
  # if(is.null(seed)){
  #   seed <- round(1e6 * runif(n = 1))
  # }
  # 
  rowData(data)$feature_symbol <- rownames(data)
  
  if(expr.meas == "umi"){
    logcounts(data) <- log1p(1e6 * proportions(counts(data), margin = 2))
  }else if(expr.meas == "tpm" | expr.meas == "cpm"){
    logcounts(data) <- log1p(counts(data))
  }else if(expr.meas == "mnn"){
    logcounts(data) <- assay(data, "corrected")
  }
  
  trajFeatureNames(data) <- tfeat <- filterTrajFeaturesByFF(data, threshold = 1.7, show_plot = T)
  data.sub <- data[tfeat, ]
  # Manual check and subset for zero count cells after filtering
  data <- data[ , colSums(data.sub@assays@data$counts) != 0]
  data.sub <- data.sub[ , colSums(data.sub@assays@data$counts) != 0]
  se <- embedSamples(data.sub)
  d <- findSpectrum(se$eigenvalues, frac = 100)
  latentSpace(data) <- se$components[ , d]
  states(data) <- clusters <- findStates(data, min_size=0.01, min_feat=5, max_pval=1e-4, min_fc=2)
  
  return(value = list(object = data, seed = NA))
}

Run.CIDR <- function(data, expr.meas = "umi"){

  # library(cidr)
  # library(SingleCellExperiment)
  # library(tidyverse)

  if(expr.meas == "umi") {
    expr.meas <- "raw"
  } else if(expr.meas == "mnn") {
    unnormed.data <- assay(data, "corrected") %>% expm1() %>% '/'(1e6) 
    counts(data) <- apply(unnormed.data[1:1000, 1:10],
                          MARGIN = 2,
                          FUN = function(x) { 
                            if( sum(x>0) ) {
                              x / min( x[ x>0 ] ) 
                            } else {
                              x / 1
                            }  
                          }) #This reverses 'proportions', probably
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
  
  print("CIDR Ran")

  return(object)
}

Run.IKAP <- function(data, file.path = "IKAP/", expr.meas = "umi", seed = NULL){

#   library(Seurat)
#   library(SingleCellExperiment)
#   library(tidyverse)

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
  print("IKAP Ran")
  return(value = list(object = object, seed = seed))
}

Run.RaceID <- function(data, seed = NULL){

  # library(RaceID)
  # library(SingleCellExperiment)
  # library(tidyverse)

  if(is.null(seed)){
    seed <- round(1e6 * runif(n = 1))
  }

  object <- SCseq(expdata = counts(data))

  object <- filterdata(object = object, mintotal = 1) %>%
    compdist() %>%
    clustexp(rseed = seed) %>%
    findoutliers()

  object@cluster$labels <- colData(data)$labels
  print("RACEID Ran")

  return(value = list(object = object, seed = seed))
}

Run.SC3 <- function(data, expr.meas = "umi", seed = NULL){

  # library(SC3)
  # library(SingleCellExperiment)
  # library(tidyverse)

  if(is.null(seed)){
    seed <- round(1e6 * runif(n = 1))
  }

  rowData(data)$feature_symbol <- rownames(data)

  if(expr.meas == "umi"){
    logcounts(data) <- log1p(1e6 * proportions(counts(data), margin = 2))
  }else if(expr.meas == "tpm" | expr.meas == "cpm"){
    logcounts(data) <- log1p(counts(data))
  } else if(expr.meas == "mnn") {
    counts(data)[counts(data) == 1] <- 0 # Converting from corrected to counts should be done via expm1, but because 
                                        # corrected contains 0 < values < 1, subtracting one causes errors down the line
    logcounts(data) <- assay(data, "corrected")
  }

  object <- sc3_prepare(object = data, rand_seed = seed) %>%
    sc3_estimate_k() %>%
    sc3_calc_dists() %>%
    sc3_calc_transfs()

  object <- sc3_kmeans(object = object, ks = metadata(object)$sc3$k_estimation) %>%
    sc3_calc_consens()

  print("SC3 Ran")

  if(dim(data)[[2]] > 5000){
    object <- sc3_run_svm(object = object, ks = metadata(object)$sc3$k_estimation)
  }

  return(value = list(object = object,
                      seed = seed))
}

Run.Seurat <- function(data, expr.meas = "umi", seed = NULL){

  # library(Seurat)
  # library(SingleCellExperiment)
  # library(tidyverse)

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

  return(value = list(object = object, seed = seed))
}