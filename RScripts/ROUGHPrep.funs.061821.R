Remove.ERCCs <- function(data){
  
  # library(SingleCellExperiment)
  
  features <- grep(pattern = "^ERCC-", rownames(data), ignore.case = TRUE, invert = TRUE)
  data <- data[features, ]
  
  return(value = data)
}

Filter.cells <- function(data, outlier.thresh = 3){
  
  # library(SingleCellExperiment)
  
  n.features <- colSums(counts(data) > 0)
  mito.genes <- grep(pattern = "^MT-", rownames(data), ignore.case = TRUE)
  pct.mito <- 100 * colSums(counts(data)[mito.genes, ]) / colSums(counts(data))
  
  colData(data)$n.features <- n.features
  colData(data)$pct.mito <- pct.mito
  
  log.features <- log1p(n.features)
  log.mito <- log1p(pct.mito)
  
  min.features <- exp(median(log.features) - outlier.thresh * mad(log.features)) - 1
  max.features <- exp(median(log.features) + outlier.thresh * mad(log.features)) - 1
  max.mito <- exp(median(log.mito) + outlier.thresh * mad(log.mito)) - 1
  
  data <- data[, n.features >= min.features & n.features <= max.features & pct.mito <= max.mito]
  
  return(value = data)
}

Remove.zeroes <- function(data){
  
  # library(SingleCellExperiment)
  
  n.cells <- rowSums(counts(data) > 0)
  rowData(data)$n.cells <- n.cells
  data <- data[n.cells >= 1, ]
  
  return(value = data)
  
}

Prep.data <- function(data) {
  
  data <- Remove.ERCCs(data) %>%
    Filter.cells() %>%
    Remove.zeroes()
  
  return(data)
  
}
