Baron.to.SCexp <- function(file.path){

# library(SingleCellExperiment)

  data <- read.csv(file = file.path)
  rownames(data) <- data[, 1]
  labels <- data[, 3]
  data <- t(as.matrix(data[, -(1:3)]))
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Goolam.to.SCexp <- function(file.path){

# library(SingleCellExperiment)
# library(tidyverse)

  data <- read.delim(file = file.path)
  colnames(data) <- str_replace(string = colnames(data), pattern = "X", replacement = "X_")
  data <- as.matrix(data)
  labels <- sapply(str_split(string = colnames(data), pattern = "_"), "[", 2)
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Kolodz.to.SCexp <- function(file.path){

# library(SingleCellExperiment)
# library(tidyverse)

  data <- read.delim(file = file.path, sep = " ")
  data <- as.matrix(data)
  labels <- sapply(str_split(string = colnames(data), pattern = "_"), "[", 3)
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Loh.to.SCexp <- function(file.path){

# library(SingleCellExperiment)
# library(tidyverse)

  data <- read.delim(file = file.path)
  rownames(data) <- data[, 2]
  data <- as.matrix(data[, -c(1, 2)])
  labels <- sapply(str_split(string = colnames(data), pattern = "\\."), "[", 1)
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Menon.to.SCexp <- function(data.file.path, feat.file.path, cell.file.path, sample = c("MR", "PR", "MR2", "PR2", "MR3", "PR3")){

# library(Matrix)
# library(SingleCellExperiment)

  data <- readMM(file = data.file.path)
  features <- read.delim(file = feat.file.path, header = FALSE)
  cells <- read.delim(file = cell.file.path)
  dimnames(data) <- list(features[, 1], cells[, 1])
  data <- as.matrix(data[, cells[, 2] == sample])
  labels <- cells[cells[, 2] == sample, 49]
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Pollen.to.SCexp <- function(file.path){

# library(SingleCellExperiment)
# library(tidyverse)

  data <- read.delim(file = file.path)
  data <- as.matrix(data)
  labels <- sapply(str_split(string = colnames(data), pattern = "_"), "[", 2)
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Ranum.to.SCexp <- function(file.path){

# library(SingleCellExperiment)
# library(tidyverse)

  data <- read.csv(file = file.path)

  data[21806, 1] <- paste0(data[21806, 1], "-1")
  data[44664, 1] <- paste0(data[44664, 1], "-2")
  data[44666, 1] <- paste0(data[44666, 1], "-1")
  data[46464, 1] <- paste0(data[46464, 1], "-2")

  rownames(data) <- data[, 1]
  data <- as.matrix(data[, -1])
  labels <- sapply(str_split(string = colnames(data), pattern = "_"), "[", 1)
  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}

Zeisel.to.SCexp <- function(mrna.file.path, mito.file.path){

# library(SingleCellExperiment)

  mrna.data <- read.delim(file = mrna.file.path)
  mrna.data <- mrna.data[, -2]
  mrna.data <- mrna.data[, order(mrna.data[7, ])]
  features <- mrna.data[-(1:10), 1]
  cells <- mrna.data[7, -1]
  labels <- as.character(mrna.data[1, -1])
  mrna.data <- as.matrix(mrna.data[-(1:10), -1])

  mito.data <- read.delim(file = mito.file.path)
  mito.data <- mito.data[, -2]
  mito.data <- mito.data[, order(mito.data[7, ])]
  features <- append(features, mito.data[-(1:10), 1])
  mito.data <- as.matrix(mito.data[-(1:10), -1])
  
  data <- rbind(mrna.data, mito.data)
  data <- apply(data, 2, as.numeric)
  dimnames(data) <- list(features, cells)

  data <- SingleCellExperiment(assays = list(counts = data), colData = list(labels = labels))

  return(value = data)
}