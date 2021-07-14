#ICVI Datasets
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)


#%%%%%%%%%%%%%%%%%%%%%%   Pancreatic/Baron Datasets %%%%%%%%%%%%%%%%%%%%%%%%%%%%
prepBaron <- function(file.path){
  panc <- read.csv(file.path)
  panc <- t(panc)
  colnames(panc) <- panc[2, ]
  panc.annotation <- panc[3,]
  panc <- panc[-c(1,2,3), ]
  panc.seurat <- CreateSeuratObject(panc, project = "HumanPancreas2")
  panc.seurat$annotation <- panc.annotation
  panc <- panc.seurat
  log.features <- log(panc$nFeature_RNA + 1)
  min.features <- exp(mean(log.features) - 3 * sd(log.features)) - 1
  max.features <- exp(mean(log.features) + 3 * sd(log.features)) - 1
  panc <- subset(panc, subset = nFeature_RNA > min.features & nFeature_RNA < max.features)
  panc <- SCTransform(panc.seurat)
  panc <- RunPCA(panc)
  panc@active.ident <- as.factor(panc$annotation)
  panc$na <- !is.na(panc$annotation)
  panc <- subset(panc, na)
  return(panc)
}

panc.1 <- prepBaron("C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet1/GSM2230757_human1_umifm_counts.csv/GSM2230757_human1_umifm_counts.csv")
panc.2 <- prepBaron("C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet2/GSM2230758_human2_umifm_counts.csv/GSM2230758_human2_umifm_counts.csv")


#%%%%%%%%%%%%%%%%%   Mouse Brain/Linnarsson Datasets  %%%%%%%%%%%%%%%%%%%%%%%%%%
prepZeisel <- function(file.path){
  mb <- as.data.frame(fread(file.path)) 
  mb.counts <- mb[-(1:11), -2]
  rownames(mb.counts) <- mb.counts[[1]]
  mb.counts <- mb.counts[-1]
  mb.meta <- mb[1:10, -1]
  names(mb.counts) <- mb.meta[8, -1]
  mb.seurat <- CreateSeuratObject(mb.counts, project = "Zeisel 2015")
  mb.seurat$annotation <- as.factor(mb.meta[2, -1])
  log.features <- log(mb.seurat$nFeature_RNA + 1)
  min.features <- exp(mean(log.features) - 3 * sd(log.features)) - 1
  max.features <- exp(mean(log.features) + 3 * sd(log.features)) - 1
  panc <- subset(mb.seurat, subset = nFeature_RNA > min.features & nFeature_RNA < max.features)
  mbs <- SCTransform(mb.seurat)
  mbs <- RunPCA(mbs)
  mbs@active.ident <- as.factor(mbs$annotation)
  return(mbs)
}

zeisel <- prepZeisel("C:/Users/asolivai/Desktop/R_Files/AutoClustR/MouseBrainCells/expression_mRNA_17-Aug-2014.txt")




#%%%%%%%%%%%%%%%%%%%%%%%%%%   Menon/Retinal Datasets %%%%%%%%%%%%%%%%%%%%%%%%%%%
filterMenon <- function(object)  {
  log.features <- log(object$nFeature_RNA + 1)
  min.features <- exp(mean(log.features) - 3 * sd(log.features)) - 1
  max.features <- exp(mean(log.features) + 3 * sd(log.features)) - 1
  object <- subset(object, subset = nFeature_RNA > min.features & nFeature_RNA < max.features)
  object <- SCTransform(object)
  object <- RunPCA(object)
  object@active.ident <- as.factor(object$Labels)
  return(object)
}

counts <- readMM("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/outs/GSE137537_counts.mtx")
features <- fread("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/outs/GSE137537_gene_names.txt", header = F) %>%
  as.data.frame()
annotation <- fread("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/outs/GSE137537_sample_annotations.tsv") %>%
  as.data.frame()

dimnames(counts) <- list(features[[1]], annotation$Barcode)
rownames(annotation) <- annotation[[1]]

menon <- CreateSeuratObject(counts = counts, project = "Menon", meta.data = annotation[-1])
table(menon$Labels)

#Got it created, now just have to subset PR and MR
pr <- subset(menon, tissue == "PR") %>%
  filterMenon()
mr <- subset(menon, tissue == "MR") %>%
  filterMenon()

pr$annotation <- pr$Labels
mr$annotation <- mr$Labels

rm(annotation, counts, features, menon)

datasets <- list(Panc1 = panc.1,
                 Panc2 = panc.2,
                 Zeisel = zeisel,
                 PR = pr,
                 MR = mr)











