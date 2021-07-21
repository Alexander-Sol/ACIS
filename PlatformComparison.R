#Comparison Workflow

#Load platforms
library(aricode)
library(cidr)
library(RaceID)
library(SC3)
library(Seurat)
library(SingleCellExperiment)
library(IKAP)

#General
library(tidyverse)
library(Matrix)
library(clusterSim)
library(clValid)
library(fpc)
library(Rtsne)
library(data.table)
library(stringr)
library(ggplot2)
library(RColorBrewer)

#AutoClustR Specific
library(lhs)
library(lme4) # Nelder_Mead
library(ParBayesianOptimization)
library(stringi)


# Seurat Pilot

b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv")
b1 <- Prep.data(b1)

b1 <- Proc.data(b1, algorithm = "seurat", expr.meas = "umi") 

b1 <- JackStraw(b1, num.replicate = 100, dims = 35) %>%
  ScoreJackStraw(dims = 1:35)
npc.num <- min(which(b1@reductions$pca@jackstraw$overall.p.values[ , 2] > 0.05), na.rm = T) - 1
if(npc.num <= 0) { npcs <- 35 }

SI.Man <- ClustR(object = b1, n.pcs = npc.num)

# CellFindR Pilot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_loc <- "CellFindR_Output/"
proj_name <- "B1_Test"
setwd(file_loc)

res <- find_res(b1)
b1 <- FindNeighbors(b1, dims = 1:20, k.param = 20)
b1 <- FindClusters(b1, resolution = 0.6)

b1_CFR <- sub_clustering(b1, output_folder = "CellFindR_Output/", proj_name = proj_name) #CellFindR is going to take the longest

# CIDR Pilot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv") %>%
  Prep.data() %>%
  Proc.data(algorithm = "cidr", expr.meas = "umi") 
cidr.npcs <- calc_npc(b1@variation)

results <- data.frame(matrix(data = 0, nrow = 10, ncol = 3))
names(results) <- c("ARI", "Clusters", "Runtime")


for( i in 1:10 ) {
  runtime <- Sys.time()
  cidr.res <- ClustR(object = b1, n.pcs = cidr.npcs)
  runtime <- Sys.time() - runtime
  results[i, 1] <- cidr.res$ARI[[2]]
  results[i, 2] <- cidr.res$ARI[[1]]
  results[i, 3] <- runtime
}


# IKAP Pilot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv") %>%
  Prep.data() %>%
  Proc.data(algorithm = "IKAP", expr.meas = "umi")

runtime <- Sys.time()
I.Res <- IKAP(b1, out.dir = "IKAP_Output")
runtime <- Sys.time() - runtime






