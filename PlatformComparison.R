#Comparison Workflow

#Load platforms
library(aricode)
library(cidr)
library(RaceID)
library(SC3)
library(Seurat)
library(SingleCellExperiment)

#General
library(tidyverse)
library(Matrix)
library(clusterSim)
library(clValid)
library(fpc)
library(Rtsne)
library(data.table)

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



