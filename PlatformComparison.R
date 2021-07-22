#Comparison Workflow

#Load platforms
library(aricode)
library(cidr)
library(RaceID)
library(SC3)
library(Seurat)
library(SingleCellExperiment)
source("RScripts/IKAP.R")
library(CellTrails)
library(scran)

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
library(wrapr)

#AutoClustR Specific
library(lhs)
library(lme4) # Nelder_Mead
library(ParBayesianOptimization)
library(stringi)

#Comparison Scripts
source("RScripts/ROUGHPrep.funs.061821.R")
source("RScripts/Clust.funs.061821.R")
source("RScripts/Load.funs.0618.R")

#Expression Key
expr.key <- list(
  Baron = "umi",
  Goolam = "cpm",
  Kolodz = "cpm",
  Loh = "tpm",
  Menon = "umi",
  Pollen = "tpm",
  Ranum = "cpm",
  Zeisel = "umi"
)


results.table <- data.frame(matrix(data = 0, nrow = 10, ncol = 4))
names(results.table) <- c("ARI", "n.clusters", "runtime", "seed")
Baron.1.results <- list(
  AutoClustR = results.table,
  CellTrails = results.table,
  CIDR = results.table,
  IKAP = results.table,
  RaceID = results.table,
  SC3 = results.table,
  Seurat = results.table
)



b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv")
b1 <- Prep.data(b1)

for(i in 1) {
  Baron.1.results$CellTrails[ i, ] <- CellTrails.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$CIDR[ i, ] <- CIDR.flow(data = b1, expr.meas = expr.key$Baron)
  Baron1.results$IKAP[ i, ] <- IKAP.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$RaceID[ i, ] <- CIDR.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$SC3[ i, ] <- CIDR.flow(data = b1, expr.meas = expr.key$Baron)
  Baron1.results$Seurat[ i, ] <- Seurat.flow(data = b1, expr.meas = expr.key$Baron)
}

saveRDS(Baron.1.results, file = "Results/Baron1/Baron1.rds")
for(algo in names(Baron.1.results)) {
  write.csv(Baron.1.results[[algo]], file = paste0("Results/Baron1/Baron1_", algo, ".csv" ))
}

CellTrails.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.CellTrails(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  eval <- Test.method(results$object, method = "celltrails")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

CIDR.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.CIDR(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  eval <- Test.method(results, method = "cidr")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      NA
    )
  )
  
}

IKAP.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.IKAP(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  eval <- Test.method(results$object, method = "ikap")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

RaceID.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.RaceID(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  eval <- Test.method(results$object, method = "raceid")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

SC3.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.SC3(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  eval <- Test.method(results$object, method = "sc3")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

Seurat.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.Seurat(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  eval <- Test.method(results$object, method = "seurat")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}









