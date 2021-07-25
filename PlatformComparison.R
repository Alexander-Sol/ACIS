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
source("RScripts/RoughTest.funs.R")
source("RScripts/Workflows.R")

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

# Generate results template -----
results.table <- data.frame(matrix(data = 0, nrow = 10, ncol = 4))
names(results.table) <- c("ARI", "n.clusters", "runtime", "seed")

res.template <- list(
  AutoClustR = results.table,
  CellTrails = results.table,
  CIDR = results.table,
  IKAP = results.table,
  RaceID = results.table,
  SC3 = results.table,
  Seurat = results.table
)

# Baron 1 ----
b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv")
b1 <- Prep.data(b1)
Baron.1.results <- res.template
for(i in 1:10) {
  Baron.1.results$CellTrails[ i, ] <- CellTrails.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$CIDR[ i, ] <- CIDR.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$IKAP[ i, ] <- IKAP.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$RaceID[ i, ] <- RaceID.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$SC3[ i, ] <- SC3.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$Seurat[ i, ] <- Seurat.flow(data = b1, expr.meas = expr.key$Baron)
}
saveRDS(Baron.1.results, file = "Results/Baron1/Baron1.rds")
for(algo in names(Baron.1.results)) {
  write.csv(Baron.1.results[[algo]], file = paste0("Results/Baron1/Baron1_", algo, ".csv" ))
}

# Goolam ----
goolam <- Goolam.to.SCexp("Data/Goolam/Goolam_et_al_2015_count_table.tsv")
Goolam.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Goolam.results$CellTrails[ i, ] <- CellTrails.flow(data = goolam, expr.meas = expr.key$goolam)
  Goolam.results$CIDR[ i, ] <- CIDR.flow(data = goolam, expr.meas = expr.key$goolam)
  Goolam.results$IKAP[ i, ] <- IKAP.flow(data = goolam, expr.meas = expr.key$goolam,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Goolam.results$RaceID[ i, ] <- RaceID.flow(data = goolam, expr.meas = expr.key$goolam,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Goolam.results$SC3[ i, ] <- SC3.flow(data = goolam, expr.meas = expr.key$goolam,
                                       seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Goolam.results$Seurat[ i, ] <- Seurat.flow(data = goolam, expr.meas = expr.key$goolam,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}

saveRDS(Goolam.results, file = "Results/Goolam/Goolam.rds")
for(algo in names(Goolam.results)) {
  write.csv(Goolam.results[[algo]], file = paste0("Results/Goolam/Goolam_", algo, ".csv" ))
}









