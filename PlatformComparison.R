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
source("RScripts/Load.funs.0618.R")
source("RScripts/ROUGHPrep.funs.061821.R")
source("RScripts/Clust.funs.061821.R")
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
seeds <- runif(40)
for(i in 1:10) {
  Baron.1.results$CellTrails[ i, ] <- CellTrails.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$CIDR[ i, ] <- CIDR.flow(data = b1, expr.meas = expr.key$Baron)
  Baron.1.results$IKAP[ i, ] <- IKAP.flow(data = b1, expr.meas = expr.key$Baron,
                                          seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Baron.1.results$RaceID[ i, ] <- RaceID.flow(data = b1, expr.meas = expr.key$Baron,
                                              seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Baron.1.results$SC3[ i, ] <- SC3.flow(data = b1, expr.meas = expr.key$Baron,
                                        seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Baron.1.results$Seurat[ i, ] <- Seurat.flow(data = b1, expr.meas = expr.key$Baron,
                                              seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}
saveRDS(Baron.1.results, file = "Results/Baron1/Baron1.rds")
for(algo in names(Baron.1.results)) {
  write.csv(Baron.1.results[[algo]], file = paste0("Results/Baron1/Baron1_", algo, ".csv" ))
}

# Goolam ----
goolam <- Goolam.to.SCexp("Data/Goolam/Goolam_et_al_2015_count_table.tsv")
goolam <- Prep.data(goolam)
Goolam.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Goolam.results$CellTrails[ i, ] <- CellTrails.flow(data = goolam, expr.meas = expr.key$Goolam)
  Goolam.results$CIDR[ i, ] <- CIDR.flow(data = goolam, expr.meas = expr.key$Goolam)
  Goolam.results$IKAP[ i, ] <- IKAP.flow(data = goolam, expr.meas = expr.key$Goolam,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Goolam.results$RaceID[ i, ] <- RaceID.flow(data = goolam, expr.meas = expr.key$Goolam,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Goolam.results$SC3[ i, ] <- SC3.flow(data = goolam, expr.meas = expr.key$Goolam,
                                       seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Goolam.results$Seurat[ i, ] <- Seurat.flow(data = goolam, expr.meas = expr.key$Goolam,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}

saveRDS(Goolam.results, file = "Results/Goolam/Goolam.rds")
for(algo in names(Goolam.results)) {
  write.csv(Goolam.results[[algo]], file = paste0("Results/Goolam/Goolam_", algo, ".csv" ))
}

# Kolodz ----
Kolodz <- Kolodz.to.SCexp("Data/Kolodz/counttable_es.csv") %>%
  Prep.data()
Kolodz.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Kolodz.results$CellTrails[ i, ] <- CellTrails.flow(data = Kolodz, expr.meas = expr.key$Kolodz)
  Kolodz.results$CIDR[ i, ] <- CIDR.flow(data = Kolodz, expr.meas = expr.key$Kolodz)
  Kolodz.results$IKAP[ i, ] <- IKAP.flow(data = Kolodz, expr.meas = expr.key$Kolodz,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Kolodz.results$RaceID[ i, ] <- RaceID.flow(data = Kolodz, expr.meas = expr.key$Kolodz,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Kolodz.results$SC3[ i, ] <- SC3.flow(data = Kolodz, expr.meas = expr.key$Kolodz,
                                       seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Kolodz.results$Seurat[ i, ] <- Seurat.flow(data = Kolodz, expr.meas = expr.key$Kolodz,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}
saveRDS(Kolodz.results, file = "Results/Kolodz/Kolodz.rds")
for(algo in names(Kolodz.results)) {
  write.csv(Kolodz.results[[algo]], file = paste0("Results/Kolodz/Kolodz_", algo, ".csv" ))
}
rm(Kolodz)

# Baron 2 ----
b2 <- Baron.to.SCexp("Data/Baron-2/GSM2230758_human2_umifm_counts.csv") %>%
  Prep.data()
Baron.2.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Baron.2.results$CellTrails[ i, ] <- CellTrails.flow(data = b2, expr.meas = expr.key$Baron)
  Baron.2.results$CIDR[ i, ] <- CIDR.flow(data = b2, expr.meas = expr.key$Baron)
  Baron.2.results$IKAP[ i, ] <- IKAP.flow(data = b2, expr.meas = expr.key$Baron,
                                          seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Baron.2.results$RaceID[ i, ] <- RaceID.flow(data = b2, expr.meas = expr.key$Baron,
                                              seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Baron.2.results$SC3[ i, ] <- SC3.flow(data = b2, expr.meas = expr.key$Baron,
                                        seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Baron.2.results$Seurat[ i, ] <- Seurat.flow(data = b2, expr.meas = expr.key$Baron,
                                              seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}
saveRDS(Baron.2.results, file = "Results/Baron2/Baron2.rds")
for(algo in names(Baron.2.results)) {
  write.csv(Baron.2.results[[algo]], file = paste0("Results/Baron2/Baron2_", algo, ".csv" ))
}
rm(b2)

# Loh ----
Loh <- Loh.to.SCexp("Data/Loh/GSM2257302_All_samples_sc_tpm.txt") %>%
  Prep.data()
Loh.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Loh.results$CellTrails[ i, ] <- CellTrails.flow(data = Loh, expr.meas = expr.key$Loh)
  Loh.results$CIDR[ i, ] <- CIDR.flow(data = Loh, expr.meas = expr.key$Loh)
  Loh.results$IKAP[ i, ] <- IKAP.flow(data = Loh, expr.meas = expr.key$Loh,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Loh.results$RaceID[ i, ] <- RaceID.flow(data = Loh, expr.meas = expr.key$Loh,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Loh.results$SC3[ i, ] <- SC3.flow(data = Loh, expr.meas = expr.key$Loh,
                                       seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Loh.results$Seurat[ i, ] <- Seurat.flow(data = Loh, expr.meas = expr.key$Loh,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}
saveRDS(Loh.results, file = "Results/Loh/Loh.rds")
for(algo in names(Loh.results)) {
  write.csv(Loh.results[[algo]], file = paste0("Results/Loh/Loh_", algo, ".csv" ))
}
rm(Loh)

# Menon PR Combined ---- 
MenonPR <- Menon.to.SCexp(data.file.path = "Data/Menon-P/GSE137537_counts.mtx",
                          feat.file.path = "Data/Menon-P/GSE137537_gene_names.txt",
                          cell.file.path = "Data/Menon-P/GSE137537_sample_annotations.tsv",
                          sample = c("PR", "PR2", "PR3")) %>%
  Prep.data()
MenonPR.results <- res.template
seeds <- runif(40)
for(i in 1) {
  MenonPR.results$CellTrails[ i, ] <- CellTrails.flow(data = MenonPR, expr.meas = expr.key$Menon)
  # MenonPR.results$CIDR[ i, ] <- CIDR.flow(data = MenonPR, expr.meas = expr.key$Menon)
  # MenonPR.results$IKAP[ i, ] <- IKAP.flow(data = MenonPR, expr.meas = expr.key$Menon,
  #                                        seed = (Sys.time() %>% as.numeric()) * seeds[i])
  # MenonPR.results$RaceID[ i, ] <- RaceID.flow(data = MenonPR, expr.meas = expr.key$Menon,
  #                                            seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  MenonPR.results$SC3[ i, ] <- SC3.flow(data = MenonPR, expr.meas = expr.key$Menon,
                                       seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  MenonPR.results$Seurat[ i, ] <- Seurat.flow(data = MenonPR, expr.meas = expr.key$Menon,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
  saveRDS(MenonPR.results, file = "Results/MenonPR/MenonPR.rds")
}
saveRDS(MenonPR.results, file = "Results/MenonPR/MenonPR.rds")
for(algo in names(MenonPR.results)) {
  write.csv(MenonPR.results[[algo]], file = paste0("Results/MenonPR/MenonPR_", algo, ".csv" ))
}
rm(MenonPR)

# Pollen ----
Pollen <- Pollen.to.SCexp("Data/Pollen/NBT_hiseq_linear_tpm_values.txt") %>%
  Prep.data()
Pollen.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Pollen.results$CellTrails[ i, ] <- CellTrails.flow(data = Pollen, expr.meas = expr.key$Pollen)
  Pollen.results$CIDR[ i, ] <- CIDR.flow(data = Pollen, expr.meas = expr.key$Pollen)
  Pollen.results$IKAP[ i, ] <- IKAP.flow(data = Pollen, expr.meas = expr.key$Pollen,
                                      seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Pollen.results$RaceID[ i, ] <- RaceID.flow(data = Pollen, expr.meas = expr.key$Pollen,
                                          seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Pollen.results$SC3[ i, ] <- SC3.flow(data = Pollen, expr.meas = expr.key$Pollen,
                                    seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Pollen.results$Seurat[ i, ] <- Seurat.flow(data = Pollen, expr.meas = expr.key$Pollen,
                                          seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}
saveRDS(Pollen.results, file = "Results/Pollen/Pollen.rds")
for(algo in names(Pollen.results)) {
  write.csv(Pollen.results[[algo]], file = paste0("Results/Pollen/Pollen_", algo, ".csv" ))
}
rm(Pollen)

# Ranum ----
Ranum <- Ranum.to.SCexp("Data/Ranum/GSE114157_p15_Expression_Matrix.csv") %>%
  Prep.data()
Ranum.results <- res.template
seeds <- runif(40)
for(i in 1:10) {
  Ranum.results$CellTrails[ i, ] <- CellTrails.flow(data = Ranum, expr.meas = expr.key$Ranum)
  Ranum.results$CIDR[ i, ] <- CIDR.flow(data = Ranum, expr.meas = expr.key$Ranum)
  Ranum.results$IKAP[ i, ] <- IKAP.flow(data = Ranum, expr.meas = expr.key$Ranum,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  Ranum.results$RaceID[ i, ] <- RaceID.flow(data = Ranum, expr.meas = expr.key$Ranum,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  Ranum.results$SC3[ i, ] <- SC3.flow(data = Ranum, expr.meas = expr.key$Ranum,
                                       seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  Ranum.results$Seurat[ i, ] <- Seurat.flow(data = Ranum, expr.meas = expr.key$Ranum,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
}
saveRDS(Ranum.results, file = "Results/Ranum/Ranum.rds")
for(algo in names(Ranum.results)) {
  write.csv(Ranum.results[[algo]], file = paste0("Results/Ranum/Ranum_", algo, ".csv" ))
}
rm(Ranum)