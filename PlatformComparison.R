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
library(batchelor)

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
source("RScripts/nPC.funs.R")
source("RScripts/AutoClustR.funs.071321.R")

#Expression Key
expr.key <- list(
  Baron = "umi",
  Goolam = "cpm",
  Kolodz = "cpm",
  Loh = "tpm",
  Menon = "umi",
  MenonCorrected = "mnn",
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

seeds <- runif(40)

# Ad Hoc AutoClustR ----

#Baron 1
b1 <- Baron.to.SCexp("Data/Baron-1/GSM2230757_human1_umifm_counts.csv")
b1 <- Prep.data(b1)
AutoClustR.results.B1 <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.B1[[i]]  <- AutoClustR.flow(b1, expr.meas = expr.key$Baron, seed = seed)
  AutoClustR.results.B1[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.B1, file = "Results/AutoClustR/Baron1.rds")
rm(b1, AutoClustR.results.B1)

# Baron 2
b2 <- Baron.to.SCexp("Data/Baron-2/GSM2230758_human2_umifm_counts.csv") %>%
  Prep.data()
AutoClustR.results.b2 <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.b2[[i]]  <- AutoClustR.flow(b2, expr.meas = expr.key$Baron, seed = seed)
  AutoClustR.results.b2[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.b2, file = "Results/AutoClustR/Baron2.rds")
rm(b2, AutoClustR.results.b2)

# Goolam
goolam <- Goolam.to.SCexp("Data/Goolam/Goolam_et_al_2015_count_table.tsv") %>%
  Prep.data()
AutoClustR.results.goolam <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.goolam[[i]]  <- AutoClustR.flow(goolam, expr.meas = expr.key$Goolam, seed = seed)
  AutoClustR.results.goolam[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.goolam, file = "Results/AutoClustR/Goolam.rds")
rm(goolam, AutoClustR.results.goolam)

# Kolodz
kolodz <- Kolodz.to.SCexp("Data/Kolodz/counttable_es.csv") %>%
  Prep.data()
AutoClustR.results.kolodz <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.kolodz[[i]]  <- AutoClustR.flow(kolodz, expr.meas = expr.key$Kolodz, seed = seed)
  AutoClustR.results.kolodz[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.kolodz, file = "Results/AutoClustR/Kolodz.rds")
rm(kolodz, AutoClustR.results.kolodz)

# Loh
loh <- Loh.to.SCexp("Data/Loh/GSM2257302_All_samples_sc_tpm.txt") %>%
  Prep.data()
AutoClustR.results.loh <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.loh[[i]]  <- AutoClustR.flow(loh, expr.meas = expr.key$Loh, seed = seed)
  AutoClustR.results.loh[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.loh, file = "Results/AutoClustR/Loh.rds")
rm(loh, AutoClustR.results.loh)

# MenonPR
MenonPR <- Menon.to.SCexp(data.file.path = "Data/Menon-P/GSE137537_counts.mtx",
                                 feat.file.path = "Data/Menon-P/GSE137537_gene_names.txt",
                                 cell.file.path = "Data/Menon-P/GSE137537_sample_annotations.tsv",
                                 sample = c("PR", "PR2", "PR3")) %>%
  Prep.data()
AutoClustR.results.MenonPR <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.MenonPR[[i]]  <- AutoClustR.flow(MenonPR, expr.meas = expr.key$MenonPR, seed = seed)
  AutoClustR.results.MenonPR[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.MenonPR, file = "Results/AutoClustR/MenonPR.rds")
rm(MenonPR, AutoClustR.results.MenonPR)

# MenonMR
MenonMR <- Menon.to.SCexp(data.file.path = "Data/Menon-P/GSE137537_counts.mtx",
                          feat.file.path = "Data/Menon-P/GSE137537_gene_names.txt",
                          cell.file.path = "Data/Menon-P/GSE137537_sample_annotations.tsv",
                          sample = c("MR", "MR2", "MR3")) %>%
  Prep.data()
AutoClustR.results.MenonMR <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.MenonMR[[i]]  <- AutoClustR.flow(MenonMR, expr.meas = expr.key$Menon, seed = seed)
  AutoClustR.results.MenonMR[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.MenonMR, file = "Results/AutoClustR/MenonMR.rds")
rm(MenonMR, AutoClustR.results.MenonMR)

# MenonMR
MenonMR <- Menon.to.SCexp(data.file.path = "Data/Menon-P/GSE137537_counts.mtx",
                          feat.file.path = "Data/Menon-P/GSE137537_gene_names.txt",
                          cell.file.path = "Data/Menon-P/GSE137537_sample_annotations.tsv",
                          sample = c("MR", "MR2", "MR3")) %>%
  Prep.data()
AutoClustR.results.MenonMR <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.MenonMR[[i]]  <- AutoClustR.flow(MenonMR, expr.meas = expr.key$Menon, seed = seed)
  AutoClustR.results.MenonMR[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.MenonMR, file = "Results/AutoClustR/MenonMR.rds")
rm(MenonMR, AutoClustR.results.MenonMR)

# Pollen
Pollen <- Pollen.to.SCexp("Data/Pollen/NBT_hiseq_linear_tpm_values.txt") %>%
  Prep.data()
AutoClustR.results.Pollen <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.Pollen[[i]]  <- AutoClustR.flow(Pollen, expr.meas = expr.key$Pollen, seed = seed)
  AutoClustR.results.Pollen[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.Pollen, file = "Results/AutoClustR/Pollen.rds")
rm(Pollen, AutoClustR.results.Pollen)

# Ranum
Ranum <- Ranum.to.SCexp("Data/Ranum/GSE114157_p15_Expression_Matrix.csv") %>%
  Prep.data()
AutoClustR.results.Ranum <- list()
for (i in 1:10){
  #Seed has to be integer value
  seed <- as.numeric( Sys.time() ) * seeds[i] 
  AutoClustR.results.Ranum[[i]]  <- AutoClustR.flow(Ranum, expr.meas = expr.key$Ranum, seed = seed)
  AutoClustR.results.Ranum[[i]][["seed"]] <- seed
}

saveRDS(AutoClustR.results.Ranum, file = "Results/AutoClustR/Ranum.rds")
rm(Ranum, AutoClustR.results.Ranum)

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
for(i in 1:3) {
  MenonPR.results$CellTrails[ i, ] <- CellTrails.flow(data = MenonPR, expr.meas = expr.key$Menon)
  MenonPR.results$CIDR[ i, ] <- CIDR.flow(data = MenonPR, expr.meas = expr.key$Menon)
  MenonPR.results$IKAP[ i, ] <- IKAP.flow(data = MenonPR, expr.meas = expr.key$Menon,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  MenonPR.results$RaceID[ i, ] <- RaceID.flow(data = MenonPR, expr.meas = expr.key$Menon,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
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

# Menon MR Combined ---- 
MenonMR <- Menon.to.SCexp(data.file.path = "Data/Menon-P/GSE137537_counts.mtx",
                          feat.file.path = "Data/Menon-P/GSE137537_gene_names.txt",
                          cell.file.path = "Data/Menon-P/GSE137537_sample_annotations.tsv",
                          sample = c("MR", "MR2", "MR3")) %>%
  Prep.data()
MenonMR.results <- res.template
seeds <- runif(40)
for(i in 1:3) {
  MenonMR.results$CellTrails[ i, ] <- CellTrails.flow(data = MenonMR, expr.meas = expr.key$Menon)
  MenonMR.results$CIDR[ i, ] <- CIDR.flow(data = MenonMR, expr.meas = expr.key$Menon)
  MenonMR.results$IKAP[ i, ] <- IKAP.flow(data = MenonMR, expr.meas = expr.key$Menon,
                                          seed = (Sys.time() %>% as.numeric()) * seeds[i])
  MenonMR.results$RaceID[ i, ] <- RaceID.flow(data = MenonMR, expr.meas = expr.key$Menon,
                                              seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
  MenonMR.results$SC3[ i, ] <- SC3.flow(data = MenonMR, expr.meas = expr.key$Menon,
                                        seed = (Sys.time() %>% as.numeric()) * seeds[i+20])
  MenonMR.results$Seurat[ i, ] <- Seurat.flow(data = MenonMR, expr.meas = expr.key$Menon,
                                              seed = (Sys.time() %>% as.numeric()) * seeds[i+30])
  saveRDS(MenonMR.results, file = "Results/MenonMR/MenonMR.rds")
}
saveRDS(MenonMR.results, file = "Results/MenonMR/MenonMR.rds")
for(algo in names(MenonMR.results)) {
  write.csv(MenonMR.results[[algo]], file = paste0("Results/MenonMR/MenonMR_", algo, ".csv" ))
}
rm(MenonMR)
# Menon PR Corrected ----


if(file.exists("Data/Menon-P/mnnCorrected.rds")) {
  MenonPR <- readRDS("Data/Menon-P/mnnCorrected.rds")
} else {
  MenonPR <- Menon.to.SCexp.mnn(data.file.path = "Data/Menon-P/GSE137537_counts.mtx",
                                feat.file.path = "Data/Menon-P/GSE137537_gene_names.txt",
                                cell.file.path = "Data/Menon-P/GSE137537_sample_annotations.tsv",
                                sample = c("PR", "PR2", "PR3"))
  saveRDS(MenonPR, "Data/Menon-P/mnnCorrected.rds")
}

MenonPR.results <- res.template
seeds <- runif(40)
for(i in 1:3) {
  MenonPR.results$CellTrails[ i, ] <- CellTrails.flow(data = MenonPR, expr.meas = expr.key$Menon)
  MenonPR.results$CIDR[ i, ] <- CIDR.flow(data = MenonPR, expr.meas = expr.key$Menon)
  MenonPR.results$IKAP[ i, ] <- IKAP.flow(data = MenonPR, expr.meas = expr.key$Menon,
                                         seed = (Sys.time() %>% as.numeric()) * seeds[i])
  MenonPR.results$RaceID[ i, ] <- RaceID.flow(data = MenonPR, expr.meas = expr.key$Menon,
                                             seed = (Sys.time() %>% as.numeric()) * seeds[i+10])
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