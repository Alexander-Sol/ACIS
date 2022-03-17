# Results Plotting
library(data.table)
library(dplyr)

pal.list <- list(
  Seurat = "#B1D46D",
  SC3 = "#F8B362",
  RaceID = "80B1D3",
  IKAP = "#F37F71",
  CIDR = "#BDB9D7",
  CellTrails = "#F9F5B6",
  AutoClustR = "#8DD0C5"
)

ARI <- sapply(Baron1, FUN = function(x) pluck(x, "post.sc", "ARI"))

ac2csv <- function(results){
  ARI <- sapply(results, FUN = function(x) pluck(x, "post.sc", "ARI"))
  n.clusters <- sapply(results, FUN = function(x) pluck(x, "post.sc", "n.clusters"))
  runtime <- sapply(results, FUN = function(x) pluck(x, "totalRuntime"))
  seed <- sapply(results, FUN = function(x) pluck(x, "seed"))
  return(
    data.frame(ARI, n.clusters, runtime, seed)
  )
}

dataNames <- c("Baron2", "Goolam", "Kolodz", "Loh", "Pollen", "Ranum")
for (name in dataNames) {
  temp.csv <- readRDS(paste0("Results/AutoClustR/", name, ".rds")) %>% ac2csv()
  write.csv(temp.csv, file = paste0("Results/", name, "/", name, "_AutoClustR.csv"))
}

test <- readRDS("Results/AutoClustR/MenonPR.rds") %>% ac2csv()
write.csv(test, "Results/MenonPR/MenonPR_SmallAutoClustR.csv")


barRes <- data.frame(matrix(data = 0, nrow = 10, ncol = 7))
names(barRes) <- names(pal.list)

barRes <- data.frame(matrix(data = 0, nrow = 70, ncol = 4))
names(barRes) <- c("Platform", "ARI", "Cluster_Number", "Runtime")
test <- cbind(
  rep(algo, 10),
  fread(paste0("Results/Baron1/Baron1_", algo, ".csv"))[ , -1]
)

for (algo in names(pal.list)){
  test <- rbind(test, cbind(
    rep(algo, 10),
    fread(paste0("Results/Baron1/Baron1_", algo, ".csv"))[ , -1]
  ))
}

avgRes <- data.frame()


