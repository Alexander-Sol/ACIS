#tSNE kMeans testing for training datasets
#Automated Testing
library(stats)
library(utils)
library(ggplot2)
library(clusterCrit, aricode)
library(aricode)
library(Rtsne)

pullPCA <- function(object, dims) {
  PCA.mat <- Embeddings(object)[,dims]
  return(PCA.mat)
}

Predict_nPCs <- function(object, pc.use = 20) {
  SD.table <- data.frame(PC = 1:pc.use,
                         SD = object[["pca"]]@stdev[1:pc.use])
  PE.table <- data.frame(matrix(nrow = pc.use-4,
                                ncol = pc.use-3,
                                data = 0))
  for(i in 5:pc.use){
    for(j in (i - 3):2){
      temp.table <- SD.table[j:i,]
      model <- lm(SD ~ PC, temp.table)
      y.hat <- model$coefficients[[2]] * (j - 1) + model$coefficients[[1]]
      # This used to be wrapped in an if/else clause that would set PE to 
      # 0 if the slope was greater than 0.84 
      PE.table[j-1, i-3] <- (SD.table$SD[j-1] - y.hat) /
        sqrt(sum(residuals(model)^2) /
               (nrow(temp.table) - 2))
    }
  }
  best.npc <- apply(PE.table, 1, function(x) {mean(x[x!=0])})  %>%
    which.max()
  if(best.npc == 1){best.npc <- 2}
  return(best.npc)
}

#Load in datasets, assemble into list
datasets <- list(Zeisel = zeisel, Panc1 = panc.1, Panc2 = panc.2, PR = pr, MR = mr)
for(object in datasets){
  Predict_nPCs(object) %>%
    print()
}

#Create list of arrays to store results
index.frame <- data.frame(matrix(nrow = 24^2, ncol = 7))
names(index.frame) <- c("CH", "DB", "Sil", "ARI",
                        "K", "Perplexity", "ClusterError")
results <- vector(mode = "list", length = 5)
for (i in 1:5) {results[[i]] <- index.frame}
names(results) <- c("Zeisel", "Panc1", "Panc2", "PR", "MR")
# times <- numeric(length = 3)
# names(times) <- c("Zeisel", "PR", "MR")

# Set parameters
k_clusters <- 2:25
perplexity.space <- 1:24*5

j <- 1
for (dataset in datasets) {
  npc <- Predict_nPCs(dataset)
  embedPCA <- pullPCA(dataset, 1:npc)
  pb <- txtProgressBar(min = 1, max = 24, initial = 1, char = "+", style = 3)
  for(perp in perplexity.space) {
    embedTSNE <- Rtsne(X = embedPCA, pca = F, check_duplicates = F, dims = 3,
                       perplexity = perp)$Y
    cluster.list <- lapply(k_clusters,
                           FUN = function(k) {
                             kmeans(embedTSNE, centers = k, nstart = 25)$cluster
                             })
    icvi.list <- lapply(cluster.list,
                        FUN = function(cluster.idents) {
                          intCriteria(embedPCA, part = cluster.idents, crit = c("Cal", "Dav", "Sil")) %>% unlist()
                        })
    index.frame.row <- 1:24 + 24*(perp/5 - 1)
    index.frame[index.frame.row, 1:3] <- do.call(rbind.data.frame, icvi.list)
    index.frame[index.frame.row, 4] <- sapply(cluster.list, ARI, c2 = dataset$annotation)
    index.frame[index.frame.row, 5] <- k_clusters
    index.frame[index.frame.row, 6] <- rep(perp, 24)
    index.frame[index.frame.row, 7] <- k_clusters - length(unique(dataset$annotation))
    setTxtProgressBar(pb, perp/5)
  }
  close(pb)
  results[[j]] <- index.frame
  j <- j + 1
}

# Save Results
saveRDS(results, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/tSNE_kMeans_Results.rds")

# Find the maximum ARI returned by each function
maxARI <- function(index.frame) {
  ch.max <- index.frame$ARI[which.max(index.frame$CH)]
  db.max <- index.frame$ARI[which.min(index.frame$DB)]
  si.max <- index.frame$ARI[which.max(index.frame$Sil)]
  return(c(ch.max, db.max, si.max))
}
maxARI.res <- lapply(results, maxARI)
maxARI.table <- do.call(rbind.data.frame, maxARI.res)
names(maxARI.table) <- c("Calinski-Harabsz", "Davies-Bouldin", "Silhouette")
maxARI.table <- rbind(maxARI.table, colSums(maxARI.table)/5)
row.names(maxARI.table) <- c(names(results), "Average_ARI@Max")
write.table(maxARI.table, file = 'clipboard', sep = '\t')

# Find the relationship between ICVI and ARI for each dataset
makeRhoFrame <- function(index.frame) {
  ch.cor <- cor.test(index.frame$ARI, index.frame$CH, method = "spearman")$estimate
  db.cor <- cor.test(index.frame$ARI, index.frame$DB, method = "spearman")$estimate
  si.cor <- cor.test(index.frame$ARI, index.frame$Sil, method = "spearman")$estimate
  return(c(ch.cor, db.cor, si.cor))
}
rhos <- lapply(results, makeRhoFrame)
rho.table <- do.call(rbind.data.frame, rhos)
names(rho.table) <- c("Calinski-Harabsz", "Davies-Bouldin", "Silhouette")
avg.rho <- colSums(rho.table)/5
rho.table <- rbind(rho.table, avg.rho)
row.names(rho.table) <- c(names(results), "Average_Rho")
write.table(rho.table, file = 'clipboard', sep = '\t')


#%%%%%%%%%%%%%%%%%%%%%%%   Visualization   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Make heatmap w/ local maxima for CH Index   %%%%%%%%%%%%%
ggplot(results$MR, aes_string(x = "K", y = "Perplexity", fill = "CH")) +
  geom_tile() +
  scale_fill_gradient(low = "gold1", high = "darkred") +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18)) +
  labs(x = "K",
       y = "Perplexity", fill = "CH") + 
  ggtitle(label = "MR Dataset: CH Index")

table(grid.res.ch$best)


# Make heatmap w/ local maxima for Silhouette Index   %%%%%%%%%%%%%
ggplot(grid.res.sil, aes_string(x = "logK", y = "logRes", fill = "Sil")) +
  geom_tile() +
  scale_fill_gradient(low = "gold1", high = "darkred") +
  geom_text(aes_string(label = "best"), size = 6) +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18)) +
  labs(x = expression("\nLog"[2]*"K.param"),
       y = expression("Log"[2]*"(Resolution*100)"), fill = "Sil") + 
  ggtitle(label = "MR Dataset: Silhouette Index")

table(grid.res.sil$best)

k <- 1
for (index.frame in results) {
  index.frame$bestDB <- ""
  index.frame$bestDB[which.min(index.frame$DB)] <- "X"
  db.plot <- ggplot(index.frame, aes_string(x = "K", y = "Perplexity", fill = "DB")) +
    geom_tile() +
    scale_fill_gradient(low = "darkred", high = "gold1") +
    geom_text(aes_string(label = "bestDB"), size = 6) +
    theme(axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18)) +
    labs(x = "K", y = "Perplexity", fill = "DB") + 
    ggtitle(label = paste0(names(results)[k], " Dataset: DB Index"))
  print(db.plot)
  
  index.frame$bestSil <- ""
  index.frame$bestSil[which.max(index.frame$Sil)] <- "X"
  sil.plot <- ggplot(index.frame, aes_string(x = "K", y = "Perplexity", fill = "Sil")) +
    geom_tile() +
    scale_fill_gradient(low = "gold1", high = "darkred") +
    geom_text(aes_string(label = "bestSil"), size = 6) +
    theme(axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18)) +
    labs(x = "K", y = "Perplexity", fill = "Sil") + 
    ggtitle(label = paste0(names(results)[k], " Dataset: Silhouette Index")) %>%
    print()
  print(sil.plot)
    
  k <- k + 1
}







