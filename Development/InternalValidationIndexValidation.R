#Testing correlation between internal and external cluster validation indices

library(dplyr)
library(Seurat)
library(AutoClustR)
library(clusterCrit)
library(aricode)

panc.1 <- read.csv("C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet1/GSM2230757_human1_umifm_counts.csv/GSM2230757_human1_umifm_counts.csv")
panc.1 <- t(panc.1)
colnames(panc.1) <- panc.1[2, ]
panc.1.annotation <- panc.1[3,]
table(as.character(panc.1.annotation))
panc.1 <- panc.1[-c(1,2,3), ]
panc.1.seurat <- CreateSeuratObject(panc.1, project = "HumanPancreas1")
panc.1.seurat$annotation <- panc.1.annotation
table(panc.1.seurat$annotation)
p1s <- panc.1.seurat
p1s$na <- !is.na(p1s$annotation)
p1s <- subset(p1s, na)
rm(panc.1)

length(Cells(p1s))
p1s <- SCTransform(panc.1.seurat)
p1s <- RunPCA(p1s)
p1s <- Choose_nPCs(p1s)
#Best npcs = 5
p1s <- LoadEmbeddings(p1s, standardize = T)
p1s <- RunUMAP(p1s, n.neighbors = 20, dims = 1:5)

p1s@active.ident <- as.factor(p1s$annotation)
DimPlot(p1s)
length(Cells(p1s))
clusters <- as.integer(Idents(p1s))
class(clusters)
class(cl$cluster)
embeddings <- p1s@tools$LoadEmbeddings



clusters <- as.integer(Idents(p1s))
#Only have to load embeddings once
embeddings <- p1s@tools$LoadEmbeddings
#If this is breaking, make sure that na's have been removed before loading embeddings, defining clusters
ind <- intCriteria(embeddings,
                   part = clusters,
                   crit = "Calinski_Harabasz")

ari <- ARI(clusters, p1s$annotation)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Reverse Louvain clustering
ch.index <- numeric(length = 1905)
sil.index <- numeric(length = 1905)
db.index <- numeric(length = 1905)
ari.index <- numeric(length = 1905)
index.frame <- data.frame(ch.index, sil.index, db.index, ari.index)
names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")
new.clusters <- clusters
pb <- txtProgressBar(min = 1, max = 95, char = "+")
for (i in 1:length(clusters)){
  new.clusters[i] <- sample(1:14, 1)
  index.frame[i, 1:3] <- intCriteria(embeddings,
                     part = new.clusters,
                     crit = c("calinski", "davies", "sil")) %>%
    unlist()
  index.frame[i, 4] <- ARI(new.clusters, p1s$annotation)
  setTxtProgressBar(pb, value = i/20)
}

ggplot(index.frame, aes(ARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("Reverse Louvain Clustering, CH_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$CH_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("Reverse Louvain Clustering, DB_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$DB_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("Reverse Louvain Clustering, Sil_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$Sil_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Random Mutations
for (i in 1:150){
  new.clusters <- clusters
  reassign <- sample(1:1905, 1)
  idx <- sample(1:1905, reassign)
  clust.num <- sample(9:19, 1)
  new.clusters[idx] <- sample(1:clust.num, reassign, replace = T) 
  ch.index[i] <- intCriteria(embeddings,
                             part = as.integer(new.clusters),
                             crit = "calinski")
  db.index[i] <- intCriteria(embeddings,
                            part = as.integer(new.clusters),
                            crit = "davies")
  sil.index[i] <- intCriteria(embeddings,
                             part = as.integer(new.clusters),
                             crit = "sil")
  ari[i] <- ARI(new.clusters, p1s$annotation)
}

ggplot(index.frame, aes(ARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("Random Clustering, CH_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$CH_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("Random Clustering, DB_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$DB_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("Random Clustering, Sil_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$Sil_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#k Means simulation study
ch.index <- numeric(length = 100)
sil.index <- numeric(length = 100)
db.index <- numeric(length = 100)
ari.index <- numeric(length = 100)
index.frame <- data.frame(ch.index, sil.index, db.index, ari.index)
names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")

pb <- txtProgressBar(min = 1, max = 75, initial = 1, char = "+", style = 2)
for (i in 1:150){
  cl <- kmeans(embeddings, sample(9:19, 1))
  index.frame[i, 1:3] <- intCriteria(embeddings,
                      cl$cluster,
                      c("Calinski", "Davies", "Sil")) %>%
    unlist()
  index.frame[i, 4] <- ARI(cl$cluster, p1s$annotation)
  setTxtProgressBar(pb, i/2)
}

ggplot(index.frame, aes(ARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("K-Means Clustering, CH_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$CH_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("K-Means Clustering, DB_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$DB_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("K-Means Clustering, Sil_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$Sil_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#graph based community detection simulation
#Example, before looping
library(RANN)
library(igraph)

ch.index <- numeric(length = 100)
sil.index <- numeric(length = 100)
db.index <- numeric(length = 100)
ari.index <- numeric(length = 100)
index.frame <- data.frame(ch.index, sil.index, db.index, ari.index)
names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")

pb <- txtProgressBar(min = 1, max = 75, initial = 1, char = "+", style = 2)
k.list <- sample(10:240, 150)
for (i in 1:150){
  #graph based clustering
  knn.info <- nn2(embeddings, k = k.list[i])
  knn <- knn.info$nn.idx
  adj <- matrix(0, nrow(embeddings), nrow(embeddings))
  rownames(adj) <- colnames(adj) <- rownames(embeddings)
  for(j in seq_len(nrow(embeddings))){
    adj[j, rownames(embeddings)[knn[j, ]]] <- 1
  }
  g <- graph.adjacency(adj, mode = "undirected")
  g <- simplify(g)
  km <- cluster_walktrap(g)
  com <- km$membership
  
  index.frame[i, 1:3] <- intCriteria(embeddings,
                                     as.integer(com),
                                     c("Calinski", "Davies", "Sil")) %>%
    unlist()
  index.frame[i, 4] <- ARI(com, p1s$annotation)
  setTxtProgressBar(pb, i/2)
}


ggplot(index.frame, aes(ARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("KNN Clustering, CH_STD vs ARI: Rho = ", 
                cor.test(x = index.frame$ARI,
                         y = index.frame$CH_Std,
                         method = "spearman")$estimate %>%
                  round(digits=4)))


ggplot(index.frame, aes(ARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("KNN Clustering, DB_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$DB_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("KNN Clustering, Sil_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$Sil_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#SNN simulation study
install.packages("dbscan")
library(dbscan)
ch.index <- numeric(length = 100)
sil.index <- numeric(length = 100)
db.index <- numeric(length = 100)
ari.index <- numeric(length = 100)
index.frame <- data.frame(ch.index, sil.index, db.index, ari.index)
names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")

eps <- sample(2:10, 1)
snn.test <- sNNclust(embeddings, k = sample(10:50, 1), eps = eps, minPts = eps*2)

pb <- txtProgressBar(min = 1, max = 75, initial = 1, char = "+", style = 2)
for (i in 1:150){
  eps <- sample(2:10, 1)
  snn.test <- sNNclust(embeddings, k = sample(10:50, 1), eps = eps, minPts = eps*2)
  index.frame[i, 1:3] <- intCriteria(embeddings,
                                     as.integer(snn.test$cluster),
                                     c("Calinski", "Davies", "Sil")) %>%
    unlist()
  index.frame[i, 4] <- ARI(snn.test$cluster, p1s$annotation)
  setTxtProgressBar(pb, i/2)
}

ggplot(index.frame, aes(ARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("SNN Clustering, CH_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$CH_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("SNN Clustering, DB_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$DB_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))


ggplot(index.frame, aes(ARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("SNN Clustering, Sil_STD vs ARI: Rho = ", 
                 cor.test(x = index.frame$ARI,
                          y = index.frame$Sil_Std,
                          method = "spearman")$estimate %>%
                   round(digits=4)))



test.sim.func <- ICVI_Simulation(embeddings, clusters, "Pancreatic Islets 1", cl.range = 5)


#Do the same thing for non-standardized embeddings
p1s.ns <- LoadEmbeddings(p1s, standardize = F)
embeddings <- p1s.ns@tools$LoadEmbeddings

test.sim.func.mbs <- ICVI_Simulation(embeddings, clusters, name = "Pancreatic Islets 1 Non Standardized", cl.range = 5, iterations = 150)
test.sim.func[[1]]
