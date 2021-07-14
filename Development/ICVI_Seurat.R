#Check ICVI vs Seurat for arbitrary datasets
library(Seurat)
library(dplyr)
library(stats)
library(utils)
library(ggplot2)
library(clusterCrit, aricode)
library(aricode)



object <- p1s
folder <- "C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet2/"


object <- Choose_nPCs(object)
#Load in embeddings used to calculate ICVI's
object <- LoadEmbeddings(object, standardize = T)
std.embed <- object@tools$LoadEmbeddings
object <- LoadEmbeddings(object, standardize = F)
nonstd.embed <- object@tools$LoadEmbeddings

#Load in annotated cluster identity
clusters <- as.factor(object$annotation) %>%
  as.integer()
cluster.number <- length(unique(clusters))

#AARI can be defined as ARI*(1-cluster.error/true.cluster.number)
iterations <- 24*24
index.frame <- data.frame(matrix(nrow = length(iterations), ncol = 12))
names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std",
                        "CH", "DB", "Sil", 
                        "ARI", "Cluster_Number", "Cluster_Error",
                        "AARI", "K.Param", "Resolution")

'k.param <- rep(c(3,5,7, 10, 20, 30), each = 24)
k.param <- rev(k.param)
resolution <- rep((1:24)/20, 6)'

#Trying log2 distributed parameters

k.param <- 2^(9:32/4) %>%
  round() %>%
  rep(each = 24)
resolution <- rep(2^(9:32/4)/100, 24)





print("Performing Seurat Simulations")


start.time <- Sys.time()
pb <- txtProgressBar(min = 1, max = iterations/10, initial = 1, char = "+", style = 3)
for (i in 1:iterations){
  #Seurat updated, default nn.method is now "annoy"
  temp.object <- FindNeighbors(object,
                               k.param = k.param[i],
                               compute.SNN = TRUE,
                               prune.SNN = 1/15,
                               nn.method = "annoy",
                               nn.eps = 0.0,
                               verbose = FALSE,
                               reduction = "pca",
                               dims = 1:object@tools$Choose_nPCs) %>%
    
    # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
    FindClusters(modularity.fxn = 1,
                 resolution = resolution[i],
                 algorithm = 1,
                 group.singletons = TRUE,
                 verbose = FALSE)
  
  # Determine and store cluster number for ("K.param", "Resolution") = (i, j)
  temp.clusters <- as.integer(Idents(temp.object))
  K <- length(levels(temp.object))
  
  # Calculate and store internal clustering validation index for ("K.param", "Resolution") = (i, j)
  if(K == 1){
    index.frame[i, 1:6] <- NA
  } else{
    index.frame[i, 1:3] <- intCriteria(std.embed,
                                       temp.clusters,
                                       c("Calinski", "Davies", "Sil")) %>%
      unlist()
    index.frame[i, 4:6] <- intCriteria(nonstd.embed,
                                       temp.clusters,
                                       c("Calinski", "Davies", "Sil")) %>%
      unlist()
  }
  
  #Write Results
  index.frame[i, 7] <- ARI(temp.clusters, clusters)
  index.frame[i, 8:9] <- c(max(temp.clusters), abs(max(temp.clusters)-cluster.number))
  index.frame[i, 10] <- index.frame[i, 7]*(1-index.frame[i,9]/cluster.number) #AARI
  index.frame[i, 11] <- k.param[i]
  index.frame[i, 12] <- resolution[i]
  
  setTxtProgressBar(pb, i/10)
}

#index.frame$AARI <- index.frame$ARI*(1-index.frame$Cluster_Error/cluster.number) #AARI
close(pb)
end.time <- Sys.time() - start.time
print(end.time)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot.list <- list()

#AARI Std

rho.table <- c(cor.test(x = index.frame$AARI,
                        y = index.frame$CH_Std,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$AARI,
                        y = index.frame$DB_Std,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$AARI,
                        y = index.frame$Sil_Std,
                        method = "spearman")$estimate)

plot.list[[1]] <- ggplot(index.frame, aes(AARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, CH_STD vs AARI: Rho = ", 
                 rho.table[1] %>%
                   round(digits=4)))


plot.list[[2]] <- ggplot(index.frame, aes(AARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, DB_STD vs AARI: Rho = ", 
                 rho.table[2] %>%
                   round(digits=4)))


plot.list[[3]] <- ggplot(index.frame, aes(AARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, Sil_STD vs AARI: Rho = ", 
                 rho.table[3] %>%
                   round(digits=4)))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARI Std %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho.table <- c(cor.test(x = index.frame$ARI,
                        y = index.frame$CH_Std,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$DB_Std,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$Sil_Std,
                        method = "spearman")$estimate)


plot.list[[4]] <- ggplot(index.frame, aes(ARI, CH_Std)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, CH_STD vs ARI: Rho = ", 
                 rho.table[1] %>%
                   round(digits=4)))



plot.list[[5]] <- ggplot(index.frame, aes(ARI, DB_Std)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, DB_STD vs ARI: Rho = ", 
                 rho.table[2] %>%
                   round(digits=4)))


plot.list[[6]] <- ggplot(index.frame, aes(ARI, Sil_Std)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, Sil_STD vs ARI: Rho = ", 
                 rho.table[3] %>%
                   round(digits=4)))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AARI Non-Std %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho.table <- c(cor.test(x = index.frame$AARI,
                        y = index.frame$CH,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$AARI,
                        y = index.frame$DB,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$AARI,
                        y = index.frame$Sil,
                        method = "spearman")$estimate)

plot.list[[7]] <- ggplot(index.frame, aes(AARI, CH)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, CH vs AARI: Rho = ", 
                 rho.table[1] %>%
                   round(digits=4)))



plot.list[[8]] <- ggplot(index.frame, aes(AARI, DB)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, DB vs AARI: Rho = ", 
                 rho.table[2] %>%
                   round(digits=4)))


plot.list[[9]] <- ggplot(index.frame, aes(AARI, Sil)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, Sil vs AARI: Rho = ", 
                 rho.table[3] %>%
                   round(digits=4)))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARI Non-Std %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho.table <- c(cor.test(x = index.frame$ARI,
                        y = index.frame$CH,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$DB,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$Sil,
                        method = "spearman")$estimate)

plot.list[[10]] <- ggplot(index.frame, aes(ARI, CH)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, CH vs ARI: Rho = ", 
                 rho.table[1] %>%
                   round(digits=4)))



plot.list[[11]] <- ggplot(index.frame, aes(ARI, DB)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, DB vs ARI: Rho = ", 
                 rho.table[2] %>%
                   round(digits=4)))


plot.list[[12]] <- ggplot(index.frame, aes(ARI, Sil)) + 
  geom_point() + 
  ggtitle(paste0("Menon PR1 Seurat Clustering, Sil vs ARI: Rho = ", 
                 rho.table[3] %>%
                   round(digits=4)))



plot.list[[13]] <- ggplot(index.frame, aes(ARI, Cluster_Number)) + 
  geom_point() + 
  ggtitle("Menon PR1 Seurat Clustering, Cluster Number vs ARI")




for (i in 1:13) {
  save.plot(plot.list[[i]], folder = folder, name = paste0("Menon_PR1_SeuratSim_", i))
}

save.plot <- function(plot, folder, name){
  ggsave(filename = paste0(folder,  name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 5,
         width = 8,
         dpi = 300)
  dev.off()
}

rho.frame <- data.frame(matrix(nrow = 6, ncol = 2))
names(rho.frame) <- c("ARI", "AARI")
rownames(rho.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "CH", "DB", "Sil")
rho.frame[ , 1] <- c(cor.test(x = index.frame$ARI,
                        y = index.frame$CH_Std,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$DB_Std,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$Sil_Std,
                        method = "spearman")$estimate,
                cor.test(x = index.frame$ARI,
                        y = index.frame$CH,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$DB,
                        method = "spearman")$estimate,
               cor.test(x = index.frame$ARI,
                        y = index.frame$Sil,
                        method = "spearman")$estimate)
rho.frame[ , 2] <- c(cor.test(x = index.frame$AARI,
                              y = index.frame$CH_Std,
                              method = "spearman")$estimate,
                     cor.test(x = index.frame$AARI,
                              y = index.frame$DB_Std,
                              method = "spearman")$estimate,
                     cor.test(x = index.frame$AARI,
                              y = index.frame$Sil_Std,
                              method = "spearman")$estimate,
                     cor.test(x = index.frame$AARI,
                              y = index.frame$CH,
                              method = "spearman")$estimate,
                     cor.test(x = index.frame$AARI,
                              y = index.frame$DB,
                              method = "spearman")$estimate,
                     cor.test(x = index.frame$AARI,
                              y = index.frame$Sil,
                              method = "spearman")$estimate)

write.table(rho.frame, "clipboard", sep="\t")

index.frame <- index.frame[c(1:7, 10, 8, 9, 11, 12)]


write.csv(index.frame, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet2/Indices.csv")

max.index <- data.frame(matrix(nrow = 8, ncol = 12))
names(max.index) <- names(index.frame)
row.names(max.index) <- names(index.frame)[1:8]
for (i in 1:8) {
  max.index[i, ] <- index.frame[which.max(index.frame[[i]]),]
}
max.index[2, ] <- index.frame[which.min(index.frame[[2]]), ]
max.index[5, ] <- index.frame[which.min(index.frame[[5]]), ]

View(max.index)
write.csv(max.index, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/PancreaticIslet2/MaxIndices.csv")
write.table(max.index, "clipboard", sep = "\t")



index.frame <- read.csv("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Menon2019_GSE137537/Macular1/Indices.csv")
index.frame <- index.frame[-1]









