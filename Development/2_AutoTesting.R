#Automated Testing
library(stats)
library(utils)
library(ggplot2)
library(clusterCrit, aricode)
library(aricode)

#Load in datasets, assemble into list
datasets <- list(zeisel,pr, mr)
for(object in datasets){
  Predict_nPCs(object) %>%
    print()
}


#Create list of arrays to store results
index.frame <- data.frame(matrix(nrow = 24^2, ncol = 12))
names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std",
                        "CH", "DB", "Sil", 
                        "ARI", "AARI", "Cluster_Number",
                        "Cluster_Error", "K.Param", "Resolution")
results <- list()
for (i in 1:3) {results[[i]] <- index.frame}
names(results) <- c("Zeisel", "PR", "MR")
times <- numeric(length = 3)
names(times) <- c("Zeisel", "PR", "MR")

#Set Parameters
k.param <- 2^(9:32/4) %>%
  round() %>%
  rep(each = 24)
resolution <- rep(2^(9:32/4)/100, 24)

#Iterate through datasets 
j <- 1
for (object in datasets) {
  npc <- Predict_nPCs(object)
  embeddings <- PopulateEmbeddings(object, npc)
  clusters <- as.integer(object@active.ident)
  cluster.number <- unique(clusters) %>%
    length()
  
  print("Performing Seurat Simulations")
  
  start.time <- Sys.time()
  pb <- txtProgressBar(min = 1, max = 24, initial = 1, char = "+", style = 3)
  for (i in 1:576){
    #Seurat updated, default nn.method is now "annoy"
    temp.object <- FindNeighbors(object,
                                 k.param = k.param[i],
                                 compute.SNN = TRUE,
                                 prune.SNN = 1/15,
                                 nn.method = "annoy",
                                 nn.eps = 0.0,
                                 verbose = FALSE,
                                 reduction = "pca",
                                 dims = 1:npc) %>%
      
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
      index.frame[i, 1:3] <- intCriteria(embeddings$STD,
                                         temp.clusters,
                                         c("Calinski", "Davies", "Sil")) %>%
        unlist()
      index.frame[i, 4:6] <- intCriteria(embeddings$NonSTD,
                                         temp.clusters,
                                         c("Calinski", "Davies", "Sil")) %>%
        unlist()
    }
    
    #Write Results
    index.frame[i, 7] <- ARI(temp.clusters, clusters)
    index.frame[i, 9:10] <- c(max(temp.clusters), abs(max(temp.clusters)-cluster.number))
    index.frame[i, 8] <- index.frame[i, 7]*(1-index.frame[i,10]/cluster.number) #AARI
    index.frame[i, 11] <- k.param[i]
    index.frame[i, 12] <- resolution[i]
    
    setTxtProgressBar(pb, i/24)
  }
  close(pb)
  
  results[[j]] <- index.frame
  end.time <- Sys.time() - start.time
  print(end.time)
  times[[j]] <- end.time
  j <- j + 1
}


#Summarize results
#Summarize spearmans correlation
makeRhoFrame <- function(index.frame) {
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
  return(rho.frame)
}
rhos <- lapply(results, makeRhoFrame)
rho.ari <- lapply(rhos, pull, var = "ARI") %>%
  data.frame()
avg.rho <- rowSums(rho.ari)/ncol(rho.ari)
rho.ari <- cbind(rho.ari, avg.rho)
write.table(rho.ari, file = 'clipboard', sep = '\t')


#Summarize ARI @ maximum
getMax <- function(index.frame){
  index.frame <- index.frame[c(1:7, 10, 8, 9, 11, 12)]
  max.index <- data.frame(matrix(nrow = 8, ncol = 12))
  names(max.index) <- names(index.frame)
  row.names(max.index) <- names(index.frame)[1:8]
  for (i in 1:8) {
    max.index[i, ] <- index.frame[which.max(index.frame[[i]]),]
  }
  max.index[2, ] <- index.frame[which.min(index.frame[[2]]), ]
  max.index[5, ] <- index.frame[which.min(index.frame[[5]]), ]
  return(max.index)
}
maxed <- lapply(results, getMax)
max.ARI <- lapply(maxed, pull, var = "ARI") %>%
  data.frame()
avg.ARI <- rowSums(max.ARI)/ncol(max.ARI)
max.ARI <- cbind(max.ARI, avg.ARI)
write.table(max.ARI, file = 'clipboard', sep = '\t')

#Summarize cluster number error @ maximum
getClusterError <- function(index.frame) {
  true.clusterNumber <- max(index.frame$Cluster_Number) - index.frame[which.max(index.frame$Cluster_Number), "Cluster_Error"]
  max.error <- data.frame(Cluster_Error = matrix(nrow = 7, ncol = 1))
  for (i in 1:7){
    max.error[i,1] <- index.frame[which.max(index.frame[[i]]), "Cluster_Number"] - true.clusterNumber
  }
  return(max.error)
}
cluster.error.table <- lapply(results, getClusterError) %>%
  data.frame(row.names = names(results[[1]])[1:7])
names(cluster.error.table) <- c("Zeisel", "PR", "MR")
write.table(cluster.error.table, file = 'clipboard', sep = '\t')


#Store Results as .rds
npc_by_obj <- lapply(datasets, Predict_nPCs)
names(npc_by_obj) <- c("Zeisel", "PR", "MR")
results.list <- list(results, rhos, maxed, npc_by_obj)
names(results.list) <- c("index.frame", "Rho.Table", "MaxValues", "NPCs")
saveRDS(results.list, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/New_Predict_exp_params.rds")


#Save the the index frames together, named, and in order
index.list <- vector(mode = "list", length = 5)
index.list[c(1,4,5)] <- results
results.predict <-  readRDS("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Predict_pcs_exp_params.rds")
index.list[2:3] <- results.predict$index.frame[2:3]
names(index.list) <- c("Zeisel", "Panc1", "Panc2", "PR", "MR")
saveRDS(index.list, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/New_Predict_Index_List.rds")

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
      if(abs(model$coefficients[[2]]) < 0.84){
        PE.table[j-1, i-3] <- (SD.table$SD[j-1] - y.hat) /
          sqrt(sum(residuals(model)^2) / (nrow(temp.table) - 2))
      } else{ 
        PE.table[j-1, i-3] <- 0 
      }
    }
  }
  best.npc <- apply(PE.table, 1, function(x) {mean(x[x!=0])})  %>%
    which.max()
  if(best.npc == 1){best.npc <- 2}
  #median.period <- round(pc.use - best.npc/2)
  return(best.npc)
}


#PopulateEmbeddings
PopulateEmbeddings <- function(object, npc){
  all.embeddings <- Embeddings(object)[, 1:npc]
  std_dev <- object[["pca"]]@stdev[1:npc]
  std.embeddings <- t(t(apply(all.embeddings, MARGIN = 2,
                              FUN = function(x) {x - mean(x)})) # subtracts PC mean from PC position
                      /std_dev) # Divide by standard deviation
  embeds <- list(all.embeddings, std.embeddings)
  names(embeds) <- c("NonSTD", "STD")
  return(embeds)
}











