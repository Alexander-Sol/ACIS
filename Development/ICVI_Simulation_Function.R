#ICVI Generic
library(dbscan)
library(RANN)
library(igraph)
library(dplyr)
library(stats)
library(utils)
library(Seurat)
library(ggplot2)
library(clustercrit, aricode)

#Example Workflow
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
rm(panc.1)

#Didn't do any QC, probably should for final versions
p1s <- SCTransform(panc.1.seurat)
p1s <- RunPCA(p1s)
p1s <- Choose_nPCs(p1s)
#Best npcs = 5
p1s$na <- !is.na(p1s$annotation)
p1s <- subset(p1s, na)
p1s <- LoadEmbeddings(p1s, standardize = T)
embeddings <- p1s@tools$LoadEmbeddings
clusters <- as.factor(p1s$annotation) %>%
  as.integer()

#Clustering is performed on the standardized PCs, which should be changed. 
#Non Standardized PCs show the same trends, just with lower rho values across the board
sim.res <- ICVI_Simulation(embeddings = embeddings,
                           clusters = clusters,
                           name = "Pancreatic Islets 1",
                           cl.range = 5,
                           iterations = 100)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Functions taken from AutoClustR

Choose_nPCs <- function(object,
                        pc.use = 50,
                        pc.select = "predict",
                        pc.estimate = 25){
  #Run PCA, if required
  object <- CheckPCs(object, pc.use)
  
  if (pc.select == "predict") {
    results <- Predict_nPCs(object, pc.use = pc.use)
    object <- results[[1]]
    #Store best nPCs in correct slot
    Tool(object) <- results[[2]]
  } else if (pc.select == "segment") {
    results <- Segment_nPCs(object, pc.use = pc.use, pc.estimate = pc.estimate)
    object <- results[[1]]
    Tool(object) <- results[[2]]
  } else if (pc.select == "manual") {
    Tool(object) <- pc.use
  }
  return(object)
}

CheckPCs <- function(object, pc.use){
  if(is.null(object@reductions$pca)){
    print("Calculating principal components")
    object <- RunPCA(object, npcs = pc.use, verbose = FALSE)
  } else if(length(object[["pca"]]@stdev) < pc.use){
    print("Calculating additional principal components")
    object <- RunPCA(object, npcs = pc.use, verbose = FALSE)
  }
  return(object)
}

Predict_nPCs <- function(object, pc.use = 50) {
  
  # Load principal components and singular values
  pcs <- 1:pc.use
  std.dev <- object[["pca"]]@stdev[pcs]
  
  # Create empty table for storing SE-mod results
  npc <- data.frame()
  
  # Create progress bar
  print("Determining dimensionality.")
  prog.bar <- txtProgressBar(min = 0, max = (pc.use - 4), char = "+", style = 3)
  
  # Perform SE-mod
  for(i in 5:pc.use){
    se.mod <- data.frame()
    count <- 1
    
    for(j in (i - 3):2){
      temp.table <- as.data.frame(cbind(j:i, std.dev[j:i]))
      n <- length(temp.table[, 1])
      colnames(temp.table) <- c("PC", "SD")
      model <- lm(SD ~ PC, temp.table)
      y.hat <- model$coefficients[[2]] * (j - 1) + model$coefficients[[1]]
      se.mod[count, 1] <- j - 1
      if(abs(model$coefficients[[2]]) < 0.84){
        se.mod[count, 2] <- (std.dev[[j - 1]] - y.hat) / sqrt(sum(residuals(model)^2) / (n - 2))
      } else{se.mod[count, 2] <- 0}
      count <- count + 1
    }
    
    npc[(i - 4), 1] <- i
    if(max(se.mod[, 2], na.rm = TRUE) == 0){npc[(i - 4), 2] <- NA}else{npc[(i - 4), 2] <- se.mod[se.mod[, 2] == max(se.mod[, 2],                                                                                                                na.rm = TRUE), 1]}
    
    # Update progress bar
    setTxtProgressBar(prog.bar, (i - 4))
  }
  
  #Add newline after progress bar
  writeLines("\n")
  
  # Determine dimensionality
  npc <- npc[is.na(npc[, 2]) == FALSE, ]
  npc.choices <- unique(npc[, 2])
  choice.counts <- tabulate(match(npc[, 2], npc.choices))
  best.npc <- npc.choices[choice.counts == max(choice.counts)]
  median.period <- npc[npc[ ,2] == best.npc , 1] %>%
    median() %>%
    round()
  
  #Add plot illustrating prediction interval approach
  object <- Plot_nPCs(object,
                      npc = best.npc,
                      pcs = pcs,
                      period = median.period,
                      std.dev = std.dev)
  print(object@tools$Plot_nPCs)
  
  #Minimum number of principal components needed is two
  if(best.npc == 1){best.npc <- 2}
  
  
  return(list(object, best.npc))
}

Plot_nPCs <- function(object, npc, pcs, period, std.dev) {
  i <- period + npc
  temp.table <- as.data.frame(cbind(pcs[(i + 1 - period):i], std.dev[(i + 1 - period):i]))
  colnames(temp.table) <- c("pc", "sd")
  model <- lm(sd ~ pc, temp.table)
  prediction <- predict(model, data.frame(pc = i - period), interval = "prediction", level = 0.95)
  temp.table <- rbind(c(npc, prediction[1]), temp.table)
  npc_plot <- ElbowPlot(object, ndims = length(pcs)) +
    geom_smooth(color = "black", data = temp.table, aes (x = pc, y = sd),  se = F, method = "lm") +
    geom_segment(aes(x = npc, xend = npc, y = temp.table[1,2], yend = std.dev[npc]),
                 color = "red", linetype = "solid", size = 1.2) +
    annotate(geom = "text",  npc + 3 + length(pcs)/10,
             y = (std.dev[npc]+std.dev[(npc+1)])/2,
             color = "black",
             label = paste0("Max Prediction Error\n",
                            round(std.dev[npc] - prediction[1], digits = 2),
                            "  SDs absolute\n",
                            round((std.dev[npc] - prediction[1])/(prediction[3] - prediction[2]), digits = 2),
                            " SDs corrected"),
             size = 4)
  Tool(object) <- npc_plot
  return(object)
}

LoadEmbeddings <- function(object, standardize){
  npc <- object@tools$Choose_nPCs
  all.embeddings <- Embeddings(object)[, 1:npc]
  std_dev <- object[["pca"]]@stdev[1:npc]
  # Standardize principal components (if applicable)
  # subtracts PC mean from PC position (apply()),
  # then divides by std_dev (t(t(/std_dev)))
  if(standardize == TRUE){
    all.embeddings <- t(t(apply(all.embeddings, MARGIN = 2,
                                FUN = function(x) {x - mean(x)}))
                        /std_dev)
  }
  Tool(object) <- all.embeddings
  return(object)
}

library(segmented)
Segment_nPCs <- function(object, pc.use, pc.estimate){
  # Perform segmented regression
  pcs <- 1:pc.use
  std_dev <- object[["pca"]]@stdev[pcs]
  pc.data <- data.frame(pcs, std_dev)
  print("Choosing dimensionality.")
  model <- segmented(lm(formula = std_dev ~ pcs, data = pc.data),
                     seg.Z = ~pcs,
                     psi = pc.estimate,
                     control = seg.control(seed = 1))
  object <- Plot_nPCs_segmented(object, pc.use = pc.use, model = model)
  print(object@tools$Plot_nPCs_segmented)
  best.npc <- floor(model$psi[[2]])
  return(list(object, best.npc))
}

Plot_nPCs_segmented <- function(object, pc.use, model) {
  pcs <- 1:pc.use
  m <- slope(model)$pcs
  b <- intercept(model)$pcs
  psi <- model$psi[[2]]
  npc <- floor(psi)
  std.dev <- object[["pca"]]@stdev[pcs]
  
  npc.plot <- ggplot(as.data.frame(cbind(pcs, std.dev)), aes(x = pcs, y = std.dev)) +
    geom_point(size = 2.5) +
    geom_segment(aes(x = 1, y = (m[1, 1] * 1) + b[1, 1], xend = psi,
                     yend = (m[1, 1] * psi) + b[1, 1]), size = 1) +
    geom_segment(aes(x = psi, y = (m[2, 1] * psi) + b[2, 1], xend = pc.use,
                     yend = (m[2, 1] * pc.use) + b[2, 1]), size = 1) +
    geom_vline(xintercept = npc, color = "red3", size = 1, linetype = "longdash") +
    geom_text(aes(x = npc + (pc.use / 20), label = paste0("npc = ", npc),
                  y = max(std.dev, na.rm = TRUE) / 2), color = "red3", size = 5) +
    theme_cowplot(font_size = 16) +
    xlab("PC") +
    ylab("Standard Deviation")
  
  Tool(object) <- npc.plot
  return(object)
}


#To Do: Change progress bars to be of length iterations(/2)
#Also, this has been clustering using the standardized embeddings, which needs to be changed
ICVI_Simulation <- function(embeddings, clusters, name, cl.range, iterations){
  
  #Initialize Plot List
  plot.list <- list()
  
  #Initialize rho results table
  rho.table <- data.frame(matrix(ncol = 5, nrow = 3),
                          row.names = c("CH", "DB", "SI"))
  names(rho.table) <- c("RevLouvain", "Random", "KMeans", "KNN", "SNN")
  new.clusters <- clusters
  cl.num <- length(unique(clusters))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Reverse Louvain clustering
  
  ch.index <- sil.index <- db.index <- ari.index <- numeric(length = length(clusters))
  index.frame <- data.frame(ch.index, sil.index, db.index, ari.index)
  names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")
  new.clusters <- clusters
  
  print("Performing reverse Louvain Clustering")
  pb <- txtProgressBar(min = 1, max = 95, char = "+")
  for (i in 1:length(clusters)){
    new.clusters[i] <- sample(1:cl.num, 1)
    index.frame[i, 1:3] <- intCriteria(embeddings,
                                       part = as.integer(new.clusters),
                                       crit = c("calinski", "davies", "sil")) %>%
      unlist()
    index.frame[i, 4] <- ARI(new.clusters, clusters)
    setTxtProgressBar(pb, value = i/20)
  }
  close(pb)
  
  rho.table[ , 1] <- c(cor.test(x = index.frame$ARI,
                                y = index.frame$CH_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$ARI,
                                y = index.frame$DB_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$ARI,
                                y = index.frame$Sil_Std,
                                method = "spearman")$estimate)
  
  plot.list[[1]] <- ggplot(index.frame, aes(ARI, CH_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","Reverse Louvain Clustering, CH_STD vs ARI: Rho = ", 
                   rho.table[1,1] %>%
                     round(digits=4)))
  print(plot.list[[1]])
  
  
  plot.list[[2]] <- ggplot(index.frame, aes(ARI, DB_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","Reverse Louvain Clustering, DB_STD vs ARI: Rho = ", 
                   rho.table[2,1] %>%
                     round(digits=4)))
  print(plot.list[[2]])
  
  plot.list[[3]] <- ggplot(index.frame, aes(ARI, Sil_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","Reverse Louvain Clustering, Sil_STD vs ARI: Rho = ", 
                   rho.table[3,1] %>%
                     round(digits=4)))
  print(plot.list[[3]])
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Random Mutations
  
  ch.index <- sil.index <- db.index <- ari.index <- numeric(length = length(iterations))
  index.frame <- data.frame(ch.index, sil.index, db.index, ari.index)
  names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")
  
  print("Randomly mutating clusters")
  pb <- txtProgressBar(min = 1, max = 75, char = "+")
  for (i in 1:iterations){
    new.clusters <- clusters
    reassign <- sample(1:length(clusters), 1)
    idx <- sample(1:length(clusters), reassign)
    clust.num <- sample((cl.num-cl.range):(cl.num+cl.range), 1)
    new.clusters[idx] <- sample(1:clust.num, reassign, replace = T) 
    
    index.frame[i, 1:3] <- intCriteria(embeddings,
                               part = as.integer(new.clusters),
                               crit = c("calinski", "davies", "sil")) %>%
      unlist
    index.frame[i, 4] <- ARI(new.clusters, clusters)
    setTxtProgressBar(pb, i/2)
  }
  close(pb)
  
  rho.table[ , 2] <- c(cor.test(x = index.frame$ARI,
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
    ggtitle(paste0(name, ", ","Random Clustering, CH_STD vs ARI: Rho = ", 
                   rho.table[1,2] %>%
                     round(digits=4)))
  print(plot.list[[4]])
  
  plot.list[[5]] <- ggplot(index.frame, aes(ARI, DB_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","Random Clustering, DB_STD vs ARI: Rho = ", 
                   rho.table[2,2] %>%
                     round(digits=4)))
  print(plot.list[[5]])
  
  plot.list[[6]] <- ggplot(index.frame, aes(ARI, Sil_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","Random Clustering, Sil_STD vs ARI: Rho = ", 
                   rho.table[3,2] %>%
                     round(digits=4)))
  print(plot.list[[6]])
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #k Means simulation study
  
  #bsi is bull shit index, trying to adjust ARI by cluster number error
  ch.index <- sil.index <- db.index <- ari.index <- nmi  <- cluster.number <- cl.error <- bsi <- numeric(length = length(iterations))
  index.frame <- data.frame(ch.index, sil.index, db.index, ari.index, nmi, cluster.number, cl.error, bsi)
  names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI", "NMI", "Cluster_Number", "Cluster_Error", "AARI")
  
  #Basically just trying to make the ARI better
  index.frame[i, 8] <- index.frame[i, 4]*(1-index.frame[i, 7]/14)
  
  print("k-means simulation")
  pb <- txtProgressBar(min = 1, max = iterations/5, initial = 1, char = "+", style = 2)
  for (i in 1:iterations){
    #cl <- kmeans(embeddings, sample((cl.num-cl.range):(cl.num+cl.range), 1))
    cl <- kmeans(embeddings, sample(5:25, 1))
    index.frame[i, 1:3] <- intCriteria(embeddings,
                                       cl$cluster,
                                       c("Calinski", "Davies", "Sil")) %>%
      unlist()
    index.frame[i, 4] <- ARI(cl$cluster, clusters)
    index.frame[i, 5] <- NMI(cl$cluster, clusters)
    index.frame[i, 6:7] <- c(max(cl$cluster), abs(max(cl$cluster)-14))
    #Basically just trying to make the ARI better
    index.frame[i, 8] <- index.frame[i, 4]*(1-index.frame[i, 7]/14)
    
    setTxtProgressBar(pb, i/5)
  }
  close(pb)
  
  rho.table[ , 3] <- c(cor.test(x = index.frame$AARI,
                                y = index.frame$CH_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$AARI,
                                y = index.frame$DB_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$AARI,
                                y = index.frame$Sil_Std,
                                method = "spearman")$estimate)
  
  plot.list[[7]] <- ggplot(index.frame, aes(AARI, CH_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","K-Means Clustering, CH_STD vs AARI: Rho = ", 
                   rho.table[1,3] %>%
                     round(digits=4)))
  print(plot.list[[7]])
  
  
  plot.list[[8]] <- ggplot(index.frame, aes(AARI, DB_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","K-Means Clustering, DB_STD vs AARI: Rho = ", 
                   rho.table[2,3] %>%
                     round(digits=4)))
  print(plot.list[[8]])
  
  
  plot.list[[9]] <- ggplot(index.frame, aes(AARI, Sil_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","K-Means Clustering, Sil_STD vs AARI: Rho = ", 
                   rho.table[3,3] %>%
                     round(digits=4)))
  print(plot.list[[9]])
  
  ordered.frame <- index.frame[order(-index.frame$`Cluster Number`),]
  index.frame.new <- cbind(index.frame, abs(14-index.frame$`Cluster Number`))
  names(index.frame.new) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI", "Cluster_Number", "Cluster_Error")
  
  plot(index.frame$BSI, index.frame$CH_Std)
  plot(index.frame$Cluster_Error, index.frame$NMI)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #K Nearest Neighbors clustering
  
  ch.index <- sil.index <- db.index <- ari.index <- nmi  <- cluster.number <- cl.error <- bsi <- numeric(length = length(iterations))
  index.frame <- data.frame(ch.index, sil.index, db.index, ari.index, nmi, cluster.number, cl.error, bsi)
  names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI", "NMI", "Cluster_Number", "Cluster_Error", "BSI")
  

  print("Running KNN")
  pb <- txtProgressBar(min = 1, max = 50, initial = 1, char = "+", style = 2)
  k.list <- sample(10:240, iterations, replace = T)
  for (i in 1:iterations){
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
    com <- km$membership %>% as.integer()
    
    index.frame[i, 1:3] <- intCriteria(embeddings,
                                       com,
                                       c("Calinski", "Davies", "Sil")) %>%
      unlist()
    index.frame[i, 4] <- ARI(com, clusters)
    index.frame[i, 5] <- NMI(com, clusters)
    index.frame[i, 6:7] <- c(max(com), abs(max(com)-14))
    #Basically just trying to make the ARI better
    index.frame[i, 8] <- index.frame[i, 4]*(1-index.frame[i, 7]/14)
    setTxtProgressBar(pb, i/5)
  }
  close(pb)
  
  rho.table[ , 4] <- c(cor.test(x = index.frame$BSI,
                                y = index.frame$CH_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$BSI,
                                y = index.frame$DB_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$BSI,
                                y = index.frame$Sil_Std,
                                method = "spearman")$estimate)
  
  plot.list[[10]] <- ggplot(index.frame, aes(BSI, CH_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","KNN Clustering, CH_STD vs ARI: Rho = ", 
                   rho.table[1, 4] %>%
                     round(digits=4)))
  print(plot.list[[10]])
  
  plot.list[[11]] <- ggplot(index.frame, aes(BSI, DB_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","KNN Clustering, DB_STD vs ARI: Rho = ", 
                   rho.table[2, 4] %>%
                     round(digits=4)))
  print(plot.list[[11]])
  
  plot.list[[12]] <- ggplot(index.frame, aes(BSI, Sil_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","KNN Clustering, Sil_STD vs ARI: Rho = ", 
                   rho.table[3, 4] %>%
                     round(digits=4)))
  print(plot.list[[12]])
  plot(index.frame$BSI, index.frame$DB_Std)
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #SNN simulation study
  # only gives decent ARI values when cluster number is wayyyyy too high
  
  
  ch.index <- sil.index <- db.index <- ari.index <- numeric(length = length(iterations))
  names(index.frame) <- c("CH_Std", "DB_Std", "Sil_Std", "ARI")
  
  
  print("Running SNN")
  pb <- txtProgressBar(min = 1, max = 75, initial = 1, char = "+", style = 2)
  for (i in 1:iterations){
    k.samp <- sample(20:100, 1)
    knn <- kNN(embeddings, k = k.samp, search = 'dist')
    snn.test <- sNNclust(knn, k = sample(5:k.samp, 1), eps = 0, minPts = 1)
    index.frame[i, 1:3] <- intCriteria(embeddings,
                                       as.integer(snn.test$cluster),
                                       c("Calinski", "Davies", "Sil")) %>%
      unlist()
    index.frame[i, 4] <- ARI(snn.test$cluster, clusters)
    print(c(length(unique(snn.test$cluster)), index.frame[i,4]))
    setTxtProgressBar(pb, i/2)
  }
  close(pb)
  max(index.frame[4])
  
  rho.table[ , 5] <- c(cor.test(x = index.frame$ARI,
                                y = index.frame$CH_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$ARI,
                                y = index.frame$DB_Std,
                                method = "spearman")$estimate,
                       cor.test(x = index.frame$ARI,
                                y = index.frame$Sil_Std,
                                method = "spearman")$estimate)
  max(index.frame$ARI)
  
  plot.list[[13]] <- ggplot(index.frame, aes(ARI, CH_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","SNN Clustering, CH_STD vs ARI: Rho = ", 
                   rho.table[1,5] %>%
                     round(digits=4)))
  print(plot.list[[13]])
  
  plot.list[[14]] <- ggplot(index.frame, aes(ARI, DB_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","SNN Clustering, DB_STD vs ARI: Rho = ", 
                   rho.table[2,5] %>%
                     round(digits=4)))
  print(plot.list[[14]])
  
  plot.list[[15]] <- ggplot(index.frame, aes(ARI, Sil_Std)) + 
    geom_point() + 
    ggtitle(paste0(name, ", ","SNN Clustering, SI_STD vs ARI: Rho = ", 
                   rho.table[3,5] %>%
                     round(digits=4)))
  print(plot.list[[15]])
  
  return(list(rho.table, plot.list))
}
