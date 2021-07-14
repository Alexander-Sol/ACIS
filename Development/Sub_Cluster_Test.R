#Subclustering Test
#These parameters maximize SI using a exp.distributed param space and predicted nPCs

#         Resolution   K.Param     nPC
#Zeisel   0.2690869       128       12
#Panc.1   0.22627417       76        5
#Panc.2   0.04756828       13        8
#PR       0.3805463        91        6
#MR       0.04756828        5        9
library(Seurat)
library(aricode)
library(stats)
library(clusterCrit)

as.int.factor <- function(x) {as.integer(levels(x))[x]}

#Calculate Silhouette scores per cluster
getSilScores <- function(object, npc) {
  embeds <- Embeddings(object)[ , 1:npc]
  embed.table <- cbind(embeds,
                       as.integer(levels(object@active.ident))[object@active.ident])
  sil.scores <- cluster::silhouette(as.integer(embed.table[ ,ncol(embed.table)]),
                                    dist(embed.table[, 1:(ncol(embed.table-1))]))
  sil <- data.frame(cluster = sil.scores[,1], si_width = sil.scores[,3])
  scores.by.cluster <- aggregate(sil$si_width, list(sil$cluster), mean) %>%
    arrange(x)
  names(scores.by.cluster) <- c("Seurat_cluster", "avg.SI.score")
  return(scores.by.cluster)
}

# SubClustering across Param Space
subCluster <- function(object, 
                       npc,
                       si.max.ident){
  temp.object <- subset(x = object, idents = si.max.ident)
  
  #Construct param space
  max.k <- (ncol(temp.object)/2) %>%
    round()
  k.param <- 2^(7:32/4)
  k.param <- k.param[k.param <= max.k]
  k.param.space <- k.param %>%
    round() %>%
    rep(each = length(k.param))
  res.space <- rep(k.param[k.param <= max.k]/ 100, length(k.param))
  iterations <- length(k.param)^2
  embeds <- Embeddings(temp.object)[,1:npc]
  
  
  #Create table to store results
  index.frame <- data.frame(matrix(nrow = iterations, ncol = 4))
  names(index.frame) <- c("Sil", 
                          "Cluster_Number",
                          "K.Param", "Resolution")
  
  #Iterate through datasets 
  start.time <- Sys.time()
  pb <- txtProgressBar(min = 1, max = length(k.param), initial = 1, char = "+", style = 3)
  for (i in 1:iterations){
    temp.object <- FindNeighbors(temp.object,
                                 k.param = k.param.space[i],
                                 compute.SNN = TRUE,
                                 prune.SNN = 1/15,
                                 nn.method = "annoy",
                                 nn.eps = 0.0,
                                 verbose = FALSE,
                                 force.recalc = F,
                                 reduction = "pca",
                                 dims = 1:npc) %>%
      FindClusters(modularity.fxn = 1,
                   resolution = res.space[i],
                   algorithm = 1,
                   group.singletons = TRUE,
                   verbose = FALSE)
    
    # Determine and store cluster number for ("K.param", "Resolution") = (i, j)
    temp.clusters <- as.integer(Idents(temp.object))
    K <- length(levels(temp.object))
    
    # Calculate and store internal clustering validation index for ("K.param", "Resolution") = (i, j)
    if(K == 1){
      index.frame[i, 1] <- NA
    } else{
      index.frame[i, 1] <- intCriteria(embeds,
                                       temp.clusters,
                                       "Sil")[[1]]
    }
    index.frame[i, 2:4] <- c(K, k.param.space[i], res.space[i])
    setTxtProgressBar(pb, i/length(k.param))
  }
  close(pb)
  end.time <- Sys.time() - start.time
  print(end.time)
  
  #Sanitize to remove NAs
  index.frame <- index.frame[!is.na(index.frame$Sil),]
  
  temp.object <- FindNeighbors(temp.object,
                               k.param = median(index.frame[index.frame$Sil == max(index.frame$Sil),
                                                            "K.Param"]),
                               dims = 1:npc) %>%
    FindClusters(resolution = median(index.frame[index.frame$Sil == max(index.frame$Sil),
                                                 "Resolution"]), 
                 verbose = F)
  
  object <- adjustIdents(object = object, 
                         temp.object = temp.object,
                         npc = npc)
  return(object)
}

#Test new clusters, accept or reject
adjustIdents <- function(object, temp.object, npc) {
  embeds <- Embeddings(object)[ , 1:npc]
  new.idents <- as.int.factor(temp.object@active.ident) + 
    1 + 
    max(as.int.factor(object@active.ident))
  names(new.idents) <- Cells(temp.object)
  combined.new <- as.int.factor(object@active.ident)
  names(combined.new) <- Cells(object)
  combined.new[names(new.idents)] <- new.idents
  original.SI <- intCriteria(embeds, as.int.factor(object@active.ident), crit = "dav")[[1]]
  new.SI <- intCriteria(embeds, as.integer(combined.new), "dav")[[1]]
  if(original.SI < new.SI) {
    print(paste(object@project.name, "REJECTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
  } else if (original.SI >= new.SI) {
    object@active.ident <- as.factor(combined.new)
    print(paste(object@project.name, "ACCEPTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
    accepted <<- accepted + 1
  }
  return(object)
}

datasets <- list(Zeisel = zeisel,
                 Panc.1 = panc.1,
                 Panc.2 = panc.2,
                 PR =pr,
                 MR = mr)
sc.datasets <- list()

optimal.params <- data.frame(dataset = c("Zeisel", "Panc.1", "Panc.2", "PR", "MR"),
                             resolution = c(0.2690869,0.22627417,0.04756828,0.3805463,0.04756828),
                             k.param = c(128, 76, 13, 91, 5),
                             npc = c(12, 5, 8, 6, 9),
                             row.names = c("Zeisel", "Panc.1", "Panc.2", "PR", "MR"))

sc.results <- data.frame(matrix(nrow = 5, ncol = 6))
names(sc.results) <- c("Datasets", "InitialCluster#", "FinalCluster#",
                       "InitialARI", "FinalARI", "SCsAccepted")
#SubCluster All Datasets
for (i in 1:5){
  object <- datasets[[i]]
  sc.results[i, "Datasets"] <- data.name <- names(datasets)[i]
  npc <- optimal.params[data.name, "npc"]
  k.param.full <- optimal.params[data.name, "k.param"]
  res.full <- optimal.params[data.name, "resolution"]
  accepted <- 0
  
  object <- FindNeighbors(object, dims = 1:npc, k.param = k.param.full) %>%
    FindClusters(resolution = res.full,
                 verbose = F)
  sc.results[i, "InitialCluster#"] <- length(levels(object))
  sc.results[i, "InitialARI"] <- ARI(as.int.factor(object@active.ident),
                                     as.factor(object$annotation) %>%
                                       as.integer())
  
  scores <- getSilScores(object = object, npc = npc)
  for ( j in scores[1:2, 1] ) {
    object <- subCluster(object = object,
                         npc = npc, 
                         si.max.ident = j)
  }
  
  sc.results[i, "FinalCluster#"] <- length(levels(object))
  sc.results[i, "FinalARI"] <- ARI(as.int.factor(object@active.ident),
                                   as.factor(object$annotation) %>%
                                     as.integer())
  sc.results[i, "SCsAccepted"] <- accepted
  sc.datasets[[i]] <- object
}

#Modified pc prediction method
#Relatively stable for different values of pc.use, although it infrequently selects
#extremely high elbows (~37)
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
  #median.period <- round(pc.use - best.npc/2)
  return(best.npc)
}

#Create Plots
for (i in 1:5){
  object <- sc.datasets[[i]]
  data.name <- names(datasets)[i]
  npc <- optimal.params[data.name, "npc"]
  k.param.full <- optimal.params[data.name, "k.param"]
  object <- RunUMAP(object, dims = 1:npc, n.neighbors = k.param.full)
  save.plot(DimPlot(object, pt.size = 1.4, group.by = "annotation", label = T) + sparse(),
            name = paste0(data.name, "_Annotation"))
  save.plot(DimPlot(object, pt.size = 1.4, label = T) + sparse(),
            name = paste0(data.name, "_DB_SubClusters"))
}

save.plot <- function(plot, name, folder = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/SubCluster_UMAPs/"){
  ggsave(filename = paste0(folder, name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
  dev.off()
}



lapply(datasets, Predict_nPCs)






