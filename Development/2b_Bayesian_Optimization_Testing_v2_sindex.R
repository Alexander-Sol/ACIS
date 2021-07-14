# Bayesian optimization testing
library(Seurat)
library(clusterCrit)
library(GPfit)
library(lhs)
library(aricode)
library(purrr)
library(clusterSim)

getICVI <- function(object, nPC = NULL) {
  icvi <- suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:nPC],
                                            method = "euclidean"),
                                   cl = as.integer(Idents(object))))
  if(is.nan(icvi)) {icvi <- 0}
  #lmao because we're minimizing here
  return(icvi)
}

AnalyzeClusters <- function(k.param, resolution, object, npcs = 10){
  # Construct SNN graph for ("K.param", "Resolution") = (i, j)
  # These need to be stored. As it is, the nn graph is being calculated multiple times for identical k.params
  temp.object <- FindNeighbors(object,
                               k.param = k.param,
                               compute.SNN = TRUE,
                               prune.SNN = 1/15,
                               nn.method = "rann",
                               nn.eps = 0.0,
                               verbose = FALSE,
                               reduction = "pca",
                               dims = 1:npcs) %>%
    
    # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
    FindClusters(modularity.fxn = 1,
                 resolution = resolution,
                 algorithm = 1,
                 group.singletons = TRUE,
                 verbose = FALSE)
  
  # Determine and store cluster number for ("K.param", "Resolution") = (i, j)
  clusters <- as.integer(Idents(temp.object))
  K <- length(levels(temp.object))
  
  # Calculate and store internal clustering validation index for ("K.param", "Resolution") = (i, j)
  if(K == 1){
    ind <- 0.2
  } else{
    #ind <- clusterCrit::intCriteria(object@reductions$pca@cell.embeddings[ , 1:npcs],
     #                               part = clusters,
      #                              crit = icvi)[[1]]
    ind <- getICVI(temp.object, nPC = npcs)
  }
  cat("+")
  return(ind)
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

scaleParams <- function(unit.value, param.range) {
  unit.value * (param.range[2] - param.range[1]) + param.range[1]
}

getARI <- function(object, k.param, resolution, npcs) {
  temp.object <- FindNeighbors(object,
                               k.param = k.param,
                               compute.SNN = TRUE,
                               prune.SNN = 1/15,
                               nn.method = "rann",
                               nn.eps = 0.0,
                               verbose = FALSE,
                               reduction = "pca",
                               dims = 1:npcs) %>%
    
    # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
    FindClusters(modularity.fxn = 1,
                 resolution = resolution,
                 algorithm = 1,
                 group.singletons = TRUE,
                 verbose = FALSE)
  
  # Determine and store cluster number for ("K.param", "Resolution") = (i, j)
  clusters <- as.integer(Idents(temp.object))
  K <- length(levels(temp.object))
  
  best.ari <- ARI(clusters, object$annotation)
  return(c(K, best.ari))
}

subCluster <- function(clustered.object,
                       cluster.ident,
                       npcs = 10,
                       force.split = F,
                       k.param.range = c(20, 500),
                       res.range = c(0.05, 2.4)) {
  sc.object <- subset(clustered.object, idents = cluster.ident)
  sc.k.param.range <- k.param.range
  if (k.param.range[2] > ncol(sc.object)/2) {sc.k.param.range[2] <- round(ncol(sc.object)/2) }
  sc.object <- bayesianClustering(object = sc.object,
                                  k.param.range = sc.k.param.range,
                                  res.range = res.range,
                                  n.starts = 25,
                                  n.iterations = 10,
                                  epsilon = 0.01,
                                  predicted.npcs = npcs)
  
  clustered.object <- adjustIdents(object = clustered.object,
                                   temp.object = sc.object,
                                   npc = npcs,
                                   force.split = force.split)
  return(clustered.object)
  
}

as.int.factor <- function(x) {as.integer(levels(x))[x]}

adjustIdents <- function(object, temp.object, npc, force.split = F) {
  embeds <- Embeddings(object)[ , 1:npc]
  new.idents <- as.int.factor(temp.object@active.ident) + 
    1 + 
    max(as.int.factor(object@active.ident))
  names(new.idents) <- Cells(temp.object)
  combined.new <- as.int.factor(object@active.ident)
  names(combined.new) <- Cells(object)
  combined.new[names(new.idents)] <- new.idents
  original.SI <- suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:npc],
                                                   method = "euclidean"),
                                          cl = as.integer(Idents(object))))
  # original.SI <- intCriteria(embeds, as.int.factor(object@active.ident), crit = "sil")[[1]]
  new.SI <-  suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:npc],
                                               method = "euclidean"),
                                      cl = combined.new %>% as.factor %>% as.numeric))
  # new.SI <- intCriteria(embeds, as.integer(combined.new), "sil")[[1]]
  if(original.SI > new.SI) {
    print(paste(object@project.name, "REJECTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
    if(force.split) {object@active.ident <- as.factor(combined.new)}
  } else if (original.SI <= new.SI) {
    object@active.ident <- as.factor(combined.new)
    print(paste(object@project.name, "ACCEPTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
  }
  return(object)
}

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



#Epsilon is a hyper parameter that controls exploration/exploitation balance.
#Higher epsilon = more exploration
bayesianClustering <- function(object,
                               k.param.range,
                               res.range,
                               n.starts = 16,
                               n.iterations = 10,
                               epsilon = 0.01, 
                               predicted.npcs = NA) {
  set.seed(Sys.time())
  unit.starts <- optimumLHS(n.starts,2) %>% data.frame() %>% setNames(c("x", "y"))
  plot(unit.starts$x, unit.starts$y)
  
  scaled.k <- scaleParams(unit.starts$x, k.param.range) %>% round()
  scaled.res <- scaleParams(unit.starts$y, res.range)
  if(is.na(predicted.npcs)) {predicted.npcs <- Predict_nPCs(object)}
  grid.space <- expand.grid(x = seq(0, 1, 0.005),
                            y = seq(0, 1, 0.01))
  
  #Need to create a progress bar for this, probably by smuggling it into AnalyzeClusters
  icvi.values <- map2_dbl(
    scaled.k, scaled.res,
    function(k, res) {
      AnalyzeClusters(k.param = k, resolution = res,
                      object = object,
                      npcs = predicted.npcs)
    }
  )  
  
  pb <- txtProgressBar(min = 1, max = n.iterations, style = 3)
  for (i in 1:n.iterations) {
    setTxtProgressBar(pb, i)
    y.max <- max(icvi.values)
    #Change this so that corr can be configured via some argument
    gpModel <- GP_fit(X = unit.starts, Y = icvi.values,
                      corr = list(type = "matern",
                                  nu = 1/2))
    plot(gpModel, range = c(0, 1), resolution = 50, surf_check = T, response = T)
    predicted.return <- predict(gpModel, xnew = grid.space)
    expected_improvement <- map2_dbl(
      predicted.return$Y_hat, sqrt(predicted.return$MSE),
      function(m, s) {
        if(s == 0) {return(0)}
        gamma <- (m - y.max - epsilon) / s
        phi <- pnorm(gamma)
        return(s * (gamma * phi + dnorm(gamma)))
      }
    )
    if(max(expected_improvement, na.rm = T) <= 0) {break}
    new.params <- grid.space[which.max(expected_improvement), ]
    unit.starts <- rbind(new.params, unit.starts)
    new.k <- scaleParams(new.params[[1]], k.param.range) %>% round()
    new.res <- scaleParams(new.params[[2]], res.range)
    new.icvi <- AnalyzeClusters(k.param = new.k, resolution = new.res,
                                object = object,
                                npcs = predicted.npcs)
    icvi.values <- append(new.icvi, icvi.values)
  }
  
  best.params <- unit.starts[which.max(icvi.values), ]
  'best.ari <- getARI(object, 
                k.param = scaleParams(best.params[[1]], k.param.range) %>% round(),
                resolution = scaleParams(best.params[[2]], res.range),
                npcs = predicted.npcs)
  print(c(max(icvi.values), which.max(icvi.values), best.ari))'
  'return(list(top.icvi = max(icvi.values),
              position = which.max(icvi.values),
              ari.at.top = best.ari,
              best.params = c(scaleParams(best.params[[1]], k.param.range) %>% round(),
                              scaleParams(best.params[[2]], res.range))
              )
         )'
  clustered.object <- FindNeighbors(object,
                                    dims = 1:predicted.npcs,
                                    k.param = scaleParams(best.params[[1]], k.param.range) %>% round(),
                                    verbose = F) %>%
    FindClusters(resolution =  scaleParams(best.params[[2]], res.range),
                 verbose = F)
  print(paste("Optimum K =", scaleParams(best.params[[1]], k.param.range) %>% round(),
              "Optimum Res =", scaleParams(best.params[[2]], res.range)))
  return(clustered.object)
}

bayesSubCluster <- function(object,
                            k.param.range,
                            res.range,
                            n.starts = 16,
                            n.iterations = 10,
                            epsilon = 0.01,
                            sc.iter = 2,
                            npcs = NA) {
  
  start.time <- Sys.time()
  
  if(is.na(npcs)) {npcs <- Predict_nPCs(object)}
  print(paste(npcs, "PCs predicted"))
  
  clustered.object <- bayesianClustering(object = object,
                                    k.param.range = k.param.range,
                                    res.range = res.range,
                                    n.starts = n.starts,
                                    n.iterations = n.iterations,
                                    epsilon = 0.01,
                                    predicted.npcs = npcs)
  
  clustered.SI <- getICVI(clustered.object, nPC = npcs)
  
  clusterSilhouettes <- getSilScores(clustered.object, npcs)
  
  if(sc.iter > 0) {
    for (i in 1:sc.iter){
      #print(clusterSilhouettes)
      sc.object <- subset(clustered.object, idents = clusterSilhouettes[i, 1])
      sc.k.param.range <- k.param.range
      if (k.param.range[2] > ncol(sc.object)/2) {sc.k.param.range[2] <- round(ncol(sc.object)/2) }
      sc.object <- bayesianClustering(object = sc.object,
                                     k.param.range = sc.k.param.range,
                                     res.range = res.range,
                                     n.starts = 25,
                                     n.iterations = 5,
                                     epsilon = 0.01,
                                     predicted.npcs = npcs)
      
      clustered.object <- adjustIdents(object = clustered.object,
                                       temp.object = sc.object,
                                       npc = npcs)
    }
  }
  
  sc.SI <- getICVI(clustered.object, nPC = npcs)
  
  end.time <- Sys.time() - start.time
  print(paste("Time elapsed:", end.time))
  
  #return(c(clustered.SI, clustered.ARI,
  #         sc.SI, sc.ARI,
  #         end.time))
  return(clustered.object)
}


'bsc.results <- lapply(datasets,
                      bayesSubCluster,
                      k.param.range = k.param.range,
                      res.range = res.range,
                      n.starts = 25,
                      n.iterations = 25,
                      epsilon = 0.01)

bsc.tab <- do.call(rbind, bsc.results)
write.table(bsc.tab, file = \'clipboard\', sep = \'\t\')
'




