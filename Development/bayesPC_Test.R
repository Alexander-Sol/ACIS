# Bayes iterating through PCs
library(lhs)
getICVI <- function(object, nPC = NULL) {
  icvi <- suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:nPC],
                                            method = "euclidean"),
                                   cl = as.integer(Idents(object))))
  if(is.nan(icvi)) {icvi <- 0}
  return(icvi)
}


AnalyzeClustersPCtest <- function(k.param, resolution, object, npcs = 10){
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
  return(ind)
}

#This would be better as functional programming. basically just pass the argument and then the function to be run after the pb is advanced
pbAnalyzeClustersPCtest <- function(k.param, resolution, object, npcs = 10){
  pb.initial$tick()
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
  return(ind)
}

scaleParams <- function(unit.value, param.range) {
  unit.value * (param.range[2] - param.range[1]) + param.range[1]
}

getARIPCtest <- function(object, k.param, resolution, npcs) {
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

#Epsilon is a hyper parameter that controls exploration/exploitation balance.
#Higher epsilon = more exploration
bayesianClusteringPCtest <- function(object,
                                 k.param.range = c(3, 340),
                                 res.range = c(0.03, 2.4),
                                 n.starts = 25,
                                 n.iterations = 25,
                                 epsilon = 0.01, 
                                 npcs = 2) {
  set.seed(Sys.time())
  unit.starts <- maximinLHS(n.starts, 2) %>% data.frame() %>% setNames(c("x", "y"))
  scaled.k <- scaleParams(unit.starts$x, k.param.range) %>% round()
  scaled.res <- scaleParams(unit.starts$y, res.range)

  grid.space <- expand.grid(x = seq(0, 1, 0.005),
                            y = seq(0, 1, 0.02))
  
  print(paste0("Calculating Silhouette scores at ", n.starts, " initial points"))
  #using global assignment is probably a bad idea here, but idk how else to get this working
  pb.initial <<- progress::progress_bar$new(total = n.starts)
  #Need to create a progress bar for this, probably by smuggling it into AnalyzeClusters
  icvi.values <- map2_dbl(
    scaled.k, scaled.res,
    function(k, res) {
      AnalyzeClustersPCtest(k.param = k, resolution = res,
                            object = object,
                            npcs = npcs)
    }
  )
  
  print(paste0("Performing Bayesian Optimization for ", n.iterations, " iterations"))
  pb <- txtProgressBar(min = 0, max = n.iterations, style = 3)
  for (i in 1:n.iterations) {
    y.max <- max(icvi.values)
    #Change this so that corr can be configured via some argument
    print(min(unit.starts))
    print(max(unit.starts))
    gpModel <- GP_fit(X = unit.starts, Y = icvi.values,
                      corr = list(type = "matern",
                                  nu = 1/2))
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
    #This should probably be changed to some threshold value
    # if(max(expected_improvement, na.rm = T) <= 0) {break}
    new.params <- grid.space[which.max(expected_improvement), ]
    unit.starts <- rbind(new.params, unit.starts)
    new.k <- scaleParams(new.params[[1]], k.param.range) %>% round()
    new.res <- scaleParams(new.params[[2]], res.range)
    new.icvi <- AnalyzeClustersPCtest(k.param = new.k, resolution = new.res,
                                      object = object,
                                      npcs = npcs)
    icvi.values <- append(new.icvi, icvi.values)
    setTxtProgressBar(pb, i)
  }
  
  best.params <- unit.starts[which.max(icvi.values), ]
  best.ari <- getARIPCtest(object, 
                       k.param = scaleParams(best.params[[1]], k.param.range) %>% round(),
                       resolution = scaleParams(best.params[[2]], res.range),
                       npcs = npcs)
  return(list(top.icvi = max(icvi.values),
              ari.at.top = best.ari[2],
              cluster.number = best.ari[1],
              best.k = scaleParams(best.params[[1]], k.param.range) %>% round(),
              best.res = scaleParams(best.params[[2]], res.range),
              npc = npcs,
              runs.till.best = 50 - which.max(icvi.values)
              )
  )
}

npcTestResults <- vector("list", 9)
for(npc in 17:25) {
  bayes3D_results <- vector("list", 10)
  for (i in 1:10) {
    bayes3D_results[[i]] <- lapply(datasets, bayesianClusteringPCtest,
                                   n.starts = 25,
                                   n.iterations = 25,
                                   npcs = npc)
  }
  npcTestResults[[npc-1]] <- bayes3D_results
  print(paste("Highest PC completed:", npc))
}


saveRDS(npcTestResults, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/Bayes_npc_test_results_17to25.rds")





# Clean Results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


npc.res <- readRDS("C:/Users/asolivai/Desktop/R_Files/AutoClustR/Bayes_npc_test_results_2to16.rds")
npc.res[16:24] <- npcTestResults[16:24]

saveRDS(npc.res, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/Bayes_npc_test_results_TrainingFull.rds")


ari.list <- list(length = 5)
for (i in 1:5) {ari.list[[i]] <- data.frame(matrix(0, nrow = 10, ncol = 24,
                                     dimnames = list(NULL, paste("nPC", 2:25))))}
names(ari.list) <- names(datasets)

for(dataset in names(datasets)){
  for (i in 1:24){
    ari.list[[dataset]][ , i] <- lapply(npc.res[[i]], pluck, dataset)  %>% 
      lapply(pluck, "ari.at.top") %>% unlist()
  }
}  


for(dataset in names(datasets)){
  d.stack <- stack(ari.list[[dataset]])
  d.stack$nPC <- (as.integer(d.stack[["ind"]]) + 1) %>% as.factor()
  plot <- ggplot(data = d.stack) + geom_boxplot(aes(x = nPC, y = values)) +
    xlab("Number of Principal Components") +
    ylab("ARI") +
    ggtitle(paste0(dataset, " dataset; ",
                   Predict_nPCs(datasets[[dataset]]),
                   " predicted nPCs"))
  save.plot(plot, dataset)
}


save.plot <- function(plot, name,
                      folder = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/PC_Testing/"){
  ggsave(filename = paste0(folder, name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 4,
         width = 8,
         dpi = 300)
  dev.off()
}
  

  
  
  
  