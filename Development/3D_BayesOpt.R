# 3D Bayesian Optimization Testing
library(Seurat)
library(clusterCrit)
library(GPfit)
library(lhs)
library(aricode)
library(purrr)
library(progress)
library(wesanderson)

getICVI <- function(object, nPC = NULL) {
  icvi <- suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:nPC],
                                            method = "euclidean"),
                                   cl = as.integer(Idents(object))))
  if(is.nan(icvi)) {icvi <- 0}
  #lmao because we're minimizing here
  return(icvi)
}

AnalyzeClusters <- function(k.param, resolution, object, npcs = 10, cluster.pc = 10){
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
                               dims = 1:cluster.pc) %>%
    
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

pbAnalyzeClusters <- function(k.param, resolution, object, npcs = 10, cluster.pc = 10){
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
                               dims = 1:cluster.pc) %>%
    
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

getARI3d <- function(object, k.param, resolution, cluster.pc) {
  temp.object <- FindNeighbors(object,
                               k.param = k.param,
                               compute.SNN = TRUE,
                               prune.SNN = 1/15,
                               nn.method = "rann",
                               nn.eps = 0.0,
                               verbose = FALSE,
                               reduction = "pca",
                               dims = 1:cluster.pc) %>%
    
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
bayesianClustering3D <- function(object,
                               k.param.range = c(3, 340),
                               res.range = c(0.03, 2.4),
                               pc.range = c(2,25),
                               n.starts = 16,
                               n.iterations = 10,
                               epsilon = 0.01, 
                               predicted.npcs = NA) {
  set.seed(Sys.time())
  unit.starts <- maximinLHS(n.starts, 3) %>% data.frame() %>% setNames(c("x", "y", "z"))
  
  scaled.k <- scaleParams(unit.starts$x, k.param.range) %>% round()
  scaled.res <- scaleParams(unit.starts$y, res.range)
  scaled.pc <- scaleParams(unit.starts$z, pc.range)
  if(is.na(predicted.npcs)) {predicted.npcs <- Predict_nPCs(object)}
  grid.space <- expand.grid(x = seq(0, 1, 0.05),
                            y = seq(0, 1, 0.05),
                            z = seq(0, 1, 0.04))
  
  print(paste0("Calculating Silhouette scores at ", n.starts, " initial points"))
  #using global assignment is probably a bad idea here, but idk how else to get this working
  pb.initial <<- progress::progress_bar$new(total = n.starts)
  #Need to create a progress bar for this, probably by smuggling it into AnalyzeClusters
  icvi.values <- pmap_dbl(
    list(scaled.k, scaled.res, scaled.pc),
    function(k, res, pcs) {
      pbAnalyzeClusters(k.param = k,
                         resolution = res,
                         object = object,
                         npcs = predicted.npcs,
                         cluster.pc = pcs)
    }
  )  
  
  print(paste0("Performing Bayesian Optimization for ", n.iterations, " iterations"))
  pb <- txtProgressBar(min = 0, max = n.iterations, style = 3)
  for (i in 1:n.iterations) {
    y.max <- max(icvi.values)
    #Change this so that corr can be configured via some argument
    if(max(unit.starts > 1)) {print(paste("1 exceeded", max(unit.starts), i))}
    if(max(unit.starts < 0)) {print(paste("0 depleted", max(unit.starts), i))}
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
    if(max(expected_improvement, na.rm = T) <= 0) {break}
    new.params <- grid.space[which.max(expected_improvement), ]
    unit.starts <- rbind(new.params, unit.starts)
    new.k <- scaleParams(new.params[[1]], k.param.range) %>% round()
    new.res <- scaleParams(new.params[[2]], res.range)
    new.pcs <- scaleParams(new.params[[3]], pc.range) %>% round()
    new.icvi <- AnalyzeClusters(k.param = new.k, resolution = new.res,
                                object = object,
                                npcs = predicted.npcs,
                                cluster.pc = new.pcs)
    icvi.values <- append(new.icvi, icvi.values)
    setTxtProgressBar(pb, i)
  }
  
  best.params <- unit.starts[which.max(icvi.values), ]
  best.ari <- getARI3d(object, 
                k.param = scaleParams(best.params[[1]], k.param.range) %>% round(),
                resolution = scaleParams(best.params[[2]], res.range),
                cluster.pc = scaleParams(best.params[[3]], pc.range) %>% round())
  return(list(top.icvi = max(icvi.values),
              ari.at.top = best.ari[2],
              cluster.number = best.ari[1],
              best.k = scaleParams(best.params[[1]], k.param.range) %>% round(),
              best.res = scaleParams(best.params[[2]], res.range),
              best.npc = scaleParams(best.params[[3]], pc.range) %>% round(),
              runs.till.best = 50 - which.max(icvi.values)
             )
         )
 ' clustered.object <- FindNeighbors(object,
                                    dims = 1:predicted.npcs,
                                    k.param = scaleParams(best.params[[1]], k.param.range) %>% round(),
                                    verbose = F) %>%
    FindClusters(resolution =  scaleParams(best.params[[2]], res.range),
                 verbose = F)
  return(clustered.object)'
}



bayes3D_results <- vector("list", 25)
for (i in 1:25) {
  bayes3D_results[[i]] <- lapply(datasets, bayesianClustering3D, n.starts = 25, n.iterations = 25)
  print(paste("Completed 3D Test iteration:", i))
}

saveRDS(bayes3D_results, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/25_bayes3D_results_forReal.rds")

source("C:/Users/asolivai/Desktop/R_Files/Scripts_Mac/AutoClustR Development/bayesPC_Test.R")

b3d <- readRDS("C:/Users/asolivai/Desktop/R_Files/AutoClustR/25_bayes3D_results_forReal.rds")


# Workup the results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save.plot <- function(plot, name,
                      folder = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/3D_Testing/"){
  ggsave(filename = paste0(folder, name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 4,
         width = 8,
         dpi = 300)
  dev.off()
}


for (dataset in names(datasets)){
  z.res <- res.list.3d[[dataset]]
  z.res.table <- data.frame(npc = 2:25,
                            count = rep(0, 24),
                            avg.ari = rep(0, 24),
                            lab = rep("", 24))
  z.res.table[z.res.table$npc %in% names(table(z.res$best.npc)), "count"] <- table(z.res$best.npc)
  ari.by.npc <- z.res %>% group_by(best.npc) %>% summarise(avg.ari = mean(ari.at.top))
  z.res.table[z.res.table$npc %in% ari.by.npc$best.npc, "avg.ari"] <- ari.by.npc$avg.ari
  z.res.table[z.res.table$avg.ari != 0, "lab"] <- ari.by.npc$avg.ari %>% round(3)
  z.res.table[z.res.table$avg.ari == 0, "avg.ari"] <- min(ari.by.npc$avg.ari)
  
  plot <- ggplot(data = z.res.table, aes(x = npc, y = count, fill = avg.ari)) +
    geom_bar(stat = "identity") +
    scale_fill_continuous(low = "gold1", high = "darkred") +
    geom_text(aes(label = lab), angle = 40, vjust = 2.0, hjust = 0.5, size = 4 ) +
    ggtitle(label = paste0(dataset, ": ", Predict_nPCs(datasets[[dataset]]), " PCs Predicted"))

  save.plot(plot, dataset)
}


#   Detritus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


res.list.3d <- setNames(vector(mode = "list", 5), names(datasets))
colors.list <- setNames(vector(mode = "list", 5), names(datasets))
for(dataset in names(datasets)){
  temp <- lapply(b3d, pluck, dataset)
  res.list.3d[[dataset]] <- do.call(rbind, temp) %>% as.data.frame()
  res.list.3d[[dataset]]$best.npc <- as.integer(res.list.3d[[dataset]]$best.npc)
  res.list.3d[[dataset]]$ari.at.top <- as.numeric(res.list.3d[[dataset]]$ari.at.top)
  
  sort.temp <- res.list.3d[[dataset]] %>%
    group_by(best.npc) %>%
    summarise(avg.ari = mean(ari.at.top))
  colors.list[[dataset]] <- data.frame(best.npc = 2:25,
                                       avg.ari = min(sort.temp$avg.ari))
  colors.list[[dataset]][colors.list[[dataset]]$best.npc %in% sort.temp$best.npc, 2] <- sort.temp$avg.ari
  colors.list[[dataset]] <- colors.list[[dataset]][order(colors.list[[dataset]]$avg.ari),]
  colors.list[[dataset]]$colors <- wes_palette("Zissou1", 24, type = "continuous") 
  colors.list[[dataset]] <- colors.list[[dataset]][order(colors.list[[dataset]]$best.npc),]
}

for(dataset in names(datasets)) {
  plot <- ggplot(data = res.list.3d[[dataset]], aes(x = best.npc)) +
    geom_histogram(bins = 24, fill = colors.list[[dataset]]$colors) +
    ggtitle(label = paste0(dataset, ": ", Predict_nPCs(datasets[[dataset]]), " PCs Predicted"))
  save.plot(plot, dataset)
}



test <- res.list.3d$Zeisel %>% 
  group_by(best.npc) %>% 
  summarise(avg.ari = mean(ari.at.top)) 
hist.colors <- data.frame(best.npc = 2:25,
                          avg.ari = min(test$avg.ari))
hist.colors[hist.colors$best.npc %in% test$best.npc, 2] <- test$avg.ari
hist.sort <- hist.colors[order(hist.colors$avg.ari),]
hist.sort$colors <- wes_palette("Zissou1", 24, type = "continuous") 
hist.colors <- hist.sort[order(hist.sort$best.npc),]

ggplot(data = res.list.3d[["PR"]], aes(x = best.npc)) +
  geom_histogram(bins = 24, fill = colors.list[["PR"]]$colors) 



res.list.3d[["Zeisel"]]$col.scale <- res.list.3d$Zeisel$ari.at.top
name.list <- setNames(hist.colors$best.npc, as.list(hist.colors$avg.ari))
levels(res.list.3d[["Zeisel"]]$col.scale) <- name.list





top.ari <- lapply(bayes3D_results, pluck, "ari.at.top") %>% unlist()
top.npc <- lapply(bayes3D_results, pluck, "best.params") %>% lapply(pluck, 3) %>% unlist()
top.sil <- lapply(bayes3D_results, pluck, "top.icvi") %>% unlist()

plot(y = top.sil, x = top.npc)




for(i in 1:10){
  bayesianClustering3D(zeisel, n.starts = 5, n.iterations = 5)[c("ari.at.top", "top.icvi", "best.npc")] %>%
    print()
}

