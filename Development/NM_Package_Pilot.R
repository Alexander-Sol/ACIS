library(neldermead)
library(clusterSim)
library(lhs)
library(aricode)
library(clusterCrit)

#Object is passed to cost function as an element of fmsfundata
scoreCluster <- function(x = null, index = NULL, fmsfundata = NULL) {
  #Also, despite being theoretically only searching the space [0 , 1], nm sometimes passes negative values for x[1]
  #the hacky work around is to manually set it to be positive, but longterm we should probably figure out what's going on
  #I also shouldn't be manually adding and subtracting 3 so frequently, but...
  k.max <- fmsfundata$k.max
  k.param.temp <- round(x[[1]] * (k.max + 3)) 
  if(k.param.temp < 3) {k.param.temp <- 3}  
  if(k.param.temp > (k.max + 3)) {k.param.temp <- k.max}
  res <- x[[2]] * 2.37 + 0.03
  if(res < 0.03) {res <- 0.03} 
  if(res > 2.4) {res <- 2.4}
  temp.object <- FindNeighbors(fmsfundata$object,
                          dims = 1:fmsfundata$nPC,
                          k.param = k.param.temp,
                          verbose = F) %>%
    FindClusters(resolution = res,
                 verbose = F)
  icvi <- getICVI.nm(temp.object, fmsfundata$nPC) 
  return(list(f=icvi,
              g=c(),
              c=c(),
              index = index,
              this = list(costfargument = fmsfundata)))
}

getICVI.nm <- function(object, nPC = NULL) {
  icvi <- suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:nPC],
                                            method = "euclidean"),
                                   cl = as.integer(Idents(object))))
  if(is.nan(icvi)) {icvi <- 0}
  #lmao because we're minimizing here
  return((1-icvi))
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


#   Subclustering for NM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

as.int.factor <- function(x) {as.integer(levels(x))[x]}

adjustIdents <- function(object, temp.object, npc) {
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
  } else if (original.SI <= new.SI) {
    object@active.ident <- as.factor(combined.new)
    print(paste(object@project.name, "ACCEPTED sub clustering"))
    print(paste("Original:", original.SI, "New:", new.SI))
  }
  return(object)
}

nm2d <- function(x0, fmsfundata = NULL) {
  nm <- neldermead()
  nm <- neldermead.set(nm, "numberofvariables", 2)
  nm <- neldermead.set(nm, "function", scoreCluster)
  nm <- neldermead.set(nm, "x0", x0)
  nm <- neldermead.set(nm, 'costfargument', fmsfundata)
  nm <- neldermead.set(nm, "verbose", F)
  nm <- neldermead.set(nm, "storehistory", T)
  nm <- neldermead.set(nm, "verbosetermination", T)
  nm <- neldermead.set(nm, 'method', 'box')
  nm <- neldermead.set(nm, "boundsmin", c(0,0))
  nm <- neldermead.set(nm, "boundsmax", c(1,1))
  nm <- neldermead.set(nm, "maxiter", 25)
  nm <- neldermead.search(nm)
  return(nm)
}

nmSubCluster <- function(object, nPC, sc.iter = 2) {
  
  #determine Silhouette scores for each cluster
  clusterScores <- getSilScores(object, nPC)
  
  #iteratively subcluster sc.iter clusters with lowest SilScores 
  for(j in 1:sc.iter){
    sub.sil.results <- vector(mode = "numeric", length = 0)
    sub.opt.store <- data.frame(matrix(0, nrow = 5, ncol = 2))
    names(sub.opt.store) <- c("k.param", "resolution")
    
    sc.object <- subset(object, idents = clusterScores[j, 1])
    k.max <- min(c((ncol(sc.object)/2) %>% round() - 3,
                   237)
                 )
    
    starts <- maximinLHS(5,2)
    fmsfundata <- structure(
      list(object = sc.object,
           nPC = nPC,
           k.max = k.max),
      class = 'optimbase.functionargs'
    )
    for (i in 1:5){
      nm <- nm2d(x0 = optimbase::transpose(starts[i, ]),
                 fmsfundata = fmsfundata)
      # It's 1-fopt because this function minimized the value of 1-SI
      sub.sil.results[i] <- 1 - neldermead.get(nm, "fopt")
      #print(1-neldermead.get(nm, "fopt"))
      sub.opt.store[i, ] <- neldermead.get(nm, "xopt") %>% optimbase::transpose()
    }
    x <- sub.opt.store[which.max(sub.sil.results),]
    k.param.temp <- round(x[[1]] * (k.max + 3)) 
    if(k.param.temp < 3) {k.param.temp <- 3}  
    if(k.param.temp > (k.max + 3))  {k.param.temp <- k.max}
    res <- x[2] * 2.37 + 0.03
    if(res < 0.03) {res <- 0.03} 
    if(res > 2.4) {res <- 2.4}
    
    sc.object <- FindNeighbors(sc.object,
                               dims = 1:nPC,
                               k.param = k.param.temp,
                               verbose = F) %>%
      FindClusters(resolution = res,
                   verbose = T)
    
    object <- adjustIdents(object = object,
                                 temp.object = sc.object,
                                 npc = nPC)
  }
  
  return(c(getICVI.nm(object, nPC = nPC),
           ARI(object@active.ident, object$annotation))
         )
}

nmCluster <- function(dataset){
  start.time <- Sys.time()
  nPC <- Predict_nPCs(dataset)
  
  sil.results <- vector(mode = "numeric", length = 0)
  ari.results <- vector(mode = "numeric", length = 0)
  eval.count <- vector(mode = "integer", length = 0)
  opt.store <- data.frame(matrix(0, nrow = 5, ncol = 2))
  names(opt.store) <- c("k.param", "resolution")
  
  set.seed(start.time)
  starts <- optimumLHS(5,2)
  
  fmsfundata <- structure(
    list(object = dataset,
         nPC = nPC,
         k.max = 237),
    class = 'optimbase.functionargs'
  )
  
  for(i in 1:5) {
    x0 <- optimbase::transpose(starts[i,])
    nm <- neldermead()
    nm <- neldermead.set(nm, "numberofvariables", 2)
    nm <- neldermead.set(nm, "function", scoreCluster)
    nm <- neldermead.set(nm, "x0", x0)
    nm <- neldermead.set(nm, 'costfargument', fmsfundata)
    nm <- neldermead.set(nm, "verbose", F)
    nm <- neldermead.set(nm, "storehistory", T)
    nm <- neldermead.set(nm, "verbosetermination", T)
    nm <- neldermead.set(nm, 'method', 'box')
    nm <- neldermead.set(nm, "boundsmin", c(0,0))
    nm <- neldermead.set(nm, "boundsmax", c(1,1))
    nm <- neldermead.set(nm, "maxiter", 25)
    nm <- neldermead.search(nm)
    # It's 1-fopt because this function minimized the value of 1-SI
    sil.results[i] <- 1 - neldermead.get(nm, "fopt")
    eval.count[i] <- neldermead.get(nm, "funevals")
    opt.store[i, ] <- neldermead.get(nm, "xopt") %>% optimbase::transpose()
  }
  
  xopt <- opt.store[which.max(sil.results) , ]
  k.param.temp <- round(xopt[[1]] * 237 + 3)
  if(k.param.temp < 3) {k.param.temp <- 3}  
  if(k.param.temp > 240) {k.param.temp <- 240}
  res <- xopt[[2]] * 2.37 + 0.03
  if(res < 0.03) {res <- 0.03} 
  
  dataset <- FindNeighbors(dataset,
                           dims = 1:nPC,
                           k.param = k.param.temp,
                           verbose = F) %>%
    FindClusters(resolution = res,
                 verbose = F)
  
  clustered.SI <- 1 - getICVI.nm(dataset, nPC = nPC)
  clustered.ARI <- ARI(dataset@active.ident, dataset$annotation)
  
  sc.results<- nmSubCluster(dataset, nPC = nPC, sc.iter = 2)
  sc.SI <- sc.results[1]
  sc.ARI <- sc.results[2]
  end.time <- Sys.time() - start.time
  
  return(c(clustered.SI, clustered.ARI,
           sc.SI, sc.ARI,
           end.time))
}

#test <- nmCluster(zeisel)

'
"
test.sc <- nmSubCluster(object, xopt = xopt, nPC = 12, sc.iter = 3)


results.tab <- matrix(0, nrow = 5, ncol = 5)
datasets <- list(Zeisel = zeisel, Panc1 = panc.1, Panc2 = panc.2, PR = pr, MR = mr)
iter <- 1
sc.results.tab <- matrix(0, nrow = 5, ncol = 2) %>% as.data.frame()
names(sc.results.tab) <- c("Silhouette", "ARI")
set.seed(Sys.time())

testSeed <- function() {
  optimumLHS(5,2) %>% print()
}




write.table(sc.results.tab, file = "clipboard", sep = "\t")
write.table(results.tab, file = "clipboard", sep = "\t")
write.table(rbind(sil.results, ari.results, eval.count), file = "clipboard", sep = "\t")"'