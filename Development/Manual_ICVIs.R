#Manual SI calculations for Panc.2

# Load PCA results
LoadEmbeddings <- function(object, standardize){
  npc <- object@tools$Choose_nPCs
  all.embeddings <- Embeddings(object)[, 1:npc]
  std_dev <- object[["pca"]]@stdev[1:npc]
  if(standardize == TRUE){
    all.embeddings <- t(t(apply(all.embeddings, MARGIN = 2,
                                FUN = function(x) {x - mean(x)})) # subtracts PC mean from PC position
                        /std_dev) # Divide by standard deviation
  }
  Tool(object) <- all.embeddings
  return(object)
}

#Panc.2 autoclustrd for SI
auto.obj <- day80
DimPlot(auto.obj, label = T)

#Panc.2 with original idents
object
DimPlot(object, label = T)

#positions of each cell w/in PCA space
#This is a pretty hacky way to do this, revist later
nonstd.embed <- LoadEmbeddings(auto.obj, standardize = F)
nonstd.embed <- nonstd.embed@tools$LoadEmbeddings

#Cluster Identiesnof each cell
clusters <- Idents(auto.obj)

clust.space <- cbind(nonstd.embed, clusters)

dist.test <- dist(nonstd.embed[1:5,], method = "euclidean", upper = F)
vec.test <- as.numeric(dist.test)
dist.test
sd(dist.test)
head(clust.space)

#pull out every cell in cluster 2
clust.2 <- clust.space[clust.space[ ,ncol(clust.space)] == 2, ]
head(clust.2)

#Find the distance between every cell in cluster 2
library(proxy)
dist.2 <- dist(clust.2[, 1:(ncol(clust.space)-1)], upper = T)
print(dist.2)

#Return the distance between cells as vector
distance.vec <- colSums.dist(dist.2)

#View the distance distribution
#Have to figure out some way to pick out outliers
plot(density(log(distance.vec)))
sd(distance.vec)


#Automate this for all clusters

#Lapply returns a list of matrices of length = clust.num
test <- lapply(split(clust.space[, 1:(ncol(clust.space-1))], clust.space[,ncol(clust.space)]), matrix, ncol = ncol(clust.space-1))

#Get the distance vector for each
getDistance <- function(mat) {
  distance.vec <- dist(mat) %>%
    colSums.dist()/((nrow(mat)*(nrow(mat)-1))/2)
}

getDistSD <- function(mat) {
  dist(mat) %>% 
    sd()
}

getSquaredDistSD <- function(mat) {
  dist(mat, method = "Euclidean", upper = F)^2 %>% 
    sd()
}

test.vec <- lapply(test, getDistance)

lapply(test, getSquaredDistSD)
lapply(test.vec, mean)
z.score.test <- lapply(test.vec, zScore) 
plot(density(z.score.test$`9`))
lapply(z.score.test, range)

zScore <- function(vec) {
  (vec - mean(vec))/sd(vec)
}

#Cluster 6 has the highest SD, and is also three clusters merged into one. 
#What is a principled way of 

cluster.5 <- subset(auto.obj, idents = 5)
DimPlot(cluster.5)
cluster.5 <- FindNeighbors(cluster.5,
                            k.param = 7,
                            compute.SNN = TRUE,
                            prune.SNN = 1/15,
                            nn.method = "annoy",
                            nn.eps = 0.0,
                            verbose = FALSE,
                            reduction = "pca",
                            dims = 1:object@tools$Choose_nPCs) %>%
  
  # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
  FindClusters(modularity.fxn = 1,
               resolution = 0.1,
               algorithm = 1,
               group.singletons = TRUE,
               verbose = FALSE)

DimPlot(cluster.5, pt.size = 3)
DimPlot(object, label = T)


#   clusterCrit   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clusters <- as.integer(Idents(day80.hires))
ind.hires <- clusterCrit::intCriteria(nonstd.embed,
                                part = clusters,
                                crit = c("cal", "dav", "sil"))


clusters <- as.integer(Idents(day80))
ind.lowres <- clusterCrit::intCriteria(nonstd.embed,
                                      part = clusters,
                                      crit = c("cal", "dav", "sil"))


ind.lowres

#So for this, only the CH index rates the hires solution better (I believe hires is the "correct solution")


library(proxy)
auto.obj <- LoadEmbeddings(day80, standardize = F)
embeddings <- auto.obj@tools$LoadEmbeddings
clusters <- Idents(auto.obj)
clust.space <- cbind(embeddings, clusters)

#Return list of matrices containing embeddings for each cluster
embeddings.by.cluster <- lapply(split(clust.space[, 1:(ncol(clust.space-1))],
                                       clust.space[,ncol(clust.space)]),
                                 matrix, ncol = ncol(clust.space-1))

intraClusterVariance <- lapply(embeddings.by.cluster, getDistSD)
intraClusterSquared <- lapply(embeddings.by.cluster, getSquaredSD)
variance.plots <- lapply(embeddings.by.cluster, getDistplot)
intraClusterVariance

getDistSD <- function(mat) {
  dist(mat) %>% 
    sd()
}

getSquaredSD <- function(mat) {
  dist(mat)^2 %>% 
    sd()
}

getDistplot <- function(mat) {
  dist(mat) %>% 
    density() %>%
    plot(main = paste("Cluster ", (mat[1, ncol(mat)] - 1)))
}

install.packages("cluster")
library(cluster)

sil.test <- cluster::silhouette(as.integer(clust.space[ ,ncol(clust.space)]), dist(clust.space[, 1:(ncol(clust.space-1))]))
plot(density(sil.test[, 3]), main = "Silhouette Scores Day 80") + abline(v = 0, col = "red", lwd = 3, lty = 2)
sil.processed <- cluster::sortSilhouette(sil.test)
sil <- as.data.frame(cbind(cluster = sil.processed[,1], neighbor = sil.processed[,2], si_width = sil.processed[,3]))

new.sil <- split(sil, sil$cluster)

sil.plots <- lapply(new.sil, plotDensity)

plotDensity <- function(sil) {
  density(sil[ , 3]) %>%
    plot(main = paste("Cluster ", sil[1,1] - 1)) +
    abline(v = 0, col = "red", lwd = 3, lty = 2)
}

mean(sil[ , 3])

select.cells <- CellSelector(DimPlot(auto.obj))

auto.obj <- day80

auto.obj$clusters <- as.integer(auto.obj$seurat_clusters)

auto.obj$clusters[select.cells] <- 12

clusters <- auto.obj$clusters
table(clusters)

clusters <- auto.obj@active.ident

clust.space <- cbind(embeddings, clusters)

sil.2 <- cluster::silhouette(as.integer(clust.space[ ,ncol(clust.space)]), dist(clust.space[, 1:(ncol(clust.space-1))]))
sil.2.processed <- cluster::sortSilhouette(sil.2)
sil.2 <- as.data.frame(cbind(cluster = sil.2.processed[,1], neighbor = sil.2.processed[,2], si_width = sil.2.processed[,3]))
sil.2.split <- split(sil.2, sil.2$cluster)
mean(sil.2[ , 3])

lapply(sil.2.split, plotDensity)
lapply(new.sil, plotDensity)

auto.obj$sil <- sil.2[ , 3]
FeaturePlot(auto.obj, features = "sil", max.cutoff = 0.1, min.cutoff = -0.1, cols = feature.cols, pt.size = 2)


cluster.8 <- subset(day80, ident = 8)
length(Cells(cluster.8))
cluster.8 <- FindNeighbors(cluster.8, dims = 1:10, k.param = 20)
cluster.8 <- FindClusters(cluster.8, resolution = 0.3)

DimPlot(cluster.8)



# P1S again %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

object <- FindNeighbors(object,
                     k.param = 7,
                     compute.SNN = TRUE,
                     prune.SNN = 1/15,
                     nn.method = "annoy",
                     nn.eps = 0.0,
                     verbose = FALSE,
                     reduction = "pca",
                     dims = 1:object@tools$Choose_nPCs) %>%
  
  # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
  FindClusters(modularity.fxn = 1,
               resolution = 0.04756828,
               algorithm = 1,
               group.singletons = TRUE,
               verbose = FALSE)

DimPlot(object)



clust.space <- cbind(nonstd.embed, object$seurat_clusters)


sil.initial <- cluster::silhouette(as.integer(clust.space[ ,ncol(clust.space)]), dist(clust.space[, 1:(ncol(clust.space-1))]))
sil <- cluster::sortSilhouette(sil.initial)
sil <- as.data.frame(cbind(cluster = sil[,1], neighbor = sil[,2], si_width = sil[,3]))
sil.2.split <- split(sil, sil$cluster)
mean(sil[ , 3])
plot(sil.initial)

plot(density(sil.2.split[[2]][[,3]]))


scores.by.cluster <- aggregate(sil$si_width, list(sil$cluster), mean) %>%
  arrange(x)

subcluster.candidates <- scores.by.cluster[1:2, 1] - 1

intCriteria(nonstd.embed, as.integer(p1s$seurat_clusters), "sil")

lapply(sil.2.split, plotDensity)


DimPlot(p1s, pt.size = 2, label = T)

cluster.3 <- subset(p1s, idents = 3)

Choose_nPCs(cluster.3, pc.select = "segment")

cluster.3 <-    FindNeighbors(cluster.3,
                  k.param = 7,
                  compute.SNN = TRUE,
                  prune.SNN = 1/15,
                  nn.method = "annoy",
                  nn.eps = 0.0,
                  verbose = FALSE,
                  reduction = "pca",
                  dims = 1:8) %>%
      
      # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
      FindClusters(modularity.fxn = 1,
                   resolution = 0.1,
                   algorithm = 1,
                   group.singletons = TRUE,
                   verbose = FALSE)

DimPlot(cluster.3)



DimPlot(p1s, group.by = "annotation", pt.size = 2)

class(cluster.3$seurat_clusters)

new.idents <- as.numeric(cluster.3$seurat_clusters) + 7
names(new.idents) <- Cells(cluster.3)

combined.new <- p1s$seurat_clusters %>%
  as.integer() 
names(combined.new) <- Cells(p1s)

head(combined.new)

class(combined.new)

combined.new[names(new.idents)] <- new.idents

intCriteria(nonstd.embed, as.integer(p1s$seurat_clusters), c("cal", "sil", "dav"))
intCriteria(nonstd.embed, as.integer(combined.new), c("cal", "sil", "dav"))

intCriteria(std.embed, as.integer(p1s$seurat_clusters), c("cal", "sil", "dav"))
intCriteria(std.embed, as.integer(combined.new), c("cal", "sil", "dav"))

p1s <- Choose_nPCs(p1s, pc.use = 8, pc.select = "manual")
object <- p1s

object <- LoadEmbeddings(object, standardize = T)
std.embed <- object@tools$LoadEmbeddings
object <- LoadEmbeddings(object, standardize = F)
nonstd.embed <- object@tools$LoadEmbeddings










#   TSCAN segmented regression test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sdev <- object@reductions$pca@stdev[1:50]
x <- 1:50
optpoint <- which.min(sapply(2:40, function(i) {
  x2 <- pmax(0, x - i)
  sum(lm(sdev ~ x + x2)$residuals^2)
}))
pcadim = optpoint + 1


cluster.1 <- subset(p1s, idents = 1)
DimPlot(cluster.1)



cluster.1 <-    FindNeighbors(cluster.1,
                              k.param = 7,
                              compute.SNN = TRUE,
                              prune.SNN = 1/15,
                              nn.method = "annoy",
                              nn.eps = 0.0,
                              verbose = FALSE,
                              reduction = "pca",
                              dims = 1:8) %>%
  
  # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
  FindClusters(modularity.fxn = 1,
               resolution = 0.1,
               algorithm = 1,
               group.singletons = TRUE,
               verbose = FALSE)

DimPlot(cluster.1)




#   AutoClustR RunTimes for Recursive Clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cluster.1 <- subset(object, idents = 1)
length(Cells(cluster.1))

cluster.1.auto <- AutoClustR(cluster.1, pc.use = 8, pc.select = "manual", merge.clusters = F, standardize = F)

