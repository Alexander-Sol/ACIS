library(clusterSim)

unit.starts <- maximinLHS(16, 2) %>% data.frame() %>% setNames(c("x", "y"))
plot(unit.starts$x, unit.starts$y)

scaleParams <- function(unit.value, param.range) {
  unit.value * (param.range[2] - param.range[1]) + param.range[1]
}

#I think the minimum k.param is 2. Definitely throws an error if you pass 1 to Find Neighbors
Run.NM <- function(object, nPC = 10,
                   vars.to.tune = c("ndims", "k.param", "resolution"),
                   param.space = list(pc.space = c(2, 25),
                                      k.space = c(3, 240),
                                      res.space = c(0.0, 2.4)),
                   alpha = 1.0,
                   gamma = 2.0,
                   beta = 0.5,
                   delta = 0.5, 
                   tol = c(0, 0, 0.02), 
                   max.iter = 10) {


  
  dims <- length(vars.to.tune)
  vertices <- matrix(runif((dims + 1) * dims), nrow = dims + 1, ncol = dims)
  vertices <- cbind(vertices, rep(NA, dims + 1))
  iter <- 0
  
  # Calculate Silhouette index for initial vertices
  for(j in 1:(dims + 1)){
    Seurat.object <- FindNeighbors(Seurat.object,
                                   dims = 1:round(scaleParams(vertices[j, 1], param.space$pc.space)),
                                   k.param = scaleParams(vertices[j, 2], param.space$k.space) %>% round(),
                                   verbose = F) %>%
      FindClusters(resolution = scaleParams(vertices[j, 3], param.space$res.space),
                   verbose = F)
    
    vertices[j, ] <- appendICVI(vertices[j, 1:3], Seurat.object, nPC = nPC) 
  }

  # Run loop until max iterations reached
  while(iter <= max.iter){
    
    
    print(paste0("Iteration: ", iter))
    iter <- iter + 1

    # Determine best, worst, and second worst vertices
    # Holy shit it should only choose one vertex those. Picking multiple will obviously fuck everything up. 
    # Dan this is so broken I can't believe you said this was working
    best.vertex <- which.max(vertices[ , 4])[1]
    worst.vertex <- which.min(vertices[ , 4])[1]
    secworst.vertex <- which(vertices[, dims + 1] == min(vertices[-worst.vertex, dims + 1]))[1]
    print(paste0("Best SI: ", vertices[best.vertex, dims+1]))

    # Calculate centroid and axis for reflection, expansion, contraction
    centroid <- colMeans(vertices[-worst.vertex, -(dims + 1)])
    axis <- vertices[worst.vertex, -(dims + 1)] - centroid

    # Try reflection
    ref.vertex <- centroid - (axis * alpha)
    for(i in 1:length(ref.vertex)){
      if(ref.vertex[[i]] > 1){ref.vertex[[i]] <- 1
      }else if(ref.vertex[[i]] < 0){
        ref.vertex[[i]] <- 0
      }
    }
    
    Seurat.object <- FindNeighbors(Seurat.object,
                                   dims = 1:round(scaleParams(ref.vertex[[1]], param.space$pc.space)), 
                                   k.param = scaleParams(ref.vertex[[2]], param.space$k.space) %>% round(),
                                   verbose = F) %>%
      FindClusters(resolution = scaleParams(ref.vertex[[3]], param.space$res.space), verbose = F)
    
    ref.vertex <- appendICVI(ref.vertex, Seurat.object, nPC = nPC)
  
    # Try expansion if reflection is best
    if(ref.vertex[[dims + 1]] < vertices[best.vertex, dims + 1]){
      exp.vertex <- centroid - (axis * gamma)
      
      for(i in 1:length(exp.vertex)){
        if(exp.vertex[[i]] > 1){exp.vertex[[i]] <- 1
        }else if(exp.vertex[[i]] < 0){
          exp.vertex[[i]] <- 0
        }
      }
      
      Seurat.object <- FindNeighbors(Seurat.object,
                                     dims = 1:as.integer(exp.vertex[[1]] * (param.space[[1]][[2]] - param.space[[1]][[1]]) + param.space[[1]][[1]]),
                                     k.param = as.integer(exp.vertex[[2]] * (param.space[[2]][[2]] - param.space[[2]][[1]]) + param.space[[2]][[1]])) %>%
        FindClusters(resolution = exp.vertex[[3]] * (param.space[[3]][[2]] - param.space[[3]][[1]]) + param.space[[3]][[1]])
      exp.vertex <- appendICVI(exp.vertex, Seurat.object, nPC = nPC)
      
      if(exp.vertex[[dims + 1]] < ref.vertex[[dims + 1]]){
        vertices[worst.vertex, ] <- exp.vertex
      }else{
        vertices[worst.vertex, ] <- ref.vertex
      }
    
    # Use reflection if better than second worst point
    } else if(ref.vertex[[dims + 1]] < vertices[secworst.vertex, dims + 1]){
      vertices[worst.vertex, ] <- ref.vertex
    
    # Otherwise, try contraction
    } else{
      con.vertex <- centroid + (axis * beta)

      Seurat.object <- FindNeighbors(Seurat.object, dims = 1:as.integer(con.vertex[[1]] * (param.space[[1]][[2]] - param.space[[1]][[1]]) + param.space[[1]][[1]]), k.param = as.integer(con.vertex[[2]] * (param.space[[2]][[2]] - param.space[[2]][[1]]) + param.space[[2]][[1]])) %>%
        FindClusters(resolution = con.vertex[[3]] * (param.space[[3]][[2]] - param.space[[3]][[1]]) + param.space[[3]][[1]])
      con.vertex <- appendICVI(con.vertex, Seurat.object, nPC = nPC)
    
      if(con.vertex[[dims + 1]] < vertices[worst.vertex, dims + 1]){
        vertices[worst.vertex, ] <- con.vertex
      
      # As last resort, shrink simplex
      }else{
        shrink.vertices <- vertices[-best.vertex, -(dims + 1)] %>% as.matrix()
        shrink.axes <- matrix(0, nrow = dims, ncol = dims)
        shrink.scores <- vector(mode = "numeric", length = dims)
        for(m in 1:dims){
          shrink.axes[m, ] <- shrink.vertices[m, ] - vertices[best.vertex, -(dims + 1)]
          shrink.vertices[m, ] <- shrink.vertices[m, ] - (shrink.axes[m, ] * delta)
          Seurat.object <- FindNeighbors(Seurat.object, dims = 1:as.integer(shrink.vertices[m, 1] * (param.space[[1]][[2]] - param.space[[1]][[1]]) + param.space[[1]][[1]]), k.param = as.integer(shrink.vertices[m, 2] * (param.space[[2]][[2]] - param.space[[2]][[1]]) + param.space[[2]][[1]])) %>%
            FindClusters(resolution = shrink.vertices[m, 3] * (param.space[[3]][[2]] - param.space[[3]][[1]]) + param.space[[3]][[1]])
          shrink.scores[m] <- index.S(d = dist(Embeddings(Seurat.object, reduction = "pca")[, 1:nPC], method = "euclidean"),
                                      cl = as.integer(Idents(Seurat.object)))
        }
        shrink.vertices <- cbind(shrink.vertices, shrink.scores)
        colnames(shrink.vertices) <- NULL
        vertices <- rbind(shrink.vertices, vertices[best.vertex, ])
      }
    }
  }

  return(value = list(best.vertex))
}

test <- Run.NM(panc.1)

appendICVI <- function(vertex, object, nPC) {
  icvi <- suppressWarnings(index.S(d = dist(Embeddings(object, reduction = "pca")[, 1:nPC],
                                    method = "euclidean"),
                                    cl = as.integer(Idents(object))))
  if(is.nan(icvi)) {icvi <- 0}
  return(append(vertex, icvi))
}






