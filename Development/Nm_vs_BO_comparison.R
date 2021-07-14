#NM vs. BO
datasets <- list(Zeisel = zeisel, Panc1 = panc.1, Panc2 = panc.2, PR = pr, MR = mr)
k.param.range <- c(3, 300)
res.range <- c(0.03, 2.4)

nmAllRuns <- vector(mode = "list", length = 100)
boAllRuns <- vector(mode = "list", length = 100)

nm.onerun <- vector(mode = "list", length = 5)
bo.onerun <- vector(mode = "list", length = 5)

for(i in 1:100) {
  nm.onerun <- lapply(datasets, nmCluster)
  nmAllRuns[[i]] <- nm.onerun
  
  bo.onerun <- lapply(datasets,
                      bayesSubCluster,
                      k.param.range = k.param.range,
                      res.range = res.range,
                      n.starts = 25,
                      n.iterations = 25,
                      epsilon = 0.01)
  boAllRuns[[i]] <- bo.onerun
  
  print(paste("Highest Iteration:", i))
}

pullResults <- function(list, dataset) {
  dataset.results <- lapply(list, pluck, dataset)
  result.table <- do.call(rbind, dataset.results) %>% data.frame()
  names(result.table) <- c("clustered.SI", "clustered.ARI",
                           "sc.SI", "sc.ARI",
                           "end.time")
  return(result.table)
}

bo.results.frame <- bo.sd.frame <- data.frame(matrix(0, nrow = 5, ncol = 5), row.names = names(datasets))
names(bo.results.frame) <- names(bo.sd.frame) <- c("clustered.SI", "clustered.ARI",
                                                   "sc.SI", "sc.ARI",
                                                   "end.time")
for (dataset in names(datasets)){
  bo.results.frame[dataset, ] <- pullResults(boAllRuns[1:16], dataset) %>% colMeans()
  bo.sd.frame[dataset, ] <- pullResults(boAllRuns[1:16], dataset) %>% sapply(sd)
}

nm.results.frame <- nm.sd.frame <- data.frame(matrix(0, nrow = 5, ncol = 5), row.names = names(datasets))
names(nm.results.frame) <- names(nm.sd.frame) <- c("clustered.SI", "clustered.ARI",
                                                   "sc.SI", "sc.ARI",
                                                   "end.time")
for (dataset in names(datasets)){
  nm.results.frame[dataset, ] <- pullResults(nmAllRuns[1:16], dataset) %>% colMeans()
  nm.sd.frame[dataset, ] <- pullResults(nmAllRuns[1:16], dataset) %>% sapply(sd)
}

write.table(nm.results.frame, sep = '\t', file = 'clipboard')

BOvsNM.trial <- list(BO.Results = boAllRuns[1:16],
                     NM.Results = nmAllRuns[1:16])

saveRDS(BOvsNM.trial, file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/BOvsNM.rds")
BOvsNM.trial <- readRDS("C:/Users/asolivai/Desktop/R_Files/AutoClustR/BOvsNM.rds")
rm(bbb, bo.onerun, bo.results.frame, bo.sd.frame, bo.test, boAllRuns, nm.onerun, nm.results.frame, nm.sd.frame, nmAllRuns)




