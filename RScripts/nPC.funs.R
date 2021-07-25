nPC.nom <- function(x, max.pcs = min(100, length(x)), ignore.PC1 = TRUE){
  
  x <- as.data.frame(cbind(1:max.pcs, x[1:max.pcs]))
  colnames(x) <- c("PC", "Root")
  if(ignore.PC1){
    x <- x[-1, ]
  }
  
  npc.candidates <- data.frame()
  
  for(i in 1:(nrow(x)- 3)){
    
    npc.candidates[i, 1] <- x[i, 1] + 3
    max.score <- -Inf
    
    for(j in 1:i){
      sample.points <- x[(j + 1):(i + 3), ]
      SLR.fit <- lm(formula = Root ~ PC, data = sample.points)
      y.hat <- SLR.fit$coefficients[[2]] * x[j, 1] + SLR.fit$coefficients[[1]]
      score <- (x[j, 2] - y.hat) / sqrt(sum(residuals(SLR.fit)^2) / (nrow(sample.points) - 2))
      if(is.finite(score) & score > max.score){
        max.score <- score
        npc.candidates[i, 2] <- x[j, 1]
      }
    }

  }

  colnames(npc.candidates) <- c("pcs.use", "n.pcs")

  return(value = npc.candidates)
}

nPC.vote <- function(x, pcs.use = min(100, max(x[, 1])), ignore.single.counts = TRUE){

  x <- x[x[, 1] <= pcs.use, ]

  count.table <- rle(sort(x[, 2]))
  count.table <- data.frame(n.pc = count.table$values, count = count.table$lengths)
  if(ignore.single.counts & max(count.table[, 2]) > 1){
    count.table <- count.table[count.table[, 2] > 1, ]
  }
  count.table[, 3] <- count.table[, 2] / (max(x[, 1]) - count.table[, 1] - 2)
  colnames(count.table)[[3]] <- "freq"
  
  n.pcs <- count.table[count.table[, 3] == max(count.table[, 3]), 1]

  return(value = n.pcs[[1]])
}

Select.nPC <- function(x, file.path, max.pcs = min(100, length(x)), pcs.use = max.pcs, ignore.PC1 = TRUE, ignore.single.counts = TRUE, do.plot = TRUE){

  library(ggplot2)

  npc.candidates <- nPC.nom(x, max.pcs = max.pcs, ignore.PC1 = ignore.PC1)

  npc.tab <- data.frame()
  for(i in 1:nrow(npc.candidates)){
    npc.tab[i, 1] <- npc.candidates[i, 1]
    npc.tab[i, 2] <- nPC.vote(npc.candidates, pcs.use = npc.candidates[i, 1], ignore.single.counts = ignore.single.counts)
  }
  
  colnames(npc.tab) <- c("max.pcs", "n.pcs")

  if(do.plot){
    npc.plot <- ggplot(data = npc.tab, mapping = aes(x = npc.tab[, 1], y = npc.tab[, 2])) + geom_step(size = 0.75) + theme_classic() + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) + scale_y_log10(breaks = sort(unique(npc.tab[, 2])), limits = c(1, 100)) + xlab("\nMax.PC") + ylab("n.PC\n")
    ggsave(filename = "nPC.plot.pdf", plot = npc.plot, path = file.path, width = 8.5, height = 5.0)
  }
  
  n.pcs <- npc.tab[npc.tab[, 1] == pcs.use, 2]

  return(value = list(npc.tab, n.pcs))
}