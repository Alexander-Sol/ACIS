#Exploring paramater space
library(dplyr)
library(purrr)
library(stats)
library(Matrix)
library(tidyverse)
library(scales)


constructGrid <- function(){
  p.grid <- matrix(nrow = 24, ncol = 24, data = 0) %>% data.frame()
  # For actual names, exp. names are easier to parse
  'names(p.grid) <- 2^(9:32/4) %>% round()
  row.names(p.grid) <- 2^(9:32/4)/100'
  names(p.grid) <- 9:32/4
  row.names(p.grid) <- 9:32/4
  return(p.grid)
}

pluckICVI <- function(index.frame, k.param, ICVI = "Sil") {
  index.frame[which(index.frame$K.Param == k.param),] %>% 
    arrange(Resolution) %>%
    pluck(ICVI)
}

extractGrid <- function(index.frame) {
  res.grid <- constructGrid()
  for (i in 1:24){
    res.grid[i] <- pluckICVI(index.frame = index.frame, 
                             k.param = round(2^(as.numeric(names(res.grid)[i]))))
  }
  return(res.grid)
}

index.list <- readRDS("C:/Users/asolivai/Desktop/R_Files/AutoClustR/New_Predict_Index_List.rds")

grid.list <- lapply(index.list, extractGrid)
names(grid.list) <- names(index.list)

#Latin Squares function
#Breaks if n.starts exceeds 24
constructLS <- function(n.starts) {
  sq.size <- c(24, 24) %/% n.starts
  unit.grid.starts <- list(sample(1:n.starts),
                           sample(1:n.starts))
  i <- j <- vector(mode = "integer", length = n.starts)
  i <- sq.size[1] * (unit.grid.starts[[1]] - 1) + sample(1:sq.size[1], n.starts, replace = T)
  j <- sq.size[2] * (unit.grid.starts[[2]] - 1) + sample(1:sq.size[2], n.starts, replace = T)
  latin.square <- list(row = i, col = j)
}

#Function for getting nonzero adjacents
adj <- function(x) {
  y <- c(-1,1) + x
  return(ifelse(y > 0 & y < 25, y, NA))
}

# Hill climb function. Input: ICVI grid, n.starts. Output: # of solutions "calculated", solution quality (x/max(grid))
hill.climb <- function(n.starts, icvi.grid) {
  latin.square <- constructLS(n.starts)
  step.grid <- matrix(data = 0, nrow = 24, ncol = 24)
  max.index <- vector(mode = "numeric", length = n.starts)
  icvi.matrix <- as.matrix(icvi.grid)
  
  for (i in 1:n.starts){
    icvi.row <- latin.square$row[i]
    icvi.col <- latin.square$col[i]
    index <- icvi.matrix[icvi.row, icvi.col]
    #Count the steps
    step.grid[icvi.row, icvi.col] <- 1
    
    repeat {
      if(is.na(index)) {
        break
      }
      neighbors <- c(icvi.matrix[adj(icvi.row), icvi.col],
                     icvi.matrix[icvi.row, adj(icvi.col)])
      step.grid[adj(icvi.row), icvi.col] <- 1
      step.grid[icvi.row, adj(icvi.col)] <- 1
      if (all(is.na(neighbors)) || is.na(index)) {
        break
      }
      if (!(index >= max(neighbors, na.rm = T))) {
        if(which.max(neighbors) %in% c(1,2)){
          icvi.row <- icvi.row + which.max(neighbors)%%2*-2 + 1
        } else if (which.max(neighbors) %in% c(3,4)) {
          icvi.col <- icvi.col + which.max(neighbors)%%2*-2 + 1
        }
        index <- icvi.matrix[icvi.row, icvi.col]
      } else {
        break
      }
    }
    
    max.index[i] <- index
  }
  'max.pct <- max.index / max(icvi.matrix, na.rm = T)
  return(data.frame(total.starts = rep(n.starts, n.starts),
                    max.index = max.index,
                    max.pct = max.pct,
                    start.num = 1:n.starts,
                    step.count = sum(step.grid)))'
  return(max(max.index, na.rm = T) / max(icvi.grid, na.rm = T))
}

# Hill climb function. Input: ICVI grid, n.starts. Output: # of solutions "calculated", solution quality (x/max(grid))
# Same as previous, steepest ascent and random restart
la.hill.climb <- function(n.starts, icvi.grid, lag.length = 5) {
  latin.square <- constructLS(n.starts)
  step.grid <- matrix(data = 0, nrow = 24, ncol = 24)
  max.index <- vector(mode = "numeric", length = n.starts)
  icvi.matrix <- as.matrix(icvi.grid)
  index <- vector(mode = "numeric", length = lag.length)
  
  for (i in 1:n.starts){
    icvi.row <- latin.square$row[i]
    icvi.col <- latin.square$col[i]
    index[1] <- icvi.matrix[icvi.row, icvi.col]
    #Count the steps
    step.grid[icvi.row, icvi.col] <- 1
    
    repeat {
      if(all(is.na(index))) {
        break
      }
      neighbors <- c(icvi.matrix[adj(icvi.row), icvi.col],
                     icvi.matrix[icvi.row, adj(icvi.col)])
      step.grid[adj(icvi.row), icvi.col] <- 1
      step.grid[icvi.row, adj(icvi.col)] <- 1
      if (all(is.na(neighbors)) || all(is.na(index))) {
        break
      }
      if (!(min(index, na.rm = T) >= max(neighbors, na.rm = T))) {
        if(which.max(neighbors) %in% c(1,2)){
          icvi.row <- icvi.row + which.max(neighbors)%%2*-2 + 1
        } else if (which.max(neighbors) %in% c(3,4)) {
          icvi.col <- icvi.col + which.max(neighbors)%%2*-2 + 1
        }
        index[1] <- icvi.matrix[icvi.row, icvi.col]
        index[2:lag.length] <- index[1:(lag.length-1)]
      } else {
        break
      }
    }
    
    max.index[i] <- max(index, na.rm = T)
  }
  'max.pct <- max.index / max(icvi.matrix, na.rm = T)
  return(data.frame(total.starts = rep(n.starts, n.starts),
                    max.index = max.index,
                    max.pct = max.pct,
                    start.num = 1:n.starts,
                    step.count = sum(step.grid)))'
  return(max(max.index, na.rm = T) / max(icvi.matrix, na.rm = T))
}

#Random/greedy late acceptance hill climb. Doesn't do steepest ascent, just picks a direction at random
#The next iteration would probably be continue in direction if step or previous step was positive
#Currently broken because there are no checks to maintain index values within bounds of icvi.grid
rla.hill.climb <- function(n.starts, icvi.grid, lag.length = 5) {
  latin.square <- constructLS(n.starts)
  step.grid <- matrix(data = 0, nrow = 24, ncol = 24)
  max.index <- vector(mode = "numeric", length = n.starts)
  icvi.matrix <- as.matrix(icvi.grid)
  index <- vector(mode = "numeric", length = lag.length)
  
  for (i in 1:n.starts){
    icvi.row <- latin.square$row[i]
    icvi.col <- latin.square$col[i]
    index[1] <- icvi.matrix[icvi.row, icvi.col]
    #Count the steps
    step.grid[icvi.row, icvi.col] <- 1
    
    repeat {
      if(all(is.na(index))) {
        break
      }
      neighbor <- sample(1:4, 1)
      if(neighbor %in% c(1,2)){
        icvi.row <- icvi.row + which.max(neighbors)%%2*-2 + 1
      } else if (neighbor %in% c(3,4)) {
        icvi.col <- icvi.col + which.max(neighbors)%%2*-2 + 1
      }
      step.grid[icvi.row, icvi.col] <- 1
      neighbor.score <- icvi.matrix[icvi.row, icvi.col]
      if (is.na(neighbor.score) || all(is.na(index))) {
        break
      }
      if (!(min(index, na.rm = T) >= neighbor.score)) {
        index[1] <- icvi.matrix[icvi.row, icvi.col]
        index[2:lag.length] <- index[1:(lag.length-1)]
      } else {
        break
      }
    }
    
    max.index[i] <- max(index, na.rm = T)
  }
  max.pct <- max.index / max(icvi.matrix, na.rm = T)
  return(data.frame(total.starts = rep(n.starts, n.starts),
                    max.index = max.index,
                    max.pct = max.pct,
                    start.num = 1:n.starts,
                    step.count = sum(step.grid)))
}





#Summary of what i'm trying to do
#I think this can actually be a recursive function. We're always restricting the grid size by 6, picking the middle
#of each grid in the squre, finding the best subset, then repeating. The only question is how to generalize the size of the subset
#idk this should be doable
grid.search <- function(n.starts = 4, icvi.grid) {
  sq.size <- c(24, 24) %/% n.starts
  #This should only randomly sample if square size is even. if odd, should just always pick the middle
  row.starts <- sq.size[1]*1:4 - sample(3:4, 4, replace = T) 
  col.starts <- sq.size[2]*1:4 - sample(3:4, 4, replace = T) 
  mesh <- icvi.grid[row.starts, col.starts]
  
  # Trying to reduce the solutions space. If it was a grid of 4 by  4, want to go to 3 by 3 (24 -> 18)
  # Basically, sum up every 4 adjacent points, choose the greatest
  # obviously is should run from 1 to some variable. I think 1:(n.starts-1)
  mesh.sums <- vector(mode = "numeric", length = 2^2)
  for (i in 1:2) {
    for (j in 1:2) {
      mesh.sums[(i-1)*2+j] <- mesh[i:(i+2), j:(j+2)] %>% sum()
    }
  } 
  mesh.row <- ceiling(which.max(mesh.sums) / 2)
  mesh.col <- which.max(mesh.sums) %% 2
  if (mesh.col == 0) {mesh.col <- 2} #There has to be a better way to do this

  #restricted grid size
  res.rows <- 1:18 + (mesh.row-1)*sq.size[1]
  res.cols <- 1:18 + (mesh.col - 1)*sq.size[2]
  res.grid <- icvi.grid[res.rows, res.cols]
  
  #to make this more generalizable, we should probably to prime factorization to determine how to subdivide space 
  #in general though, should subdivide into grid w/square length != sq.size
  prime.factor.placeholder <- 6
  second.sq.size <- c(18, 18) %/% prime.factor.placeholder
  #should randomly sample if square size is even
  row.starts <- second.sq.size[1]*1:prime.factor.placeholder - 1
  col.starts <- second.sq.size[2]*1:prime.factor.placeholder - 1
  res.mesh <- res.grid[row.starts, col.starts]
  
  mesh.sums <- vector(mode = "numeric", length = 3^2)
  for (i in 1:3) {
    for (j in 1:3) {
      mesh.sums[(i-1)*3+j] <- res.mesh[i:(i+3), j:(j+3)] %>% sum()
    }
  } 
  mesh.row <- ceiling(which.max(mesh.sums) / 3)
  mesh.col <- which.max(mesh.sums) %% 3
  if (mesh.col == 0) {mesh.col <- 3} #There has to be a better way to do this
  
  res.rows <- 1:12 + (mesh.row-1)*second.sq.size[1]
  res.cols <- 1:12 + (mesh.col-1)*second.sq.size[2]
  res.grid <- res.grid[res.rows, res.cols]
  
  sq.size <- c(12,12) %/% 6
  row.starts <- sq.size[1]*1:6 - sample(0:1, 6, replace = T)
  col.starts <- sq.size[2]*1:6 - sample(0:1, 6, replace = T)
  #Will have to add a check so that the same solutions aren't recalculated multiple times
  res.mesh <- res.grid[row.starts, col.starts]
  mesh.sums <- vector(mode = "numeric", length = 4^2)
  for (i in 1:4) {
    for (j in 1:4) {
      mesh.sums[(i-1)*4+j] <- res.mesh[i:(i+1), j:(j+1)] %>% sum()
    }
  } 
  mesh.row <- ceiling(which.max(mesh.sums) / 4)
  mesh.col <- which.max(mesh.sums) %% 4
  if (mesh.col == 0) {mesh.col <- 4} #There has to be a better way to do this

  res.rows <- 1:6 + (mesh.row-1)*sq.size[1]
  res.cols <- 1:6 + (mesh.col-1)*sq.size[2]
  res.grid <- res.grid[res.rows, res.cols]
  
  if(max(icvi.grid, na.rm = T) == max(res.grid, na.rm = T)) {
    return(1)
  } else {return(max(res.grid, na.rm = T)/max(icvi.grid, na.rm = T))}
}

grid.search(grid.list$MR)

sapply(rep(4, 100), grid.search, icvi.grid = grid.list$PR)
grid.frame <- data.frame(Zeisel = sapply(rep(4, 100), grid.search, icvi.grid = grid.list$Zeisel),
                         Panc1 = sapply(rep(4, 100), grid.search, icvi.grid = grid.list$Panc1),
                         Panc2 = sapply(rep(4, 100), grid.search, icvi.grid = grid.list$Panc2),
                         PR = sapply(rep(4, 100), grid.search, icvi.grid = grid.list$PR),
                         MR = sapply(rep(4, 100), grid.search, icvi.grid = grid.list$MR))
boxplot(grid.frame)

grid.collapse <- grid.frame %>%
  select(Zeisel, Panc1, Panc2, PR, MR) %>%
  gather(key = "dataset", value = "Solution_Quality")
ggplot(grid.collapse, aes(x = dataset, y = Solution_Quality, fill = dataset)) + geom_violin(scale = "width") + 
  scale_fill_brewer(palette = "Blues") + theme_minimal()

ggplot(la.max.collapse, aes(x = dataset, y = Solution_Quality, fill = dataset)) + geom_violin(scale = "width") + ylim(c(0.85, 1)) +
  scale_fill_brewer(palette = "Blues") + theme_minimal()

ggplot(hc.max.collapse, aes(x = dataset, y = Solution_Quality, fill = dataset)) + geom_violin(scale = "width") + ylim(c(0.85, 1)) +
  scale_fill_brewer(palette = "Blues") + theme_minimal()

#I think for this, n.starts should be set to 24. At least initially
expected.steps <- function(n.starts, icvi.grid, threshold = 1) {
  latin.square <- constructLS(n.starts)
  max.index <- vector(mode = "numeric", length = n.starts)
  step.grid <- matrix(data = 0, nrow = 24, ncol = 24)
  icvi.matrix <- as.matrix(icvi.grid)
  
  for (i in 1:n.starts){
    icvi.row <- latin.square$row[i]
    icvi.col <- latin.square$col[i]
    index <- icvi.matrix[icvi.row, icvi.col]
    #Count the steps
    step.grid[icvi.row, icvi.col] <- 1
    
    repeat {
      if(is.na(index)) {
        break
      }
      if(index >= threshold*max(icvi.matrix, na.rm = T)) {
        return(sum(step.grid))
      }
      neighbors <- c(icvi.matrix[adj(icvi.row), icvi.col],
                     icvi.matrix[icvi.row, adj(icvi.col)])
      step.grid[adj(icvi.row), icvi.col] <- 1
      step.grid[icvi.row, adj(icvi.col)] <- 1
      if (all(is.na(neighbors)) || is.na(index)) {
        break
      }
      if (!(index >= max(neighbors, na.rm = T))) {
        if(which.max(neighbors) %in% c(1,2)){
          icvi.row <- icvi.row + which.max(neighbors)%%2*-2 + 1
        } else if (which.max(neighbors) %in% c(3,4)) {
          icvi.col <- icvi.col + which.max(neighbors)%%2*-2 + 1
        }
        index <- icvi.matrix[icvi.row, icvi.col]
      } else {
        break
      }
    }
  }
  return(NA)
}

la.expected.steps <- function(n.starts, icvi.grid, lag.length = 5, threshold = 1) {
  latin.square <- constructLS(n.starts)
  step.grid <- matrix(data = 0, nrow = 24, ncol = 24)
  max.index <- vector(mode = "numeric", length = n.starts)
  icvi.matrix <- as.matrix(icvi.grid)
  index <- vector(mode = "numeric", length = lag.length)
  
  for (i in 1:n.starts){
    icvi.row <- latin.square$row[i]
    icvi.col <- latin.square$col[i]
    index[1] <- icvi.matrix[icvi.row, icvi.col]
    #Count the steps
    step.grid[icvi.row, icvi.col] <- 1
    
    repeat {
      if(all(is.na(index))) {
        break
      } 
      if (max(index, na.rm = T) >= threshold*max(icvi.matrix, na.rm = T)) {
        return(sum(step.grid))
      }
      neighbors <- c(icvi.matrix[adj(icvi.row), icvi.col],
                     icvi.matrix[icvi.row, adj(icvi.col)])
      step.grid[adj(icvi.row), icvi.col] <- 1
      step.grid[icvi.row, adj(icvi.col)] <- 1
      if (all(is.na(neighbors)) || all(is.na(index))) {
        break
      }
      if (!(min(index, na.rm = T) >= max(neighbors, na.rm = T))) {
        if(which.max(neighbors) %in% c(1,2)){
          icvi.row <- icvi.row + which.max(neighbors)%%2*-2 + 1
        } else if (which.max(neighbors) %in% c(3,4)) {
          icvi.col <- icvi.col + which.max(neighbors)%%2*-2 + 1
        }
        index[1] <- icvi.matrix[icvi.row, icvi.col]
        index[2:lag.length] <- index[1:(lag.length-1)]
      } else {
        break
      }
    }
  }
  return(NA)
}

#Currently broken because there are no checks to maintain index values within bounds of icvi.grid
rla.expected.steps <- function(n.starts, icvi.grid, lag.length = 5, threshold = 1) {
  latin.square <- constructLS(n.starts)
  step.grid <- matrix(data = 0, nrow = 24, ncol = 24)
  icvi.matrix <- as.matrix(icvi.grid)
  index <- vector(mode = "numeric", length = lag.length)
  
  for (i in 1:n.starts){
    icvi.row <- latin.square$row[i]
    icvi.col <- latin.square$col[i]
    index[1] <- icvi.matrix[icvi.row, icvi.col]
    #Count the steps
    step.grid[icvi.row, icvi.col] <- 1
    
    repeat {
      if(all(is.na(index))) {
        break
      } 
      if (max(index, na.rm = T) >= threshold*max(icvi.matrix, na.rm = T)) {
        return(sum(step.grid))
      }
      neighbor <- sample(1:4, 1)
      if(neighbor %in% c(1,2)){
        icvi.row <- icvi.row + neighbor + 1
      } else if (neighbor %in% c(3,4)) {
        icvi.col <- icvi.col + neighbor + 1
      }
      step.grid[icvi.row, icvi.col] <- 1
      neighbor.score <- icvi.matrix[icvi.row, icvi.col]
      if (is.na(neighbor.score) || all(is.na(index))) {
        break
      }
      if (!(min(index, na.rm = T) >= neighbor.score)) {
        index[1] <- icvi.matrix[icvi.row, icvi.col]
        index[2:lag.length] <- index[1:(lag.length-1)]
      } else {
        break
      }
    }
  }
  return(NA)
}


#Cumulative Probability Distribution for simulated data
sim.cdf <- function(x, steps = steps.values) {
  sum(steps[!is.na(steps)] < x) / length(steps)
}

#Cumulative Probability Distribution for random draws
#Proof can be found @ https://arxiv.org/pdf/1404.1161v1.pdf
arxiv.cdf <- function(n, K){
  1 - (choose(576-n, K) / choose(576, K))
}

draws.until.success <- data.frame(A_1 = sapply(1:576, arxiv.cdf, K = 1),
                                  B_2 = sapply(1:576, arxiv.cdf, K = 2),
                                  C_10 = sapply(1:576, arxiv.cdf, K = 10),
                                  D_50 = sapply(1:576, arxiv.cdf, K = 50),
                                  Number_of_Solutions_Tested = 1:576)
draws.collapse <- draws.until.success %>%
  select(Number_of_Solutions_Tested, A_1, B_2, C_10, D_50) %>%
  gather(key = "Number of \nacceptable solutions", value = "Probability of drawing an acceptable solution", -Number_of_Solutions_Tested)
draws.plot <- ggplot(data = draws.collapse, aes(x = Number_of_Solutions_Tested, y = `Probability of drawing an acceptable solution`)) +
  geom_line(aes(color = `Number of \nacceptable solutions`, linetype = `Number of \nacceptable solutions`), size = 1.2) +
  xlim(c(1, 576)) +
  ylim(c(0,1)) +
  ggtitle("Cumulative Probability of Selecting an Acceptable Solution")
print(draws.plot)


 
#   Testing local search vs randomly checking solutions   %%%%%%%%%%%%%%%%%%%%%%
n.starts <- 24
for(dataset in names(grid.list)) {
  for(acceptance.threshold in c(0.95, 1)) {
    steps.values    <- lapply(rep(n.starts,1000),
                             expected.steps,
                             icvi.grid = grid.list[[dataset]],
                             threshold = acceptance.threshold) %>% unlist()
    
    la.steps.values <- lapply(rep(n.starts,1000),
                              la.expected.steps,
                              icvi.grid = grid.list[[dataset]],
                              threshold = acceptance.threshold,
                              lag.length = 7) %>% unlist()

    #Define parameters for hyper geometric distribution
    black.balls <- (grid.list[[dataset]] >= acceptance.threshold*max(grid.list[[dataset]], na.rm = T)) %>%
      sum(na.rm = T)
    #Summarize and plot results
    probability.frame <- data.frame(hill.climb.prob = sapply(1:576, sim.cdf, steps = steps.values),
                                    hyper.geom.prob = sapply(1:576, arxiv.cdf, K = black.balls),
                                    la.hill.climb.prob = sapply(1:576, sim.cdf, steps = la.steps.values),
                                    number.of.solutions.tested = 1:576)
    
    
    pf.collapse <- probability.frame %>%
      select(number.of.solutions.tested, hill.climb.prob, hyper.geom.prob, la.hill.climb.prob) %>%
      gather(key = "method", value = "P_select_best", -number.of.solutions.tested)
    prob.plot <- ggplot(data = pf.collapse, aes(x = number.of.solutions.tested, y = P_select_best)) +
      geom_line(aes(color = method, linetype = method), size = 1.2) +
      xlim(c(1, 576)) +
      ylim(c(0,1)) +
      ggtitle(paste0("Cumulative Probability of selecting solution ",
                     acceptance.threshold*100, "% of best possible solution: ",
                     dataset))
    print(prob.plot)
  }
}


# Find distribution of result quality
n.starts <- 24
la.steps.values <- data.frame(matrix(nrow=100,ncol=5))
names(la.steps.values) <- names(grid.list)
for(dataset in names(grid.list)) {
  la.steps.values[[dataset]] <- lapply(rep(n.starts,100),
                                  la.hill.climb,
                                  icvi.grid = grid.list[[dataset]],
                                  lag.length = 7) %>% unlist()
}

la.max.collapse <- la.steps.values %>%
  select(Zeisel, Panc1, Panc2, PR, MR) %>%
  gather(key = "dataset", value = "Solution_Quality")
ggplot(la.max.collapse, aes(x = dataset, y = Solution_Quality)) + geom_violin(scale = "width") + ylim(c(0.85,1))

n.starts <- 24
hc.steps.values <- data.frame(matrix(nrow=100,ncol=5))
names(hc.steps.values) <- names(grid.list)
for(dataset in names(grid.list)) {
  hc.steps.values[[dataset]] <- lapply(rep(n.starts,100),
                                       hill.climb,
                                       icvi.grid = grid.list[[dataset]]) %>% 
    unlist()
}

hc.max.collapse <- hc.steps.values %>%
  select(Zeisel, Panc1, Panc2, PR, MR) %>%
  gather(key = "dataset", value = "Solution_Quality")
ggplot(hc.max.collapse, aes(x = dataset, y = Solution_Quality)) + geom_violin(scale = "width") + ylim(c(0.85,1))


#   Testing different numbers of starts for random search   %%%%%%%%%%%%%%%%%%%%

for(dataset in names(grid.list)) {
  for(acceptance.threshold in c(0.95, 1)) {
    steps.by.start <- vector(mode = "list", length = "4")
    for (n.starts in c(6, 12, 18, 24)) {
      steps.by.start[[n.starts/6]] <- lapply(rep(n.starts, 500),
                                             la.expected.steps,
                                             icvi.grid = grid.list[[dataset]],
                                             threshold = acceptance.threshold,
                                             lag.length = 7) %>% unlist()
    }
    #Define parameters for hyper geometric distribution
    black.balls <- (grid.list[[dataset]] >= acceptance.threshold*max(grid.list[[dataset]], na.rm = T)) %>%
      sum(na.rm = T)
    #Summarize and plot results
    probability.frame <- data.frame(la.6.starts = sapply(1:576, sim.cdf, steps = steps.by.start[[1]]),
                                    la.12.starts = sapply(1:576, sim.cdf, steps = steps.by.start[[2]]),
                                    la.18.starts = sapply(1:576, sim.cdf, steps = steps.by.start[[3]]),
                                    la.24.starts = sapply(1:576, sim.cdf, steps = steps.by.start[[4]]),
                                    hyper.geom.prob = sapply(1:576, arxiv.cdf, K = black.balls),
                                    number.of.solutions.tested = 1:576)
    
    
    pf.collapse <- probability.frame %>%
      select(number.of.solutions.tested, la.6.starts, la.12.starts, la.18.starts, la.24.starts, hyper.geom.prob) %>%
      gather(key = "method", value = "P_select_best", -number.of.solutions.tested)
    prob.plot <- ggplot(data = pf.collapse, aes(x = number.of.solutions.tested, y = P_select_best)) +
      geom_line(aes(color = method, linetype = method), size = 1.2) +
      xlim(c(1, 576)) +
      ylim(c(0,1)) +
      ggtitle(paste0("Cumulative Probability of selecting solution ",
                     acceptance.threshold*100, "% of best possible solution: ",
                     dataset))
    print(prob.plot)
  }
}




#  Full Summary of Hill Climbing algorithm   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

hc.results <- lapply(rep(1:24, 10), hill.climb, icvi.grid = grid.list$PR)
hc.summary.table <- bind_rows(hc.results)
#Quickly visualize results
boxplot(max.pct ~ total.starts, data = hc.summary.table, ylim = c(0.5,1))
boxplot(step.count ~ total.starts, data = hc.summary.table)



la.hc.results <- lapply(1:24, la.hill.climb, icvi.grid = grid.list$Zeisel)
la.hc.summary.table <- bind_rows(la.hc.results)
#Quickly visualize results
boxplot(max.pct ~ total.starts, data = la.hc.summary.table, ylim = c(0.5,1))
plot(step.count ~ total.starts, data = la.hc.summary.table)



#   False max stuff, all tangentially related   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constructTable <- function(index.frame, ICVI = "Sil"){
  mr.table <- data.frame(ICVI = index.frame[ICVI], K.param = index.frame$K.Param, Resolution = index.frame$Resolution)
  mr.table$logRes <-  log2(mr.table$Resolution*100)
  
  k.map <- 2^(9:32/4) %>% round()
  names(k.map) <- 9:32/4
  for(i in 1:nrow(mr.table)) {
    mr.table$logK[i] <- names(k.map)[k.map == mr.table$K.param[i]]
  }
  mr.table$best <- ifelse(sapply(1:nrow(mr.table),
                                 FUN = isfalseMax,
                                 res.frame = mr.table,
                                 icvi = ICVI), "O", "")
  mr.table$best <- ifelse(index.frame[ICVI] == max(index.frame[ICVI], na.rm = T), "X", mr.table$best)
  return(mr.table)
}

#This is cursed, has absolutely no error handling, breaks if you don't pass the correct icvi/table combo, ect.
isfalseMax <- function(row, res.frame, icvi){
  k.val <- res.frame$logK[row]
  r.val <- res.frame$logRes[row]
  if(k.val %in% c(2.25, 8) | r.val %in% c(2.25, 8)) {
    return(F)
  } else {
    p.key <- 9:32/4
    
    ind.val <- res.frame[res.frame$logK == k.val & res.frame$logRes == r.val, icvi]
    if(is.na(ind.val)){
      return(F)
    }
    
    #Retrieve scores for k.param neighbors
    k.n <- p.key[c(which(p.key == k.val)-1, which(p.key == k.val)+1)]
    k.n.val <- res.frame[res.frame$logK %in% k.n & res.frame$logRes == r.val, icvi]
    
    #Retrieve scores for resolution neighbors
    r.n <- p.key[c(which(p.key == r.val)-1, which(p.key == r.val)+1)]
    r.n.val <- res.frame[res.frame$logRes %in% r.n & res.frame$logK == k.val, icvi]
    
    if (ind.val >= max(c(k.n.val, r.n.val),
                       na.rm = T)
    ) {
      return(T)
    } else {return(F)}
  }
}


# Make heatmap w/ local maxima for CH Index   %%%%%%%%%%%%%
grid.res.ch <- constructTable(results$MR, ICVI = "CH")

ggplot(grid.res.ch, aes_string(x = "logK", y = "logRes", fill = "CH")) +
  geom_tile() +
  scale_fill_gradient(low = "gold1", high = "darkred") +
  geom_text(aes_string(label = "best"), size = 6) +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18)) +
  labs(x = expression("\nLog"[2]*"K.param"),
       y = expression("Log"[2]*"(Resolution*100)"), fill = "CH") + 
  ggtitle(label = "MR Dataset: CH Index")

table(grid.res.ch$best)


# Make heatmap w/ local maxima for Silhouette Index   %%%%%%%%%%%%%
grid.res.sil <- constructTable(results$MR, ICVI = "Sil")

ggplot(grid.res.sil, aes_string(x = "logK", y = "logRes", fill = "Sil")) +
  geom_tile() +
  scale_fill_gradient(low = "gold1", high = "darkred") +
  geom_text(aes_string(label = "best"), size = 6) +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18)) +
  labs(x = expression("\nLog"[2]*"K.param"),
       y = expression("Log"[2]*"(Resolution*100)"), fill = "Sil") + 
  ggtitle(label = "MR Dataset: Silhouette Index")

table(grid.res.sil$best)


View(results$MR)






