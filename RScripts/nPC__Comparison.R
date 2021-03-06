# Compiling nPC functions

# AutoClustR (possibly outdated)
# Only works for Seurat objects, will have to update
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


SD.table <- data.frame(PC = 1:pc.use,
                       SD = b1[["pca"]]@stdev[1:pc.use])


# CIDR
#' @title nPC Calculation
#'
#' @rdname calc_npc
#' @name calc_npc
#'
#' @description calculates nPC.
#'
#' @details
#' This function is used internally in the method \code{nPC};
#' see the help page for \code{nPC} for more details.
#'
#'
#' @param var a vector; proportion of variation explainded by each of the principal coordinates.
#'
#' @return This function is used internally in the method \code{nPC};
#' see the help page for \code{nPC} for more details.
#'
#'
#' @export calc_npc
#' 
#'

# var isn't a fixed value, but I'm pretty sure it's huge (like 1000+)
# This is only really applicable if you run the rest of CIDR's preprocessing workflow

calc_npc <- function(var, N=1, cutoff_divisor=10) {
  ## only use the first 90% of var - ignore outliers in last 10%
  l <- length(var)
  if (l > 0) {
    l <- as.integer(0.9 * l) 
    if (l == 0) {
      l <- len(var)
    }
    var <- var[1:l]
  }
  NPC_DEFAULT <- 4
  d <- var[-length(var)] - var[-1]
  descending_d <- sort(d, decreasing=T)
  max_d <- descending_d[1]
  mean_N_max <- mean(descending_d[1:N])
  ## get measure of spread of data as a proportion of the largest 3 differences
  spread <- mean(d)/mean(descending_d[1:3])*100
  if (spread > 15) {
    ## data too evenly spread for this algorithm to converge?
    return(NPC_DEFAULT)
  } else if (spread > 10) {
    cutoff_divisor <- 5
  }
  cutoff <- mean_N_max/cutoff_divisor
  groups <- list()
  groups[[1]] <- c(1)
  group_index <- 1
  for (i in 2:length(var)) {
    if (d[i-1] < cutoff) {
      ## include in current group
      groups[[group_index]] <- c(groups[[group_index]], i)
      ## check stopping criteria
      len_curr_group <- length(groups[[group_index]])
      if (len_curr_group > 7 || (i > 10 && len_curr_group > 3)) {
        ## return the last index in the previous group
        ## catch case where only 1 group
        if (group_index == 1) {
          return(NPC_DEFAULT)
        }
        prev_group <- groups[[group_index - 1]]
        nPC <- prev_group[length(prev_group)]
        ## not helpful to return nPC=1, so use default in this case
        if (nPC == 1) {
          nPC <- NPC_DEFAULT
        }
        return(nPC)
      }
    } else {
      ## start new group
      group_index <- group_index + 1
      groups[[group_index]] <- i
    } 
  }
  ## must not have converged - return default
  return(NPC_DEFAULT) 
}

# Seurat Jackstraw Example
b1 <- Proc.data(b1, algorithm = "seurat", expr.meas = "umi") 
b1 <- JackStraw(b1, num.replicate = 100, dims = 35) %>%
  ScoreJackStraw(dims = 1:35)
npc.num <- min(which(b1@reductions$pca@jackstraw$overall.p.values[ , 2] > 0.05), na.rm = T) - 1
if(npc.num <= 0) { npcs <- 35 }


#CellTrails nPC (https://github.com/dcellwanger/CellTrails/blob/master/R/dimred-methods.R)
.findSpectrum_def <- function(D, frac=100) {
  cs <- cumsum(diff(D))
  f <- ifelse(frac <= 1, length(D) * frac, frac)
  h <- head(cs, f)
  fit <- .linear_fit(x.in=seq_along(h), y.in=h)
  n <- min(which(diff(which(h > fit$y)) > 1)) + 1
  
  ggpl <- .plotSpectrum_def(list(frac=frac, n=n, cs=cs, fit=fit))
  print(ggpl)
  
  seq_len(n)
}