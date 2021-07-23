#Workflows

CellTrails.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.CellTrails(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  units(runtime) <- "mins"
  eval <- Test.method(results$object, method = "celltrails")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

CIDR.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.CIDR(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  units(runtime) <- "mins"
  eval <- Test.method(results, method = "cidr")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      NA
    )
  )
  
}

IKAP.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.IKAP(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  units(runtime) <- "mins"
  eval <- Test.method(results$object, method = "ikap")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

RaceID.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.RaceID(data)
  runtime <- Sys.time() - start.time
  units(runtime) <- "mins"
  eval <- Test.method(results$object, method = "raceid")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

SC3.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.SC3(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  units(runtime) <- "mins"
  eval <- Test.method(results$object, method = "sc3")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}

Seurat.flow <- function(data, expr.meas) {
  
  start.time <- Sys.time()
  results <- Run.Seurat(data, expr.meas = expr.meas)
  runtime <- Sys.time() - start.time
  units(runtime) <- "mins"
  eval <- Test.method(results$object, method = "seurat")
  
  return( 
    c(
      eval$ARI,
      eval$n.clusters,
      runtime,
      results$seed
    )
  )
  
}