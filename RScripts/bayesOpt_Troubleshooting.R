# Bayes Opt Troubleshooting
# Contains internal functions from parBayesianOptimization to assist in troubleshooting


randParams <- function(boundsDT, rPoints, FAIL = TRUE) {
  
  # Attempt to procure rPoints unique parameter sets by lhs.
  attempt <- 1
  newPars <- data.table()
  poi <- rPoints
  
  while(attempt <= 100) {
    
    latinCube <- data.table(lhs::improvedLHS(n = poi, k = nrow(boundsDT)))
    
    setnames(latinCube, boundsDT$N)
    
    newPars <- unique(rbind(unMMScale(latinCube, boundsDT),newPars))
    
    if (nrow(newPars) == rPoints) break else poi <- rPoints-nrow(newPars)
    
    if (attempt >= 100 & FAIL) stop("Latin Hypercube Sampling could not produce the required distinct parameter sets. \nTry decreasing gsPoints or initPoints.")
    
    attempt <- attempt + 1
    
  }
  
  setnames(newPars, boundsDT$N)
  return(newPars)
  
}

unMMScale <- function(tabl, boundsDT) {
  
  umms <- lapply(boundsDT$N, function(x) {
    
    B <- boundsDT[get("N")==x,]
    
    n <- tabl[[x]]*B$R+B$L
    
    if (B$C == "integer") n <- round(n)
    
    return(n)
    
  })
  
  setDT(umms)
  if(!identical(names(tabl),boundsDT$N)) umms <- cbind(umms, tabl[,-boundsDT$N, with = F])
  setnames(umms, names(tabl))
  return(umms)
  
}

boundsToDT <- function(bounds) {
  data.table(
    N = names(bounds)
    , L = sapply(bounds, function(x) x[1])
    , U = sapply(bounds, function(x) x[2])
    , R = sapply(bounds, function(x) x[2]) - sapply(bounds, function(x) x[1])
    , C = sapply(bounds, function(x) class(x))
  )
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x.bounds = list(k.param = c(2L, 160L),
                res = c(0.0, 2.0))

boundsDT <- boundsToDT(x.bounds)
initGrid <- randParams(boundsDT = boundsDT, n.priors)
iter = 1:nrow(initGrid)
Params <- initGrid[get("iter"), ]

nestFun <- function() {
  do.call(what = function(...){testFun(balls = balls)}, args = as.list(Params))
}

testFun <- function(k.param, res, balls = NULL){
  print(balls)
  print(k.param)
  print(res)
}