calcM<-function (x, lambdas, parallel, ncores) 
{
  if (parallel) {
    M <- mcmapply(FUN = calcMforcolumn, x = list(x = x), 
                  j = 1:ncol(x), lambdas = list(lambdas = lambdas), 
                  mc.cores = ncores)
  }
  else {
    M <- mapply(FUN = calcMforcolumn, j = 1:ncol(x), x = list(x = x), 
                lambdas = list(lambdas = lambdas))
  }
  M <- apply(M, 1, mean)
  return(M)
}