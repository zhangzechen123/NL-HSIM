score.getZforlambda<-function (x, lambda, parallel = FALSE, ncores = 8, oldschool = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  Z <- matrix(numeric(n * p), n)
  if (oldschool) {
    message("doing getZforlambda oldschool")
    for (i in 1:p) {
      glmnetfit <- glmnet(x[, -i], x[, i])
      prediction <- predict(glmnetfit, x[, -i], s = lambda)
      Z[, i] <- x[, i] - prediction
    }
  }
  else {
    if (parallel) {
      Z <- mcmapply(score.getZforlambda.unitfunction, i = 1:p, 
                    x = list(x = x), lambda = lambda, mc.cores = ncores)
    }
    else {
      Z <- mapply(score.getZforlambda.unitfunction, i = 1:p, 
                  x = list(x = x), lambda = lambda)
    }
  }
  Z <- score.rescale(Z, x)
  return(Z)
}