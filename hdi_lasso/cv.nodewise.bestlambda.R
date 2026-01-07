cv.nodewise.bestlambda<-function (lambdas, x, K = 10, parallel = FALSE, ncores = 8, oldschool = FALSE, 
          verbose = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  l <- length(lambdas)
  dataselects <- sample(rep(1:K, length = n))
  if (oldschool) {
    message("doing cv.nodewise.error oldschool")
    totalerr <- numeric(l)
    for (c in 1:p) {
      for (i in 1:K) {
        whichj <- dataselects == i
        glmnetfit <- glmnet(x[!whichj, -c, drop = FALSE], 
                            x[!whichj, c, drop = FALSE], lambda = lambdas)
        predictions <- predict(glmnetfit, x[whichj, -c, 
                                            drop = FALSE], s = lambdas)
        totalerr <- totalerr + apply((x[whichj, c] - 
                                        predictions)^2, 2, mean)
      }
    }
    totalerr <- totalerr/(K * p)
    stop("lambda.1se not implemented for oldschool cv.nodewise.bestlamba")
  }
  else {
    if (parallel) {
      totalerr <- mcmapply(cv.nodewise.err.unitfunction, 
                           c = 1:p, K = K, dataselects = list(dataselects = dataselects), 
                           x = list(x = x), lambdas = list(lambdas = lambdas), 
                           mc.cores = ncores, SIMPLIFY = FALSE, verbose = verbose, 
                           p = p)
    }
    else {
      totalerr <- mapply(cv.nodewise.err.unitfunction, 
                         c = 1:p, K = K, dataselects = list(dataselects = dataselects), 
                         x = list(x = x), lambdas = list(lambdas = lambdas), 
                         SIMPLIFY = FALSE, verbose = verbose, p = p)
    }
    err.array <- array(unlist(totalerr), dim = c(length(lambdas), 
                                                 K, p))
    err.mean <- apply(err.array, 1, mean)
    err.se <- apply(apply(err.array, c(1, 2), mean), 1, sd)/sqrt(K)
  }
  pos.min <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]
  stderr.lambda.min <- err.se[pos.min]
  list(lambda.min = lambda.min, lambda.1se = max(lambdas[err.mean < 
                                                           (min(err.mean) + stderr.lambda.min)]))
}