cv.nodewise.totalerr<-function (c, K, dataselects, x, lambdas) 
{
  totalerr <- matrix(nrow = length(lambdas), ncol = K)
  for (i in 1:K) {
    whichj <- dataselects == i
    glmnetfit <- glmnet(x = x[!whichj, -c, drop = FALSE], 
                        y = x[!whichj, c, drop = FALSE], lambda = lambdas)
    predictions <- predict(glmnetfit, newx = x[whichj, -c, 
                                               drop = FALSE], s = lambdas)
    totalerr[, i] <- apply((x[whichj, c] - predictions)^2, 
                           2, mean)
  }
  totalerr
}