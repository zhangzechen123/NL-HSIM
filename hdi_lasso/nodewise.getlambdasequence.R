nodewise.getlambdasequence<-function (x) 
{
  nlambda <- 100
  p <- ncol(x)
  lambdas <- c()
  for (c in 1:p) {
    lambdas <- c(lambdas, glmnet(x[, -c], x[, c])$lambda)
  }
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}