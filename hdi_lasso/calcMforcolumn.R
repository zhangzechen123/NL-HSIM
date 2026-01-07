calcMforcolumn<-function (x, j, lambdas) 
{
  glmnetfit <- glmnet(x[, -j], x[, j], lambda = lambdas)
  predictions <- predict(glmnetfit, x[, -j], s = lambdas)
  Zj <- x[, j] - predictions
  Znorms <- apply(Zj^2, 2, sum)
  Zxjnorms <- as.vector(crossprod(Zj, x[, j])^2)
  return(Znorms/Zxjnorms)
}