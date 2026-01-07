score.getZforlambda.unitfunction<-function (i, x, lambda) 
{
  glmnetfit <- glmnet(x[, -i], x[, i])
  prediction <- predict(glmnetfit, x[, -i], s = lambda)
  return(x[, i] - prediction)
}