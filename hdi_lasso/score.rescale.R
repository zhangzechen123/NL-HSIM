score.rescale<-function (Z, x) 
{
  scaleZ <- diag(crossprod(Z, x))/nrow(x)
  Z <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z = Z, scaleZ = scaleZ))
}