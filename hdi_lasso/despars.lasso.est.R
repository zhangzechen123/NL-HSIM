despars.lasso.est<-function (x, y, Z, betalasso) 
{
  b <- crossprod(Z, y - x %*% betalasso)/nrow(x) + betalasso
  if (ncol(b) == 1) 
    b <- as.vector(b)
  b
}