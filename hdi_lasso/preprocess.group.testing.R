preprocess.group.testing<-function (N, cov, conservative) 
{
  if (conservative) {
    NULL
  }
  else {
    zz <- mvrnorm(N, rep(0, ncol(cov)), cov)
    scale(zz, center = FALSE, scale = sqrt(diag(cov)))
  }
}