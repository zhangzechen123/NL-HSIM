est.stderr.despars.lasso<-function (x, y, Z, betalasso, sigmahat, robust, robust.div.fixed = FALSE) 
{
  if (robust) {
    stderr <- sandwich.var.est.stderr(x = x, y = y, Z = Z, 
                                      betainit = betalasso)
    n <- nrow(x)
    if (robust.div.fixed) 
      stderr <- stderr * n/(n - sum(betalasso != 0))
  }
  else {
    stderr <- (sigmahat * sqrt(diag(crossprod(Z))))/nrow(x)
  }
  stderr
}