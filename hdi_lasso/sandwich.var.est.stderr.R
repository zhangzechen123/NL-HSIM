sandwich.var.est.stderr<-function (x, y, betainit, Z) 
{
  n <- nrow(x)
  p <- ncol(x)
  if (!isTRUE(all.equal(rep(1, p), colSums(Z * x)/n, tolerance = 10^-8))) {
    rescale.out <- score.rescale(Z = Z, x = x)
    Z <- rescale.out$Z
  }
  if (length(betainit) > ncol(x)) {
    x <- cbind(rep(1, nrow(x)), x)
  }
  else {
    if (all(x[, 1] == 1)) {
    }
    else {
    }
  }
  eps.tmp <- as.vector(y - x %*% betainit)
  eps.tmp <- eps.tmp - mean(eps.tmp)
  sigmahatZ.direct <- sqrt(colSums(sweep(eps.tmp * Z, MARGIN = 2, 
                                         STATS = crossprod(eps.tmp, Z)/n, FUN = `-`)^2))
  sigmahatZ.direct/n
}