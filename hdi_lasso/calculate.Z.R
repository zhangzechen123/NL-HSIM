calculate.Z<-function (x, parallel, ncores, verbose, Z, do.ZnZ = FALSE, debug.verbose = FALSE) 
{
  if (is.null(Z)) {
    message("Nodewise regressions will be computed as no argument Z was provided.")
    message("You can store Z to avoid the majority of the computation next time around.")
    message("Z only depends on the design matrix x.")
    nodewiselasso.out <- score.nodewiselasso(x = x, wantTheta = FALSE, 
                                             parallel = parallel, ncores = ncores, cv.verbose = verbose || 
                                               debug.verbose, do.ZnZ = do.ZnZ, verbose = debug.verbose)
    Z <- nodewiselasso.out$out$Z
    scaleZ <- nodewiselasso.out$out$scaleZ
  }
  else {
    scaleZ <- rep(1, ncol(Z))
    if (!isTRUE(all.equal(rep(1, ncol(x)), colSums(Z * x)/nrow(x), 
                          tolerance = 10^-8))) {
      rescale.out <- score.rescale(Z = Z, x = x)
      Z <- rescale.out$Z
      scaleZ <- rescale.out$scaleZ
    }
  }
  list(Z = Z, scaleZ = scaleZ)
}