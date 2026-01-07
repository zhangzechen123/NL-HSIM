lasso.proj<-function (x, y, family = "gaussian", standardize = TRUE, 
          multiplecorr.method = "holm", N = 10000, parallel = FALSE, 
          ncores = getOption("mc.cores", 2L), betainit = "cv lasso", 
          sigma = NULL, Z = NULL, verbose = FALSE, return.Z = FALSE, 
          suppress.grouptesting = FALSE, robust = FALSE, do.ZnZ = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  if (standardize) 
    sds <- apply(x, 2, sd)
  else sds <- rep(1, p)
  pdata <- prepare.data(x = x, y = y, standardize = standardize, 
                        family = family)
  x <- pdata$x
  y <- pdata$y
  if (family == "binomial") 
    sigma <- 1
  Zout <- calculate.Z(x = x, parallel = parallel, ncores = ncores, 
                      verbose = verbose, Z = Z, do.ZnZ = do.ZnZ)
  Z <- Zout$Z
  scaleZ <- Zout$scaleZ
  initial.estimate <- initial.estimator(betainit = betainit, 
                                        sigma = sigma, x = x, y = y)
  betalasso <- initial.estimate$beta.lasso
  sigmahat <- initial.estimate$sigmahat
  bproj <- despars.lasso.est(x = x, y = y, Z = Z, betalasso = betalasso)
  se <- est.stderr.despars.lasso(x = x, y = y, Z = Z, betalasso = betalasso, 
                                 sigmahat = sigmahat, robust = robust)
  scaleb <- 1/se
  bprojrescaled <- bproj * scaleb
  pval <- 2 * pnorm(abs(bprojrescaled), lower.tail = FALSE)
  cov2 <- crossprod(Z)
  pcorr <- if (multiplecorr.method == "WY") {
    p.adjust.wy(cov = cov2, pval = pval, N = N)
  }
  else if (multiplecorr.method %in% p.adjust.methods) {
    p.adjust(pval, method = multiplecorr.method)
  }
  else stop("Unknown multiple correction method specified")
  if (suppress.grouptesting) {
    group.testing.function <- NULL
    cluster.group.testing.function <- NULL
  }
  else {
    pre <- preprocess.group.testing(N = N, cov = cov2, conservative = FALSE)
    group.testing.function <- function(group, conservative = TRUE) {
      calculate.pvalue.for.group(brescaled = bprojrescaled, 
                                 group = group, individual = pval, conservative = conservative, 
                                 zz2 = pre)
    }
    cluster.group.testing.function <- get.clusterGroupTest.function(group.testing.function = group.testing.function, 
                                                                    x = x)
  }
  out <- list(pval = as.vector(pval), pval.corr = pcorr,beta.statistic=bprojrescaled,beta.cov=cov2,beta.cor=cov2cor(cov2),groupTest = group.testing.function, 
              clusterGroupTest = cluster.group.testing.function, sigmahat = sigmahat, 
              standardize = standardize, sds = sds, bhat = bproj/sds, 
              se = se/sds, betahat = betalasso/sds, family = family, 
              method = "lasso.proj", call = match.call())
  if (return.Z) 
    out <- c(out, list(Z = scale(Z, center = FALSE, scale = 1/scaleZ)))
  names(out$pval) <- names(out$pval.corr) <- names(out$bhat) <- names(out$sds) <- names(out$se) <- names(out$betahat) <- colnames(x)
  class(out) <- "hdi"
  return(out)
}