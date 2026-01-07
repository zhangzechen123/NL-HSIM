## =========================================================
##  Data-generating mechanisms for NL-HSIM simulations
##  （continuous / SNP，lienar / nonlinear，CWSM / DSSM）
## =========================================================

#' Internal: AR(1) correlation matrix
#' @keywords internal
.ar1_sigma <- function(p, rho) {
  i <- seq_len(p)
  absdiff <- abs(outer(i, i, "-"))
  rho^absdiff
}


## ---------- Generate continuous X1, X2 ----------
#' Internal: generate continuous X1, X2
#' @keywords internal
.gen_continuous_X <- function(n, p1, p2, rho, seed) {
  
  Sigma1 <- .ar1_sigma(p1, rho)
  set.seed(seed)
  X1 <- MASS::mvrnorm(n, rep(0, p1), Sigma1)
  Sigma2 <- .ar1_sigma(p2, rho)
  set.seed(seed*2+1)
  X2 <- MASS::mvrnorm(n, rep(0, p2), Sigma2)
  list(X1 = X1, X2 = X2)
}


## ---------- Generate SNP-like X1, X2 ----------
#' Internal: generate discrete SNP-like X1, X2 by trichotomizing
#' @keywords internal
.gen_snp_X <- function(n, p1, p2, rho, seed) {
  set.seed(seed)
  Sigma1 <- .ar1_sigma(p1, rho)
  G1 <- MASS::mvrnorm(n, rep(0, p1), Sigma1)
  q1 <- quantile(G1, probs = c(1/3, 2/3))
  X1 <- apply(G1, 2, function(x) as.numeric(cut(x, breaks = c(-Inf, q1[1], q1[2], Inf), labels = c(0,1,2))))
  set.seed(seed*2+1)
  Sigma2 <- .ar1_sigma(p2, rho)
  G2 <- MASS::mvrnorm(n, rep(0, p2), Sigma2)
  q2 <- quantile(G2, probs = c(1/3, 2/3))
  X2 <- apply(G2, 2, function(x) as.numeric(cut(x, breaks = c(-Inf, q2[1], q2[2], Inf), labels = c(0,1,2))))
  list(X1 = as.matrix(X1), X2 = as.matrix(X2))
}
