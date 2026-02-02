# ============================================================
# environment.R
#
# This file initializes the simulation environment.
# It loads all required packages and defines core numerical
# utilities used throughout the KPCA and inference pipeline.
#
# NOTE:
# - This file must be sourced before running any simulation.
# - Functions defined here (e.g., bandwidth selection,
#   squared distance computation) are critical for numerical
#   stability and reproducibility.
# ============================================================

## 1) Install and load required packages --------------------------------------
required_pkgs <- c(
  "Matrix", "glmnet", "foreach", "MASS", "iterators",
  "caret", "mvtnorm", "sumFREGAT", "parallel", "MFSIS", "svd"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

## 2) Utility functions: stdv, sqdist, eff_kmeans -----------------------------
stdv=function (data){
  n=dim(data)[1] 
  dim=dim(data)[2] 
  dis=matrix(numeric(),n,1)
  m =colMeans(data)  
  for (i in 1:n){
    dis[i,] = norm(as.matrix(m - data[i,]),'F')^2
  } 
  b = mean(dis) 
  b
}
Kernel_G=function(x,b){
  K=matrix(numeric(),nrow(x),nrow(x))
  K
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(x)) {
      K[i,j]=exp(-t((x[i,]-x[j,]))%*%(x[i,]-x[j,])/b)}} 
  K
}

sqdist=function(a,b){
  aa=colSums(a*a)
  bb=colSums(b*b)
  ab=t(a)%*%b
  d=abs(replicate(length(bb),aa)+t(replicate(length(aa),bb))-2*ab)
  d
}

eff_kmeans <- function(data, m, Maxiter) {
  n   <- nrow(data)
  dim <- ncol(data)
  dex <- sample(n)
  center <- as.matrix(data[dex[1:m], , drop = FALSE])
  
  for (i in 1:Maxiter) {
    nul <- rep(0, m)
    d   <- sqdist(t(center), t(data))
    xx  <- apply(d, 2, min)
    idx <- apply(d, 2, which.min)
    
    for (j in 1:m) {
      dex <- which(idx == j)
      l   <- length(dex)
      cltr <- as.matrix(data[dex, , drop = FALSE])
      if (l > 1) {
        center[j, ] <- colMeans(cltr)
      } else if (l == 1) {
        center[j, ] <- cltr
      } else {
        nul[j] <- 1
      }
    }
    dex    <- which(nul == 0)
    m      <- length(dex)
    center <- as.matrix(center[dex, , drop = FALSE])
  }
  center
}

