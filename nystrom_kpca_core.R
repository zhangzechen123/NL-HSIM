# ============================================================
# Nyström KPCA core routine
#
# This function performs kernel PCA using a Nyström
# approximation to reduce computational complexity.
#
# Key implementation details:
# - A fixed proportion (5%) of the total sample size is used
#   as Nyström centers, selected via an efficient k-means
#   algorithm.
# - The Gaussian kernel bandwidth is determined in a
#   data-driven manner based on average pairwise squared
#   distances.
#
# The output eigenvectors and eigenvalues are used to construct
# KPCA-based design matrices for downstream inference.
# ============================================================

nystrom_kpca_core <- function(
    X,
    n_total
) {
  n <- nrow(X)
  if (n == 0L || ncol(X) == 0L) {
    ## Degenerate case: no data
    return(list(
      evecs = matrix(0, 0, 0),
      evals = numeric(0L),
      bandwidth = NA_real_
    ))
  }
  ## Bandwidth parameter for the Gaussian kernel (data-driven)
  b  <- stdv(X)
  # Fixed Nyström sampling ratio (5%) adopted for stability
  # and computational efficiency across simulation settings
  mm <- max(2L, floor(0.05 * n_total))
  count   <- 0L
  max_try <- 100L
  ## Robust selection of centers via k-means (eff_kmeans),
  ## with multiple attempts to avoid rare numerical failures.
  repeat {
    count <- count + 1L
    ok <- TRUE
    center <- tryCatch(
      eff_kmeans(X, m = mm, Maxiter = 5L),
      error = function(e) {
        ok <<- FALSE
        NULL
      }
    )
    if (ok && is.matrix(center) && nrow(center) >= 2L && !anyNA(center)) {
      break
    }
    if (count >= max_try) {
      stop("eff_kmeans failed after ", max_try, " attempts in nystrom_kpca_core().")
    }
  }
  ## Gaussian kernel among centers (W) and between all samples and centers (E)
  W <- exp(-sqdist(t(as.matrix(center)), t(as.matrix(center))) / b)
  E <- exp(-sqdist(t(as.matrix(X)),      t(as.matrix(center))) / b)
  ## Eigen-decomposition of the center-center kernel
  eig <- eigen(W)
  vec <- eig$vectors
  val <- eig$values
  ## Keep only positive (non-negligible) eigenvalues
  pidx <- which(val > 1e-6)
  if (length(pidx) == 0L) {
    stop("No positive eigenvalues found in nystrom_kpca_core().")
  }
  
  inva <- diag(val[pidx]^(-0.5))
  G    <- E %*% vec[, pidx, drop = FALSE] %*% inva
  H    <- diag(nrow(G)) - matrix(1 / nrow(G), nrow(G), nrow(G))
  Gc   <- H %*% G
  Kc   <- Gc %*% t(Gc)
  
  K_SVD <- svd(Kc)
  
  list(
    evecs = K_SVD$v,
    evals = K_SVD$d,
    bandwidth = b
  )
}

