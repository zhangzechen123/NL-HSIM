## --------------------------------------------------------------
## Nyström kernel PCA core routine.
##
## Input:
##  - X:           n x p data matrix (rows are samples).
##  - n_total:     total sample size used in the simulation (used
##                 to determine the number of centers).
##  - seed:        random seed controlling the Nyström subsampling
##                 via k-means centers.
##  - center_ratio: proportion of "centers" used in the Nyström
##                  approximation (e.g. 0.05 * n_total).
##
## Steps:
##  1) Compute a bandwidth parameter b via stdv(X).
##  2) Use eff_kmeans() to select mm centers in the input space.
##  3) Construct Gaussian kernel matrices W (center-center) and
##     E (sample-center).
##  4) Perform an eigen-decomposition of W, extract positive
##     eigenvalues/vectors, and obtain the approximate feature
##     map G (Nyström embedding).
##  5) Center G in feature space and compute the centered kernel
##     matrix Kc = Gc Gc^T.
##  6) Perform SVD on Kc to obtain eigenvectors/eigenvalues of
##     the centered kernel.
##
## Output:
##  - evecs:    n x r matrix of kernel principal component directions.
##  - evals:    length-r vector of corresponding eigenvalues.
##  - bandwidth: bandwidth parameter b used in the Gaussian kernel.
## --------------------------------------------------------------
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
  ## Number of Nyström centers (at least 2 to avoid degeneracy)
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

