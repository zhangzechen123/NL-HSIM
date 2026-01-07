## kpca.R
##
## This file implements the Nyström-based kernel PCA module used in NL-HSIM:
## - nystrom_kpca_core():    Nyström approximation of Gaussian kernel PCA
## - select_pc_preimage():   choose the number of kernel PCs via pre-image error
## - build_kpca_designs():   construct three KPCA-based designs (P, PA, A)
## - kpca_module_group():    wrapper to apply the above steps to one group of X


build_kpca_designs <- function(evecs, evals, r_P) {
  ## --------------------------------------------------------------------
  ## Construct three KPCA-based design matrices:
  ##
  ##  - PX_P:  uses the first r_P PCs (selected by pre-image criterion).
  ##  - PX_PA: "pre-image + average" hybrid design, using a subset of
  ##           the first r_P PCs whose eigenvalues are above the mean
  ##           eigenvalue among the first r_P.
  ##  - PX_A:  "average" design, using all PCs whose eigenvalues are
  ##           above the global mean eigenvalue.
  ##
  ## Each design matrix uses the standard KPCA embedding:
  ##   X -> evecs[, i] * sqrt(evals[i]), i in a selected index set.
  ##
  ## Input:
  ##  - evecs: n x d matrix of kernel eigenvectors.
  ##  - evals: length-d vector of kernel eigenvalues.
  ##  - r_P:   number of PCs selected by the pre-image criterion.
  ##
  ## Output:
  ##  - PX_P, PX_PA, PX_A: n x r design matrices for different KPCA
  ##    variants.
  ##  - r_P, r_PA, r_A:   effective dimensions used in each design.
  ## --------------------------------------------------------------------
  n <- nrow(evecs)
  d <- length(evals)
  
  if (d == 0L || n == 0L) {
    return(list(
      PX_P  = matrix(0, n, 0),
      PX_PA = matrix(0, n, 0),
      PX_A  = matrix(0, n, 0),
      r_P   = 0L,
      r_PA  = 0L,
      r_A   = 0L
    ))
  }
  
  r_P <- as.integer(r_P)
  if (is.na(r_P) || r_P <= 0L) {
    return(list(
      PX_P  = matrix(0, n, 0),
      PX_PA = matrix(0, n, 0),
      PX_A  = matrix(0, n, 0),
      r_P   = 0L,
      r_PA  = 0L,
      r_A   = 0L
    ))
  }
  
  r_P <- min(r_P, d)
  
  ## Hybrid PA: average eigenvalue within the first r_P eigenvalues
  if (r_P == 1L) {
    r_PA <- 1L
  } else {
    front_evals <- evals[1:r_P]
    thr_front   <- mean(front_evals, na.rm = TRUE)
    idx_PA      <- which(front_evals > thr_front)
    r_PA        <- max(length(idx_PA), 1L)
  }
  
  ## Pure average eigenvalue criterion on all eigenvalues
  thr_all <- mean(evals, na.rm = TRUE)
  idx_A   <- which(evals > thr_all)
  r_A     <- max(length(idx_A), 1L)
  
  r_PA <- min(r_PA, d)
  r_A  <- min(r_A,  d)
  
  PX_P  <- matrix(0, n, r_P)
  PX_PA <- matrix(0, n, r_PA)
  PX_A  <- matrix(0, n, r_A)
  
  for (i in seq_len(r_P)) {
    PX_P[, i] <- evecs[, i] * sqrt(evals[i])
  }
  for (i in seq_len(r_PA)) {
    PX_PA[, i] <- evecs[, i] * sqrt(evals[i])
  }
  for (i in seq_len(r_A)) {
    PX_A[, i] <- evecs[, i] * sqrt(evals[i])
  }
  
  list(
    PX_P  = PX_P,
    PX_PA = PX_PA,
    PX_A  = PX_A,
    r_P   = r_P,
    r_PA  = r_PA,
    r_A   = r_A
  )
}


kpca_module_group <- function(
    X_group,
    n_total,
    seed,
    nfold           = 10L,
    center_ratio    = 0.05,
    max_iter_preimg = 1000L,
    tol_preimg      = 1e-5
) {
  ## --------------------------------------------------------------------
  ## High-level wrapper for applying Nyström KPCA to one feature group.
  ##
  ## For a given group X_group (e.g., one omics block), this function:
  ##
  ##  1) Calls preimage() to choose the number of PCs r_P
  ##     by minimizing cross-validated pre-image error.
  ##  2) Calls nystrom_kpca_core() once on the full X_group to obtain
  ##     kernel eigenvectors/eigenvalues.
  ##  3) Calls build_kpca_designs() to construct three KPCA-based
  ##     design matrices: PX_P, PX_PA, PX_A.
  ##
  ## Input:
  ##  - X_group:       n x p matrix of one feature group.
  ##  - n_total:       total sample size in the simulation.
  ##  - seed:          base seed for KPCA and CV.
  ##  - nfold:         number of CV folds used in pre-image selection.
  ##  - center_ratio:  proportion of Nyström centers (0 < center_ratio <= 1).
  ##  - max_iter_preimg: maximum iterations in pre-image optimization.
  ##  - tol_preimg:    convergence threshold in pre-image optimization.
  ##
  ## Output:
  ##  - PX_P, PX_PA, PX_A: three n x r design matrices.
  ##  - r_P, r_PA, r_A:    corresponding numbers of PCs used.
  ##  - RE:                reconstruction error curve from pre-image
  ##                        selection (for diagnostic purposes).
  ## --------------------------------------------------------------------
  if (nrow(X_group) == 0L || ncol(X_group) == 0L) {
    return(list(
      PX_P  = matrix(0, nrow(X_group), 0),
      PX_PA = matrix(0, nrow(X_group), 0),
      PX_A  = matrix(0, nrow(X_group), 0),
      r_P   = 0L,
      r_PA  = 0L,
      r_A   = 0L,
      RE    = numeric(0L)
    ))
  }

  ## 1) Pre-image reconstruction error criterion to select r_P
  sel <- preimage(
    X_group         = X_group,
    n_total         = n_total,
    seed            = seed,
    nfold           = nfold
  )
  
  ## 2) Nyström KPCA on the whole group
  spec <- nystrom_kpca_core(
    X           = X_group,
    n_total     = n_total
  )
  
  ## 3) Construct P, A, and PA designs
  des <- build_kpca_designs(
    evecs = spec$evecs,
    evals = spec$evals,
    r_P   = sel$r_P
  )
  
  c(des, list(RE = sel$RE))
}
