# ============================================================
# KPCA design matrix construction
#
# Three KPCA-based design matrices are constructed:
# - P  : pre-image–guided truncation using r_P components
# - PA : adaptive truncation within the P set
# - A  : fully adaptive truncation based on global eigenvalues
#
# Eigenvectors are scaled by the square root of their
# corresponding eigenvalues to preserve variance structure.
# ============================================================


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
    tol_preimg      = 1e-5,
    method          = c("ALL", "P", "A", "PA")
) {
  ## --------------------------------------------------------------------
  ## High-level wrapper for applying Nyström KPCA to one feature group.
  ##
  ## method:
  ##  - "P"  : only pre-image–guided truncation (PX_P, r_P)
  ##  - "A"  : only average-eigenvalue truncation (PX_A, r_A)
  ##  - "PA" : pre-image first, then adaptive truncation within first r_P (PX_PA, r_PA)
  ##  - "ALL": return PX_P, PX_PA, PX_A (backward-compatible default)
  ## --------------------------------------------------------------------
  method <- match.arg(method)
  
  n <- nrow(X_group)
  if (n == 0L || ncol(X_group) == 0L) {
    return(list(
      PX_P  = matrix(0, n, 0),
      PX_PA = matrix(0, n, 0),
      PX_A  = matrix(0, n, 0),
      r_P   = 0L,
      r_PA  = 0L,
      r_A   = 0L,
      RE    = numeric(0L)
    ))
  }
  
  ## Helper: build A-only design
  build_A_only <- function(evecs, evals) {
    thr_all <- mean(evals, na.rm = TRUE)
    idx_A   <- which(evals > thr_all)
    r_A     <- max(length(idx_A), 1L)
    r_A     <- min(r_A, length(evals))
    
    PX_A <- matrix(0, nrow(evecs), r_A)
    for (i in seq_len(r_A)) {
      PX_A[, i] <- evecs[, i] * sqrt(evals[i])
    }
    list(PX_A = PX_A, r_A = as.integer(r_A))
  }
  
  ## A-only: no pre-image
  if (method == "A") {
    spec <- nystrom_kpca_core(
      X       = X_group,
      n_total = n_total
    )
    
    if (is.null(spec$evecs) || is.null(spec$evals) || length(spec$evals) == 0L) {
      return(list(
        PX_P  = matrix(0, n, 0),
        PX_PA = matrix(0, n, 0),
        PX_A  = matrix(0, n, 0),
        r_P   = 0L,
        r_PA  = 0L,
        r_A   = 0L,
        RE    = numeric(0L)
      ))
    }
    Ares <- build_A_only(spec$evecs, spec$evals)
    return(list(
      PX_P  = matrix(0, n, 0),
      PX_PA = matrix(0, n, 0),
      PX_A  = Ares$PX_A,
      r_P   = 0L,
      r_PA  = 0L,
      r_A   = Ares$r_A,
      RE    = numeric(0L)
    ))
  }
  
  ## Pre-image selection (needed for P / PA / ALL)
  sel <- preimage(
    X_group = X_group,
    n_total = n_total,
    seed    = seed,
    nfold   = nfold
  )
  
  ## Nyström KPCA on the whole group
  spec <- nystrom_kpca_core(
    X       = X_group,
    n_total = n_total
  )
  
  if (is.null(spec$evecs) || is.null(spec$evals) || length(spec$evals) == 0L) {
    return(list(
      PX_P  = matrix(0, n, 0),
      PX_PA = matrix(0, n, 0),
      PX_A  = matrix(0, n, 0),
      r_P   = 0L,
      r_PA  = 0L,
      r_A   = 0L,
      RE    = numeric(0L)
    ))
  }
  
  ## Construct designs (P, PA, A) from the same KPCA spectrum
  des <- build_kpca_designs(
    evecs = spec$evecs,
    evals = spec$evals,
    r_P   = sel$r_P
  )
  
  ## Return only what is requested (stable output schema)
  if (method == "P") {
    return(list(
      PX_P  = des$PX_P,
      PX_PA = matrix(0, n, 0),
      PX_A  = matrix(0, n, 0),
      r_P   = des$r_P,
      r_PA  = 0L,
      r_A   = 0L,
      RE    = sel$RE
    ))
  }
  
  if (method == "PA") {
    return(list(
      PX_P  = matrix(0, n, 0),
      PX_PA = des$PX_PA,
      PX_A  = matrix(0, n, 0),
      r_P   = des$r_P,
      r_PA  = des$r_PA,
      r_A   = 0L,
      RE    = sel$RE
    ))
  }
  
  ## method == "ALL"
  c(des, list(RE = sel$RE))
}
