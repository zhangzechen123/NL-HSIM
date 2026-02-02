# ============================================================
# Pre-image reconstruction–based dimension selection
#
# This function selects the truncation level of KPCA
# components using cross-validated pre-image reconstruction
# errors.
#
# IMPORTANT:
# - This procedure is NOT intended for reconstructing
#   individual observations.
# - It is used solely as a criterion to select the number
#   of KPCA components (r_P).
#
# Reconstruction error is evaluated at the level of
# test-sample means to stabilize dimension selection.
# ============================================================

preimage <- function(
    X_group,
    n_total,
    seed,
    nfold = 10L
) {
  n <- nrow(X_group)
  if (n == 0L || ncol(X_group) == 0L) {
    return(list(
      r_P = 0L,
      RE  = numeric(0L)
    ))
  }
  ## Create CV folds on the sample indices
  set.seed(seed)
  folds <- caret::createFolds(seq_len(n), k = nfold)
  ## The maximum feasible rank is determined by the smallest training set 
  max_r <- n - max(vapply(folds, length, integer(1L)))
  if (max_r <= 0L) {
    return(list(
      r_P = 0L,
      RE  = numeric(0L)
    ))
  }
  
  RE_vec <- rep(NA_real_, max_r)
  ## Loop over candidate number of PCs (num)
  for (num in seq_len(max_r)) {
    Ei <- numeric(nfold)
    ## Cross-validation over folds
    for (m in seq_len(nfold)) {
      traindata <- as.matrix(X_group[-folds[[m]], , drop = FALSE])
      testdata  <- as.matrix(X_group[ folds[[m]], , drop = FALSE])
      ## Run Nyström KPCA on the training data
      spec <- nystrom_kpca_core(
        X           = traindata,
        n_total     = n_total
      )
      eigvector <- spec$evecs
      eigvalue  <- spec$evals
      b         <- spec$bandwidth
      ## If the requested num exceeds the available rank, skip
      if (num > length(eigvalue)) {
        Ei[m] <- NA_real_
        next
      }
      ## Reduced coordinates in feature space for the first 'num' PCs
      reduceX <- matrix(0, nrow(eigvector), num)
      for (i in seq_len(num)) {
        reduceX[, i] <- eigvector[, i] * sqrt(eigvalue[i])
      }
      
      iter <- 1000L
      N    <- nrow(testdata)
      gamma <- numeric(N)
      for (i in seq_len(N)) {
        gamma[i] <- t(as.matrix(eigvector[i, 1:num, drop = FALSE])) %*% colMeans(reduceX)
      }
      
      z <- colMeans(testdata)
      for (s in seq_len(iter)) {
        prez <- z
        z    <- as.numeric(z)
        xx   <- t(testdata) - z
        xx   <- xx^2
        xx   <- -colSums(xx) / b
        xx   <- exp(xx)
        xx   <- xx * gamma
        
        z <- xx %*% testdata / sum(xx)
        z <- t(z)
        if (is.nan(z[1])) {
          next
        }
        denom <- norm(z, "F")
        if (denom == 0 || is.nan(denom)) {
          break
        }
        if (norm(prez - z, "F") / denom < 1e-5) {
          break
        }
      }
      
      prez[is.nan(prez)] <- 0
      diff_vec <- colMeans(testdata) - prez
      Ei[m]    <- as.numeric(t(diff_vec) %*% diff_vec)
    }
    
    RE_vec[num] <- mean(Ei, na.rm = TRUE)
  }
  
  if (all(is.na(RE_vec))) {
    return(list(
      r_P = 0L,
      RE  = RE_vec
    ))
  }
  
  r_P_opt <- which.min(RE_vec)[1]
  list(
    r_P = as.integer(r_P_opt),
    RE  = RE_vec
  )
}
# Reconstruction error is computed between the mean of the
# test samples and the estimated pre-image, rather than
# pointwise reconstruction, to reduce variability
