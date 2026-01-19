## lasso_model_omnibus.R
##
## For a given design matrix X_design and response Y, this function:
##  1) runs de-sparsified lasso (lasso.proj),
##  2) computes group-wise MinP (Westfall–Young) and iART-A statistics,
##  3) combines them via ACATO to obtain omnibus p-values for Group 1 and Group 2.
##
## Additional safeguards:
##  - If ncol(X_design) < 3, it skips lasso.proj() to avoid glmnet errors
##    (glmnet requires at least 2 predictors in nodewise regressions).
##  - Parallel settings are automatically adjusted depending on the OS
##    (Windows is forced to run serially).
get_repo_root <- function() {
  # 1) If using RStudio project, assume working dir is project root
  if (file.exists("hdi_lasso")) return(".")
  # 2) Otherwise, try location of this file (when sourced)
  this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA)
  if (!is.na(this_file)) {
    cand <- dirname(this_file)
    if (file.exists(file.path(cand, "hdi_lasso"))) return(cand)
  }
  # 3) fallback: current dir
  return(".")
}

repo_root <- get_repo_root()
hdi_dir   <- file.path(repo_root, "hdi_lasso")

hdi_files <- c(
  "ART.A.R",
  "lasso.proj.R",
  "prepare.data.R",
  "calculate.Z.R",
  "score.nodewiselasso.R",
  "nodewise.getlambdasequence.R",
  "cv.nodewise.bestlambda.R",
  "cv.nodewise.err.unitfunction.R",
  "cv.nodewise.totalerr.R",
  "score.getZforlambda.R",
  "score.getZforlambda.unitfunction.R",
  "score.rescale.R",
  "initial.estimator.R",
  "do.initial.fit.R",
  "despars.lasso.est.R",
  "est.stderr.despars.lasso.R",
  "preprocess.group.testing.R",
  "get.clusterGroupTest.function.R",
  "sandwich.var.est.stderr.R",
  "improve.lambda.pick.R",
  "calcM.R",
  "calcMforcolumn.R",
  "p.adjust.wy.R",
  "ridge.proj.R"
)

for (f in hdi_files) {
  fp <- file.path(hdi_dir, f)
  if (!file.exists(fp)) stop("Missing HDI-Lasso file: ", fp)
  source(fp)
}


lasso_model_omnibus <- function(
    X_design,
    Y,
    n_pc_g1,
    use_parallel = TRUE,                   # allow parallelization on non-Windows systems
    ncores      = getOption("mc.cores", 2L)
) {
  p_one <- 1.0
  
  ## Case 1: empty design matrix
  if (is.null(X_design) || ncol(X_design) == 0L) {
    return(list(
      fit       = NULL,
      PG1       = p_one,
      PG2       = p_one,
      P_ART1    = p_one,
      P_ART2    = p_one,
      Omnibus1  = p_one,
      Omnibus2  = p_one
    ))
  }
  
  ## Case 2: too few columns (< 3), safeguard against glmnet errors
  if (ncol(X_design) < 3L) {
    warning("X_design has less than 3 columns; skipping lasso.proj and returning p=1.")
    return(list(
      fit       = NULL,
      PG1       = p_one,
      PG2       = p_one,
      P_ART1    = p_one,
      P_ART2    = p_one,
      Omnibus1  = p_one,
      Omnibus2  = p_one
    ))
  }
  
  ## Adjust parallel settings based on OS and available cores
  os_type <- .Platform$OS.type
  
  if (!use_parallel || os_type == "windows") {
    # On Windows or when parallelization is disabled, force serial execution
    use_parallel <- FALSE
    ncores       <- 1L
  } else {
    # On non-Windows systems (Linux/macOS), check the actual number of cores
    if (!requireNamespace("parallel", quietly = TRUE)) {
      use_parallel <- FALSE
      ncores       <- 1L
    } else {
      max_cores <- parallel::detectCores(logical = TRUE)
      if (is.na(max_cores) || max_cores < 2L) {
        use_parallel <- FALSE
        ncores       <- 1L
      } else {
        ncores <- min(as.integer(ncores), max_cores)
        if (ncores <= 1L) {
          use_parallel <- FALSE
          ncores       <- 1L
        }
      }
    }
  }
  
  ## Run de-sparsified lasso
  X_design <- scale(X_design)
  
  fit <- lasso.proj(
    x                   = as.matrix(X_design),
    y                   = Y,
    multiplecorr.method = "WY",
    parallel            = use_parallel,
    ncores              = ncores,
    robust              = TRUE
  )
  
  pvals <- fit$pval
  cov   <- fit$beta.cov
  L_tot <- length(pvals)
  
  if (L_tot == 0L) {
    return(list(
      fit       = fit,
      PG1       = p_one,
      PG2       = p_one,
      P_ART1    = p_one,
      P_ART2    = p_one,
      Omnibus1  = p_one,
      Omnibus2  = p_one
    ))
  }
  
  ## Split into Group 1 and Group 2
  L1 <- as.integer(n_pc_g1)
  if (is.na(L1) || L1 < 0L) {
    L1 <- 0L
  }
  if (L1 > L_tot) {
    L1 <- L_tot
  }
  L2 <- L_tot - L1
  
  ## ===================== MinP (Westfall–Young) =====================
  ## group 1
  if (L1 == 0L) {
    PG1 <- p_one
  } else {
    idx1 <- seq_len(L1)
    cov1 <- cov[idx1, idx1, drop = FALSE]
    PG1  <- min(p.adjust.wy(cov = cov1, pval = pvals[idx1]))
  }
  
  ## group 2
  if (L2 == 0L) {
    PG2 <- p_one
  } else {
    if (L1 == 0L) {
      idx2 <- seq_len(L_tot)
    } else {
      idx2 <- seq.int(L1 + 1L, L_tot)
    }
    cov2 <- cov[idx2, idx2, drop = FALSE]
    P2   <- pvals[idx2]
    PG2  <- min(p.adjust.wy(cov = cov2, pval = P2))
  }
  
  ## ===================== iART-A (within-group) + ACATO =====================
  ## group 1
  if (L1 == 0L) {
    P_ART1 <- p_one
  } else if (L1 == 1L) {
    P_ART1 <- pvals[1L]
  } else {
    P1 <- sort(pvals[seq_len(L1)])
    k1 <- 2L
    k2 <- L1
    P_arta_1 <- numeric(length = k2 - k1 + 1L)
    for (k in k1:k2) {
      P_arta_1[k - k1 + 1L] <- ART.A(P1, k = k, L = L1)[1L]
    }
    P_ART1 <- ACATO(P_arta_1)
    P_ART1 <- ifelse(P_ART1 == 1, 1 - 1 / (L1 - 1L), P_ART1)
  }
  
  ## group 2
  if (L2 == 0L) {
    P_ART2 <- p_one
  } else if (L2 == 1L) {
    if (L1 == 0L) {
      P_ART2 <- pvals[1L]
    } else {
      P_ART2 <- pvals[L1 + 1L]
    }
  } else {
    if (L1 == 0L) {
      idx2 <- seq_len(L_tot)
    } else {
      idx2 <- seq.int(L1 + 1L, L_tot)
    }
    P2     <- sort(pvals[idx2])
    L2_eff <- length(P2)
    k1     <- 2L
    k2     <- L2_eff
    P_arta_2 <- numeric(length = k2 - k1 + 1L)
    for (k in k1:k2) {
      P_arta_2[k - k1 + 1L] <- ART.A(P2, k = k, L = L2_eff)[1L]
    }
    P_ART2 <- ACATO(P_arta_2)
    P_ART2 <- ifelse(P_ART2 == 1, 1 - 1 / (L2_eff - 1L), P_ART2)
  }
  
  ## ===================== Omnibus = ACATO(MinP, iART-A) =====================
  Omnibus1 <- ACATO(c(PG1, P_ART1))
  Omnibus2 <- ACATO(c(PG2, P_ART2))
  
  list(
    fit       = fit,
    PG1       = PG1,
    PG2       = PG2,
    P_ART1    = P_ART1,
    P_ART2    = P_ART2,
    Omnibus1  = Omnibus1,
    Omnibus2  = Omnibus2
  )
}
