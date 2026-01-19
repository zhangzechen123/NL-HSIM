# ==============================================================
# PCA baseline module (comparison method)
# - Input: x1, x2 (screened test matrices), Y (centered test outcome),
#          lasso_model_omnibus() must be available (already sourced).
# - Output: list(res_lasso_pca, num_PC1_pca, num_PC2_pca, PX_pca, PC1, PC2)
# ==============================================================
run_pca_baseline<-function(x1,
                      x2,
                      Y,
                      var_threshold = 0.85,
                      use_parallel = TRUE,
                      ncores = 2L){
  # ---------- PCA for Group 1 ----------
if (ncol(x1) == 0L) {
  PC1          <- matrix(0, nrow(x1), 0)
  num_PC1_pca  <- 0L
} else {
  pc_ridge1   <- princomp(x1, cor = TRUE)
  lambda1     <- pc_ridge1$sdev^2
  var_per1    <- cumsum(lambda1) / sum(lambda1)
  id_var1     <- which(var_per1 > var_threshold )[1L]
  if (is.na(id_var1)) id_var1 <- ncol(x1)  # if never exceeds threshold, take all PCs
  n_pc1       <- 1:id_var1
  Z1          <- scale(x1) %*% pc_ridge1$loadings
  PC1         <- Z1[, n_pc1, drop = FALSE]
  num_PC1_pca <- length(n_pc1)
}

## PCA for Group 2
if (ncol(x2) == 0L) {
  PC2          <- matrix(0, nrow(x2), 0)
  num_PC2_pca  <- 0L
} else {
  pc_ridge2   <- princomp(x2, cor = TRUE)
  lambda2     <- pc_ridge2$sdev^2
  var_per2    <- cumsum(lambda2) / sum(lambda2)
  id_var2     <- which(var_per2 > var_threshold)[1L]
  if (is.na(id_var2)) id_var2 <- ncol(x2)
  n_pc2       <- 1:id_var2
  Z2          <- scale(x2) %*% pc_ridge2$loadings
  PC2         <- Z2[, n_pc2, drop = FALSE]
  num_PC2_pca <- length(n_pc2)
}

## Combine PCA scores from Group 1 and Group 2
PX_pca <- cbind(PC1, PC2)

## Run lasso + omnibus tests on PCA design
res_lasso_pca <- lasso_model_omnibus(
  X_design     = PX_pca,
  Y            = Y,
  n_pc_g1      = num_PC1_pca,
  use_parallel = TRUE,
  ncores       = 2L
)
list(
  res_lasso_pca$Omnibus1, res_lasso_pca$Omnibus2
)
}