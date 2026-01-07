## main_simulation.R
##
## Simulation driver for NL-HSIM:
##  1) Generate X1, X2 (SNP or continuous),
##  2) Generate outcome Y under a chosen scenario,
##  3) Screen features via MFSIS (DCSIS),
##  4) Build KPCA-based design matrices (pre-image / average / hybrid),
##  5) Run de-sparsified lasso with omnibus tests for each design,
##  6) Construct final omnibus statistics O1/O2,
##  7) Build a PCA-based baseline and run the same omnibus tests,
##  8) Collect all p-values and write them to a CSV file.

setwd("C:/research/paper/NL-HSIM")
source("environment.R")
source("generate_x.R")
source("nystrom_kpca_core.R")
source("preimage.R")
source("kpca.R")
source("lasso_model_omnibus.R")

## ------------------------------------------------------------------
## 1) Basic settings
## ------------------------------------------------------------------

n    <- 400   # sample size: n = 400 or 600
p1   <- 100   # dimension of Group 1 (low-dimensional block)
p2   <- 700   # dimension of Group 2 (high-dimensional block)
rho  <- 0.7   # AR(1) correlation: 0.2, 0.5, 0.7
seed <- 1L    # simulation seed (can be looped over)

## ------------------------------------------------------------------
## 2) Generate covariates X1, X2
## ------------------------------------------------------------------
## Case 1: SNP data
dat <- .gen_snp_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)
X1  <- dat$X1
X2  <- dat$X2

## Case 3: continuous data
#dat <- .gen_continuous_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)
#X1  <- dat$X1
#X2  <- dat$X2

## ------------------------------------------------------------------
## 3) Specify nonlinear / linear signal structure (example: Case 1, nonlinear CWSM)
## ------------------------------------------------------------------
## Case 1: SNP data – nonlinear CWSM (cosine-based)
hx1 <- apply(X1[, 1:15, drop = FALSE], 2, function(x) cos(pi * (x^2)))
hx1 <- rowSums(hx1)
hx2 <- apply(X2[, 1:50, drop = FALSE], 2, function(x) cos(pi * (x^2)))
hx2 <- rowSums(hx2)
c1  <- 1
c2  <- 1
## Alternative signal strengths (commented out):
## c1 <- 0.1; c2 <- 0.2
## c1 <- 0.2; c2 <- 0.4
## c1 <- 0.3; c2 <- 0.6

## Case 1: SNP data – nonlinear DSSM
## hx1 <- cos(pi * (X1[, 1]^2))
## hx2 <- apply(X2[, 1:3, drop = FALSE], 2, function(x) cos(pi * (x^2)))
## hx2 <- rowSums(hx2)
## c1  <- 0.5; c2 <- 1
## c1  <- 1;   c2 <- 3
## c1  <- 1.5; c2 <- 5

## Case 1: SNP data – linear CWSM
## hx1 <- rowSums(X1[, 1:15, drop = FALSE])
## hx2 <- rowSums(X2[, 1:50, drop = FALSE])
## c1  <- 0.05;  c2 <- 0.05
## c1  <- 0.075; c2 <- 0.075
## c1  <- 0.10;  c2 <- 0.10

## Case 1: SNP data – linear DSSM
## hx1 <- X1[, 1]
## hx2 <- X2[, 1]
## c1  <- 0.10; c2 <- 0.05
## c1  <- 0.30; c2 <- 0.15
## c1  <- 0.50; c2 <- 0.25


## ------------------------------------------------------------------
## 3) Specify nonlinear / linear signal structure (example: Case 3, nonlinear CWSM)
## ------------------------------------------------------------------

## Case 3: continuous X – nonlinear CWSM
## hx1 <- rowSums(X1[, 1:15, drop = FALSE]^2)
## hx2 <- rowSums(X2[, 1:50, drop = FALSE]^2)
## c1  <- 0.06; c2 <- 0.07
## c1  <- 0.07; c2 <- 0.08
## c1  <- 0.08; c2 <- 0.09

## Case 3: continuous X – nonlinear DSSM
## hx1 <- X1[, 1]^2
## hx2 <- rowSums(X2[, 1:3, drop = FALSE]^2)
## c1  <- 0.30; c2 <- 0.40
## c1  <- 0.40; c2 <- 1.20
## c1  <- 0.50; c2 <- 2.00

## Case 3: continuous X – linear CWSM
## hx1 <- rowSums(X1[, 1:15, drop = FALSE])
## hx2 <- rowSums(X2[, 1:50, drop = FALSE])
## c1  <- 0.03; c2 <- 0.03
## c1  <- 0.04; c2 <- 0.04
## c1  <- 0.05; c2 <- 0.05

## Case 3: continuous X – linear DSSM
## hx1 <- X1[, 1]
## hx2 <- X2[, 1]
## c1  <- 0.10; c2 <- 0.10
## c1  <- 0.30; c2 <- 0.30
## c1  <- 0.50; c2 <- 0.50


## ------------------------------------------------------------------
## 4) Generate outcome Y
## ------------------------------------------------------------------
set.seed(seed * 100L)
epsilon <- rnorm(n, mean = 0, sd = 1)
y       <- epsilon + c1 * hx1 + c2 * hx2

## ------------------------------------------------------------------
## 5) Train–test split
## ------------------------------------------------------------------
X_all <- cbind(X1, X2)
set.seed(seed)
idx_train <- sample(seq_len(nrow(X_all)), size = nrow(X_all) / 2)
X_train   <- X_all[idx_train, , drop = FALSE]
X_test    <- X_all[-idx_train, , drop = FALSE]
y_train   <- y[idx_train]
y_test    <- y[-idx_train]

## Centered outcome on the test set
mn <- mean(y_test)
Y  <- y_test - mn

## ------------------------------------------------------------------
## 6) Feature screening on training data, then project to test data
## ------------------------------------------------------------------
A1 <- MFSIS(X_train, y_train, method = "DCSIS")
## Optional: uncomment for debugging
## print(A1)

X_test_screened <- X_test[, A1, drop = FALSE]

## Split screened test matrix into group 1 and group 2
x1 <- as.matrix(X_test_screened[, A1 <= p1, drop = FALSE])
x2 <- as.matrix(X_test_screened[, A1  > p1, drop = FALSE])


## ------------------------------------------------------------------
## 7) KPCA modules for Group 1 and Group 2
## ------------------------------------------------------------------
## Group 1
res1 <- kpca_module_group(
  X_group         = x1,
  n_total         = n,
  seed            = seed,      # e.g. gg
  nfold           = 10L,
  center_ratio    = 0.05,
  max_iter_preimg = 1000L,
  tol_preimg      = 1e-5
)
PX_P1   <- res1$PX_P
PX_PA1  <- res1$PX_PA
PX_A1   <- res1$PX_A
num_of_PC_P1  <- res1$r_P
num_of_PC_PA1 <- res1$r_PA
num_of_PC_A1  <- res1$r_A

## Group 2
res2 <- kpca_module_group(
  X_group         = x2,
  n_total         = n,
  seed            = seed * 4L + 1L,
  nfold           = 10L,
  center_ratio    = 0.05,
  max_iter_preimg = 1000L,
  tol_preimg      = 1e-5
)
PX_P2   <- res2$PX_P
PX_PA2  <- res2$PX_PA
PX_A2   <- res2$PX_A
num_of_PC_P2  <- res2$r_P
num_of_PC_PA2 <- res2$r_PA
num_of_PC_A2  <- res2$r_A

## Combine group-level designs
PX_P  <- cbind(PX_P1,  PX_P2)
PX_PA <- cbind(PX_PA1, PX_PA2)
PX_A  <- cbind(PX_A1,  PX_A2)

## Record group-wise KPC counts
r_P1  <- ncol(PX_P1)
r_P2  <- ncol(PX_P2)
r_PA1 <- ncol(PX_PA1)
r_PA2 <- ncol(PX_PA2)
r_A1  <- ncol(PX_A1)
r_A2  <- ncol(PX_A2)

## ------------------------------------------------------------------
## 8) Lasso + omnibus tests for three KPCA-based designs
## ------------------------------------------------------------------
## 1) Pre-image based KPC design
res_lasso_P <- lasso_model_omnibus(
  X_design     = PX_P,
  Y            = Y,
  n_pc_g1      = r_P1,
  use_parallel = TRUE,   # on Windows this will automatically be overridden to FALSE
  ncores       = 2L
)
## 2) Average eigenvalue based KPC design
res_lasso_A <- lasso_model_omnibus(
  X_design     = PX_A,
  Y            = Y,
  n_pc_g1      = r_A1,
  use_parallel = TRUE,
  ncores       = 2L
)
## 3) Hybrid pre-image + average eigenvalue design
res_lasso_PA <- lasso_model_omnibus(
  X_design     = PX_PA,
  Y            = Y,
  n_pc_g1      = r_PA1,
  use_parallel = TRUE,
  ncores       = 2L
)

## ------------------------------------------------------------------
## 9) Final omnibus statistics O1, O2 (combining A and PA designs)
## ------------------------------------------------------------------
## This follows your original script:
##   O1 = ACATO(Omnibus_A1, Omnibus_PA1)
##   O2 = ACATO(Omnibus_A2, Omnibus_PA2)

O1 <- ACATO(c(res_lasso_A$Omnibus1,  res_lasso_PA$Omnibus1))
O2 <- ACATO(c(res_lasso_A$Omnibus2,  res_lasso_PA$Omnibus2))


## ------------------------------------------------------------------
## 10) PCA-based baseline (comparison method)
## ------------------------------------------------------------------
## PCA for Group 1
if (ncol(x1) == 0L) {
  PC1          <- matrix(0, nrow(X_test_screened), 0)
  num_PC1_pca  <- 0L
} else {
  pc_ridge1   <- princomp(x1, cor = TRUE)
  lambda1     <- pc_ridge1$sdev^2
  var_per1    <- cumsum(lambda1) / sum(lambda1)
  id_var1     <- which(var_per1 > 0.85)[1L]
  n_pc1       <- 1:id_var1
  Z1          <- scale(x1) %*% pc_ridge1$loadings
  PC1         <- Z1[, n_pc1, drop = FALSE]
  num_PC1_pca <- length(n_pc1)
}

## PCA for Group 2
if (ncol(x2) == 0L) {
  PC2          <- matrix(0, nrow(X_test_screened), 0)
  num_PC2_pca  <- 0L
} else {
  pc_ridge2   <- princomp(x2, cor = TRUE)
  lambda2     <- pc_ridge2$sdev^2
  var_per2    <- cumsum(lambda2) / sum(lambda2)
  id_var2     <- which(var_per2 > 0.85)[1L]
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

## ------------------------------------------------------------------
## 11) Collect all p-values and write to CSV
## ------------------------------------------------------------------
## The order follows your original script:
##   (PG_P1, PG_P2, ART_P1, ART_P2, Omnibus_P1, Omnibus_P2,
##    PG_A1, PG_A2, ART_A1, ART_A2, Omnibus_A1, Omnibus_A2,
##    PG_PA1, PG_PA2, ART_PA1, ART_PA2, Omnibus_PA1, Omnibus_PA2,
##    O1, O2,
##    PG_pca1, PG_pca2, ART_pca1, ART_pca2, Omnibus_pca1, Omnibus_pca2)

p_value_raw <- as.numeric(c(
  ## Pre-image design
  res_lasso_P$Omnibus1, res_lasso_P$Omnibus2,
  ## Average design
  res_lasso_A$Omnibus1, res_lasso_A$Omnibus2,
  ## Hybrid (pre-image + average) design
  res_lasso_PA$Omnibus1, res_lasso_PA$Omnibus2,
  ## Final omnibus O1 / O2
  O1, O2,
  ## PCA baseline
  res_lasso_pca$Omnibus1, res_lasso_pca$Omnibus2
))

## Match original formatting: truncate extremely small p-values
p_value_clean     <- ifelse(p_value_raw <= 9e-4, 0, p_value_raw)
p_value_formatted <- formatC(p_value_clean, format = "f", digits = 5)
p_value           <- as.numeric(p_value_formatted)

## Build a small result object and write to CSV
Results <- cbind("output", seed, t(p_value))
colnames(Results)[1:2] <- c("Scenario", "Seed")

write.csv(
  Results,
  file      = paste0("Pre_image_KPCA_", seed, ".csv"),
  row.names = FALSE
)