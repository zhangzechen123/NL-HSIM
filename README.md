# NL-HSIM

**NL-HSIM** is a nonlinear high-dimensional set-based inference framework designed to detect group-level effects under complex correlation structures.  
This repository provides fully reproducible code for simulation studies (Cases 1–3) and real data analysis.

The implementation follows a **transparent script-based design**: all simulation settings, signal structures, and tuning parameters are explicitly specified in scripts, with no hidden defaults inside functions.

---

## 1. Requirements

- **R ≥ 4.0** (recommended)
- Required R packages are automatically loaded in `environment.R`,which loads required packages and defines core numerical utilities, including squared-distance computation and data-driven kernel bandwidth selection. This file is essential for ensuring numerical stability and reproducibility across simulation settings.

### Quick start
- Please **download the entire `NL-HSIM` repository** and extract it to a local directory.
- Open R (RStudio is recommended) and set the working directory to the project root folder
   (i.e., the folder that contains `README.md` and `environment.R`).
```r
## After downloading "Code -> Download ZIP" and unzipping,
## GitHub typically creates a folder named "NL-HSIM-main".
## Set the working directory to the project root (the folder containing environment.R).
setwd("C:/Users/Downloads/NL-HSIM-main/NL-HSIM-main") # <-- example, adjust as needed

## Sanity check: these files should exist in the working directory
stopifnot(file.exists("environment.R"))
stopifnot(file.exists("generate_x.R"))
```
---

## 2. Repository structure

The scripts in `hdi_lasso/` are utility functions used by `lasso_model_omnibus.R`.  
They are loaded automatically via relative paths; users do **not** need to modify the working directory.

<!-- readme-tree start -->
```
.
├── .RData
├── .github
│   └── workflows
│       └── action_readme_tree.yml
├── .gitignore
├── LICENSE
├── README.md
├── aspu
│   ├── PowerUniv.R
│   ├── Sum.R
│   ├── SumSqU.R
│   ├── UminPd.R
│   └── aSPUd2.R
├── case2_example.R
├── environment.R
├── generate_x.R
├── hdi_lasso
│   ├── ART.A.R
│   ├── calcM.R
│   ├── calcMforcolumn.R
│   ├── calculate.Z.R
│   ├── cv.nodewise.bestlambda.R
│   ├── cv.nodewise.err.unitfunction.R
│   ├── cv.nodewise.totalerr.R
│   ├── despars.lasso.est.R
│   ├── do.initial.fit.R
│   ├── est.stderr.despars.lasso.R
│   ├── get.clusterGroupTest.function.R
│   ├── improve.lambda.pick.R
│   ├── initial.estimator.R
│   ├── lasso.proj.R
│   ├── nodewise.getlambdasequence.R
│   ├── p.adjust.wy.R
│   ├── prepare.data.R
│   ├── preprocess.group.testing.R
│   ├── ridge.proj.R
│   ├── sandwich.var.est.stderr.R
│   ├── score.getZforlambda.R
│   ├── score.getZforlambda.unitfunction.R
│   ├── score.nodewiselasso.R
│   └── score.rescale.R
├── kpca.R
├── lasso_model_omnibus.R
├── nystrom_kpca_core.R
├── preimage.R
├── run_pca_baseline.R
├── sorted_Group8.RData
├── tree.bak
├── treeview.sh
└── tutorial_case1.md

5 directories, 46 files
```
<!-- readme-tree end -->


---

## 3. Complete simulation example (Case 1)

This section provides a **fully runnable example** demonstrating the complete NL-HSIM pipeline.  
Users may modify only a few parameters (sample size, correlation level, signal strength) and directly run the code in **RStudio (recommended)**.

To evaluate empirical type I error and power, the same script can be repeated over multiple seeds (e.g., 1000 replications).

---

```r
## Step 0. User-defined parameters
n    <- 300   # sample size
p1   <- 100   # dimension of Group 1 (low-dimensional block)
p2   <- 700   # dimension of Group 2 (high-dimensional block)
rho  <- 0.8   # AR(1) correlation: 0.2, 0.5, 0.8
seed <- 1L    # simulation seed (can be looped over)
c1 <- 2       # signal strength for Group 1
c2 <- 0       # signal strength for Group 2

## Step 1. Load Modules
source("environment.R")
source("generate_x.R")
source("nystrom_kpca_core.R")
source("preimage.R")
source("kpca.R")
source("lasso_model_omnibus.R")
source("run_pca_baseline.R")
## Step 2. Generate Covariates
## Generate SNP-type predictors (Case 1 & Case 2)
## dat <- .gen_snp_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)
## Generate continuous predictors (Case 3)
## dat <- .gen_continuous_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)
## Case 2 uses real SNP data (sorted_Group8.RData)
## Here we use Case 1 (SNP data)
dat <- .gen_snp_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)
X1 <- dat$X1
X2 <- dat$X2
## Step 3. Specify Signal Structure
##Example: Case 1 – Nonlinear CWSM (Cosine-based)
hx1 <- apply(X1[, 1:15, drop = FALSE], 2, function(x) cos(pi * (x^2)))
hx1 <- rowSums(hx1)
hx2 <- apply(X2[, 1:50, drop = FALSE], 2, function(x) cos(pi * (x^2)))
hx2 <- rowSums(hx2)
## Step 4. Generate Outcome
set.seed(seed * 100L)
epsilon <- rnorm(n, mean = 0, sd = 1)
y <- epsilon + c1 * hx1 + c2 * hx2
## Step 5. Train–Test Split
X_all <- cbind(X1, X2)
set.seed(seed)
idx_train <- sample(seq_len(nrow(X_all)), size = nrow(X_all) / 2)
X_train <- X_all[idx_train, , drop = FALSE]
X_test  <- X_all[-idx_train, , drop = FALSE]
y_train <- y[idx_train]
y_test  <- y[-idx_train]
## Center outcome on test set
mn <- mean(y_test)
Y  <- y_test - mn
## Step 6. Feature Screening (Training → Test Projection)
A1 <- MFSIS(X_train, y_train, method = "DCSIS")
X_test_screened <- X_test[, A1, drop = FALSE]
## Split into Group 1 and Group 2
x1 <- as.matrix(X_test_screened[, A1 <= p1, drop = FALSE])
x2 <- as.matrix(X_test_screened[, A1  > p1, drop = FALSE])
## Step 7. KPCA Modules for Each Group
##Group 1
res1 <- kpca_module_group(
  X_group         = x1,
  n_total         = n,
  seed            = seed,
  nfold           = 10L,
  center_ratio    = 0.05,
  max_iter_preimg = 1000L,
  tol_preimg      = 1e-5,
  method = "ALL"
)

PX_P1  <- res1$PX_P    # Pre-image based KPCs
PX_PA1 <- res1$PX_PA   # Pre-image + average eigenvalue KPCs
PX_A1  <- res1$PX_A    # Average eigenvalue based KPCs

r_P1  <- res1$r_P
r_PA1 <- res1$r_PA
r_A1  <- res1$r_A
##Group 2
res2 <- kpca_module_group(
  X_group         = x2,
  n_total         = n,
  seed            = seed * 4L + 1L,
  nfold           = 10L,
  center_ratio    = 0.05,
  max_iter_preimg = 1000L,
  tol_preimg      = 1e-5,
  method = "ALL"
)

PX_P2  <- res2$PX_P
PX_PA2 <- res2$PX_PA
PX_A2  <- res2$PX_A

r_P2  <- res2$r_P
r_PA2 <- res2$r_PA
r_A2  <- res2$r_A

##Step 8. Combine Group-level Designs
PX_P  <- cbind(PX_P1,  PX_P2)
PX_PA <- cbind(PX_PA1, PX_PA2)
PX_A  <- cbind(PX_A1,  PX_A2)

##Step 9. LASSO + Omnibus Inference
res_lasso_P <- lasso_model_omnibus(
  X_design     = PX_P,
  Y            = Y,
  n_pc_g1      = r_P1,
  use_parallel = TRUE,
  ncores       = 2L
)

res_lasso_A <- lasso_model_omnibus(
  X_design     = PX_A,
  Y            = Y,
  n_pc_g1      = r_A1,
  use_parallel = TRUE,
  ncores       = 2L
)

res_lasso_PA <- lasso_model_omnibus(
  X_design     = PX_PA,
  Y            = Y,
  n_pc_g1      = r_PA1,
  use_parallel = TRUE,
  ncores       = 2L
)

##Step 10. Final Omnibus Combination (NL-HSIM(O))
O1 <- ACATO(c(res_lasso_A$Omnibus1,  res_lasso_PA$Omnibus1))
O2 <- ACATO(c(res_lasso_A$Omnibus2,  res_lasso_PA$Omnibus2))


##Linear Baseline (Linear-HSIM),Linear-HSIM is implemented in:
pca_out <- run_pca_baseline(
  x1            = x1,
  x2            = x2,
  Y             = Y,
  var_threshold = 0.85,
  use_parallel  = TRUE,
  ncores        = 2L
)

##Step 11. Collect and Save P-values
p_value_raw <- as.numeric(c(
  O1,
  O2,
  res_lasso_P$Omnibus1,
  res_lasso_P$Omnibus2,
  res_lasso_A$Omnibus1,
  res_lasso_A$Omnibus2,
  res_lasso_PA$Omnibus1,
  res_lasso_PA$Omnibus2,
  pca_out[[1]],
  pca_out[[2]]
))

print(p_value_raw)
```

---

## 4. Expected output

After running the code above, the main output is a numeric vector of p-values
stored in `p_value_raw`, corresponding to group-level hypothesis tests
under different NL-HSIM designs and the linear baseline.
The elements of p_value_raw are ordered as follows:
```r
 [1] 0.0000000 #NL-HSIM(O) omnibus p-value for Group 1
 [2] 0.2674444 #NL-HSIM(O) omnibus p-value for Group 2
 [3] 0.0000000 #NL-HSIM(P) p-value for Group 1 (pre-image based KPC selection)  
 [4] 0.4204262 #NL-HSIM(P) p-value for Group 2
 [5] 0.0000000 #NL-HSIM(A) p-value for Group 1 (average-eigenvalue based KPC selection)  
 [6] 0.3792151 #NL-HSIM(A) p-value for Group 2
 [7] 0.0000000 #NL-HSIM(PA) p-value for Group 1 (hybrid pre-image + average selection)  
 [8] 0.1981656 #NL-HSIM(PA) p-value for Group 2
 [9] 0.9198808 #Linear-HSIM p-value for Group 1 (PCA-based linear baseline)  
 [10]0.3251009 #Linear-HSIM p-value for Group 2

```
In this experiment, all NL-HSIM variants yield statistically significant p-values
for Group 1 (p < 0.05), while the corresponding p-values for Group 2 are non-significant.
This pattern is consistent across the pre-image (P), average-eigenvalue (A),
hybrid (PA), and omnibus (O) designs.

In contrast, the linear baseline (Linear-HSIM), which relies on standard PCA,
fails to detect the Group 1 effect, producing non-significant p-values for both groups.

These results indicate that NL-HSIM successfully captures the nonlinear group-level signal
embedded in Group 1, whereas linear PCA-based approaches lack sensitivity to such nonlinear
structures.

---
