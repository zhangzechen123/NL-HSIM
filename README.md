# NL-HSIM

**NL-HSIM** is a nonlinear high-dimensional set-based inference framework designed to detect group-level effects under complex correlation structures.  
This repository provides fully reproducible code for simulation studies (Cases 1–3) and real data analysis.

The implementation follows a **transparent script-based design**: all simulation settings, signal structures, and tuning parameters are explicitly specified in scripts, with no hidden defaults inside functions.

---

## 1. Requirements

- **R ≥ 4.0** (recommended)
- Required R packages are automatically loaded in `environment.R`

### Important notes
- Please **download the entire `NL-HSIM` repository** and extract it to your local machine.
- Run all scripts from the **project root directory**.
- **Do not use `setwd()`**. All modules are loaded via relative paths.

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
└── treeview.sh

5 directories, 44 files
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
n    <- 400   # sample size
p1   <- 100   # dimension of Group 1 (low-dimensional block)
p2   <- 700   # dimension of Group 2 (high-dimensional block)
rho  <- 0.7   # AR(1) correlation: 0.2, 0.5, 0.8
seed <- 1L    # simulation seed (can be looped over)
c1 <- 0       # signal strength for Group 1
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
  tol_preimg      = 1e-5
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
  tol_preimg      = 1e-5
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


##Linear Baseline (Linear-HSIM)，Linear-HSIM is implemented in:
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
  pca_out[[3]],
  pca_out[[4]]
))
=======
\# NL-HSIM



\*\*NL-HSIM\*\* is a nonlinear high-dimensional set-based inference framework designed to detect group-level effects under complex correlation structures.  

This repository provides fully reproducible code for simulation studies (Cases 1–3) and real data analysis.



The implementation follows a \*\*transparent script-based design\*\*: all simulation settings, signal structures, and tuning parameters are explicitly specified in scripts, with no hidden defaults inside functions.



---



\## 1. Requirements



\- \*\*R ≥ 4.0\*\* (recommended)

\- Required R packages are automatically loaded in `environment.R`



\### Important notes

\- Please \*\*download the entire `NL-HSIM` repository\*\* and extract it to your local machine.

\- Run all scripts from the \*\*project root directory\*\*.

\- \*\*Do not use `setwd()`\*\*. All modules are loaded via relative paths.



---



\## 2. Repository structure



The scripts in `hdi\_lasso/` are utility functions used by `lasso\_model\_omnibus.R`.  

They are loaded automatically via relative paths; users do \*\*not\*\* need to modify the working directory.



<!-- readme-tree start -->

<!-- readme-tree end -->





---



\## 3. Complete simulation example (Case 1)



This section provides a \*\*fully runnable example\*\* demonstrating the complete NL-HSIM pipeline.  

Users may modify only a few parameters (sample size, correlation level, signal strength) and directly run the code in \*\*RStudio (recommended)\*\*.



To evaluate empirical type I error and power, the same script can be repeated over multiple seeds (e.g., 1000 replications).



---



```r

\## Step 0. User-defined parameters

n    <- 400   # sample size

p1   <- 100   # dimension of Group 1 (low-dimensional block)

p2   <- 700   # dimension of Group 2 (high-dimensional block)

rho  <- 0.7   # AR(1) correlation: 0.2, 0.5, 0.8

seed <- 1L    # simulation seed (can be looped over)

c1 <- 0       # signal strength for Group 1

c2 <- 0       # signal strength for Group 2



\## Step 1. Load Modules

source("environment.R")

source("generate\_x.R")

source("nystrom\_kpca\_core.R")

source("preimage.R")

source("kpca.R")

source("lasso\_model\_omnibus.R")

source("run\_pca\_baseline.R")

\## Step 2. Generate Covariates

\## Generate SNP-type predictors (Case 1 \& Case 2)

\## dat <- .gen\_snp\_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)

\## Generate continuous predictors (Case 3)

\## dat <- .gen\_continuous\_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)

\## Case 2 uses real SNP data (sorted\_Group8.RData)

\## Here we use Case 1 (SNP data)

dat <- .gen\_snp\_X(n = n, p1 = p1, p2 = p2, rho = rho, seed = seed)

X1 <- dat$X1

X2 <- dat$X2

\## Step 3. Specify Signal Structure

\##Example: Case 1 – Nonlinear CWSM (Cosine-based)

hx1 <- apply(X1\[, 1:15, drop = FALSE], 2, function(x) cos(pi \* (x^2)))

hx1 <- rowSums(hx1)

hx2 <- apply(X2\[, 1:50, drop = FALSE], 2, function(x) cos(pi \* (x^2)))

hx2 <- rowSums(hx2)

\## Step 4. Generate Outcome

set.seed(seed \* 100L)

epsilon <- rnorm(n, mean = 0, sd = 1)

y <- epsilon + c1 \* hx1 + c2 \* hx2

\## Step 5. Train–Test Split

X\_all <- cbind(X1, X2)

set.seed(seed)

idx\_train <- sample(seq\_len(nrow(X\_all)), size = nrow(X\_all) / 2)

X\_train <- X\_all\[idx\_train, , drop = FALSE]

X\_test  <- X\_all\[-idx\_train, , drop = FALSE]

y\_train <- y\[idx\_train]

y\_test  <- y\[-idx\_train]

\## Center outcome on test set

mn <- mean(y\_test)

Y  <- y\_test - mn

\## Step 6. Feature Screening (Training → Test Projection)

A1 <- MFSIS(X\_train, y\_train, method = "DCSIS")

X\_test\_screened <- X\_test\[, A1, drop = FALSE]

\## Split into Group 1 and Group 2

x1 <- as.matrix(X\_test\_screened\[, A1 <= p1, drop = FALSE])

x2 <- as.matrix(X\_test\_screened\[, A1  > p1, drop = FALSE])

\## Step 7. KPCA Modules for Each Group

\##Group 1

res1 <- kpca\_module\_group(

&nbsp; X\_group         = x1,

&nbsp; n\_total         = n,

&nbsp; seed            = seed,

&nbsp; nfold           = 10L,

&nbsp; center\_ratio    = 0.05,

&nbsp; max\_iter\_preimg = 1000L,

&nbsp; tol\_preimg      = 1e-5

)



PX\_P1  <- res1$PX\_P    # Pre-image based KPCs

PX\_PA1 <- res1$PX\_PA   # Pre-image + average eigenvalue KPCs

PX\_A1  <- res1$PX\_A    # Average eigenvalue based KPCs



r\_P1  <- res1$r\_P

r\_PA1 <- res1$r\_PA

r\_A1  <- res1$r\_A

\##Group 2

res2 <- kpca\_module\_group(

&nbsp; X\_group         = x2,

&nbsp; n\_total         = n,

&nbsp; seed            = seed \* 4L + 1L,

&nbsp; nfold           = 10L,

&nbsp; center\_ratio    = 0.05,

&nbsp; max\_iter\_preimg = 1000L,

&nbsp; tol\_preimg      = 1e-5

)



PX\_P2  <- res2$PX\_P

PX\_PA2 <- res2$PX\_PA

PX\_A2  <- res2$PX\_A



r\_P2  <- res2$r\_P

r\_PA2 <- res2$r\_PA

r\_A2  <- res2$r\_A



\##Step 8. Combine Group-level Designs

PX\_P  <- cbind(PX\_P1,  PX\_P2)

PX\_PA <- cbind(PX\_PA1, PX\_PA2)

PX\_A  <- cbind(PX\_A1,  PX\_A2)



\##Step 9. LASSO + Omnibus Inference

res\_lasso\_P <- lasso\_model\_omnibus(

&nbsp; X\_design     = PX\_P,

&nbsp; Y            = Y,

&nbsp; n\_pc\_g1      = r\_P1,

&nbsp; use\_parallel = TRUE,

&nbsp; ncores       = 2L

)



res\_lasso\_A <- lasso\_model\_omnibus(

&nbsp; X\_design     = PX\_A,

&nbsp; Y            = Y,

&nbsp; n\_pc\_g1      = r\_A1,

&nbsp; use\_parallel = TRUE,

&nbsp; ncores       = 2L

)



res\_lasso\_PA <- lasso\_model\_omnibus(

&nbsp; X\_design     = PX\_PA,

&nbsp; Y            = Y,

&nbsp; n\_pc\_g1      = r\_PA1,

&nbsp; use\_parallel = TRUE,

&nbsp; ncores       = 2L

)



\##Step 10. Final Omnibus Combination (NL-HSIM(O))

O1 <- ACATO(c(res\_lasso\_A$Omnibus1,  res\_lasso\_PA$Omnibus1))

O2 <- ACATO(c(res\_lasso\_A$Omnibus2,  res\_lasso\_PA$Omnibus2))





\##Linear Baseline (Linear-HSIM)，Linear-HSIM is implemented in:

pca\_out <- run\_pca\_baseline(

&nbsp; x1            = x1,

&nbsp; x2            = x2,

&nbsp; Y             = Y,

&nbsp; var\_threshold = 0.85,

&nbsp; use\_parallel  = TRUE,

&nbsp; ncores        = 2L

)



\##Step 11. Collect and Save P-values

p\_value\_raw <- as.numeric(c(

&nbsp; O1,

&nbsp; O2,

&nbsp; res\_lasso\_P$Omnibus1,

&nbsp; res\_lasso\_P$Omnibus2,

&nbsp; res\_lasso\_A$Omnibus1,

&nbsp; res\_lasso\_A$Omnibus2,

&nbsp; res\_lasso\_PA$Omnibus1,

&nbsp; res\_lasso\_PA$Omnibus2,

&nbsp; pca\_out\[\[3]],

&nbsp; pca\_out\[\[4]]

))



>>>>>>> Stashed changes
