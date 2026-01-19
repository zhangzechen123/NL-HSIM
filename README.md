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
<!-- readme-tree end -->


---

## 3. Complete simulation example (Case 1)

This section provides a **fully runnable example** demonstrating the complete NL-HSIM pipeline.  
Users may modify only a few parameters (sample size, correlation level, signal strength) and directly run the code in **RStudio (recommended)**.

To evaluate empirical type I error and power, the same script can be repeated over multiple seeds (e.g., 1000 replications).

---

### Step 0. User-defined parameters

```r
n    <- 400   # sample size
p1   <- 100   # dimension of Group 1 (low-dimensional block)
p2   <- 700   # dimension of Group 2 (high-dimensional block)
rho  <- 0.7   # AR(1) correlation: 0.2, 0.5, 0.8
seed <- 1L    # simulation seed (can be looped over)

c1 <- 0       # signal strength for Group 1
c2 <- 0       # signal strength for Group 2

### Step 1. Load required modules

```r
source("environment.R")
source("generate_x.R")
source("nystrom_kpca_core.R")
source("preimage.R")
source("kpca.R")
source("lasso_model_omnibus.R")
source("pca_baseline.R")

