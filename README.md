# NL-HSIM

**NL-HSIM** is a nonlinear high-dimensional set-based inference framework designed to detect group-level effects under complex correlation structures.  
This repository provides fully reproducible code for：
simulation studies (Cases 1–3)
real data analysis
and a minimal end-to-end tutorial of the full NL-HSIM workflow.

**Minimal tutorial:** 
For a fully runable example, please see [`tutorial_case1.md`](./tutorial_case1.md) 
The tutorial walks through the complete NL-HSIM pipeline, including:
data generation,
feature screening,
group-wise KPCA,
de-sparsified LASSO inference,
and final omnibus p-value calculation.

The repository follows a transparent script-based design: simulation settings, signal structures, and tuning parameters are specified explicitly in scripts, with no hidden defaults inside functions.

---

## 1. Requirements

- **R ≥ 4.0** (recommended)
- Required packages are loaded in [`environment.R`](./environment.R), which also defines core numerical utilities used throughout the project, including squared-distance computation and data-driven kernel bandwidth selection.

## 2. Quick start
Download the full repository from GitHub ([`Code -> Download ZIP)`] and extract it locally.
Open R (RStudio recommended).
Set the working directory to the project root, i.e., the folder containing [`README.md`](./README.md) and [`environment.R`](./environment.R).

Example:
```r
## Adjust this path to your local project folder
setwd("C:/Users/Downloads/NL-HSIM-main")

## Sanity checks
stopifnot(file.exists("environment.R"))
stopifnot(file.exists("generate_x.R"))
```
---

## 3. Repository structure

The scripts in [`hdi_lasso/`](./hdi_lasso) are utility functions used by [`lasso_model_omnibus.R`](./lasso_model_omnibus.R).  
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
├── benchmark_runtime_memory.R
├── block test of kpcs matrix.R
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
├── tutorial_case1.md
├── tutorial_case2.md
└── tutorial_case3.md

5 directories, 49 files
```
<!-- readme-tree end -->


---

## 4. Main entry files

The main workflow is organized around the following files:

- [`environment.R`](./environment.R)
  Loads required packages and defines shared numerical utilities.
- [`generate_x.R`](./generate_x.R)
  Generates simulated predictor data.
- [`kpca.R`](./kpca.R)
  Performs kernel PCA–based feature extraction.
- [`nystrom_kpca_core.R`](./nystrom_kpca_core.R)
  Implements the Nyström approximation for scalable KPCA.
- [`preimage.R`](./preimage.R)
  Contains pre-image reconstruction routines.
- [`lasso_model_omnibus.R`](./lasso_model_omnibus.R)
  Performs de-sparsified LASSO inference and omnibus p-value aggregation.
- [`run_pca_baseline.R`](./run_pca_baseline.R)
  Runs the linear PCA-based baseline.
- [`tutorial_case1.md`](./tutorial_case1.md), [`tutorial_case2.md`](./tutorial_case2.md), [`tutorial_case3.md`](./tutorial_case3.md)
  Runnable tutorial documents for the three simulation settings.

## 5. Minimal tutorial and reproducibility
A minimal vignette-style tutorial is provided in [`tutorial_case1.md`](./tutorial_case1.md)
, which demonstrates the complete NL-HSIM workflow from simulated SNP-type data generation to final omnibus inference.

This repository is designed for **full reproducibility**.
All simulation settings, signal structures, and tuning parameters are defined explicitly in scripts rather than hidden inside package-like wrappers.

## 6. Notes
This repository emphasizes transparency, reproducibility, and script-level control.
[`tutorial_case1.md`](./tutorial_case1.md) is the recommended entry point for new users.
Additional simulation settings and real-data analysis scripts are also included in the repository.

## 7. Citation / usage note
If you use this repository in your work, please cite the corresponding NL-HSIM manuscript.

---

