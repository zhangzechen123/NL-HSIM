# NL-HSIM

**Minimal tutorial:** please see [`tutorial_case1.md`](./tutorial_case1.md) for a fully runnable end-to-end example of the NL-HSIM workflow (from data generation and feature screening to KPCA, de-sparsified LASSO inference, and final omnibus p-values).

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

## 3. Minimal tutorial

A fully runnable end-to-end tutorial is provided in [`tutorial_case1.md`](./tutorial_case1.md).

The tutorial includes:
- user-defined parameters,
- required module loading,
- SNP-type data generation,
- signal specification,
- train/test splitting,
- DC-SIS feature screening,
- group-wise KPCA,
- de-sparsified LASSO inference,
- final omnibus p-value calculation,
- and interpretation of the expected output.

This README keeps a concise overview of the repository, while the tutorial file serves as the minimal vignette-style entry point for new users.

## 4. Notes

- The repository follows a script-based design for transparency and reproducibility.
- Simulation settings and signal structures are explicitly defined in scripts.
- The tutorial file `tutorial_case1.md` is intended as the minimal entry point for new users.
- Real-data analysis scripts and additional simulation settings remain available in the repository.

---

