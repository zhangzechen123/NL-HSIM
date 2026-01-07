\# NL-HSIM:  Nonlinear high-dimensional statistical inference



This repository provides the R implementation for the simulation study in:

\*\*"Set-based nonlinear high-dimensional statistical inference model for detecting the genetic effect on complex diseases"\*\*.



\## Repository structure

\- `R/main\_simulation.R`: main driver script for simulations (data generation → screening → KPCA designs → de-sparsified lasso → omnibus p-values).

\- `R/environment.R`: required packages and utility functions.

\- `R/generate\_x.R`: data-generating mechanisms for continuous/SNP-like predictors.

\- `R/nystrom\_kpca\_core.R`: Nyström Gaussian KPCA core routine.

\- `R/preimage.R`: pre-image based selection of kernel PC number.

\- `R/kpca.R`: KPCA module wrapper and design construction (P / PA / A).

\- `R/lasso\_model\_omnibus.R`: de-sparsified lasso + MinP (WY) + iART-A + ACATO omnibus.



\## Requirements

\- R >= 4.x

\- Packages: Matrix, glmnet, caret, mvtnorm, MFSIS, parallel, etc.

> Tip: we recommend using `renv` to reproduce the exact package versions.



\## Quick start (one replicate)

```r

setwd("path/to/NL-HSIM")  # or open the project in RStudio and run from root

source("R/main\_simulation.R")



