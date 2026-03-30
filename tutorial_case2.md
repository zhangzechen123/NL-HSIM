# NL-HSIM Minimal Tutorial (Case 2: Real-data-driven simulation)

This tutorial provides a fully runnable end-to-end example of **Case 2** in NL-HSIM, where the simulation is driven by **real SNP data grouped by genes**.

Compared with Case 1, which uses synthetic grouped predictors, Case 2 starts from real SNP groups stored in `sorted_Group8.RData`, simulates linear or nonlinear phenotype signals on selected genes, and then evaluates NL-HSIM against existing set-based methods such as **SKAT** and **aSPU**.

This tutorial covers:
- loading real SNP groups,
- generating linear or nonlinear outcomes,
- train/test splitting,
- DC-SIS screening,
- mapping screened SNPs back to gene groups,
- group-wise KPCA,
- de-sparsified LASSO inference,
- omnibus aggregation with FDR correction,
- and comparison with SKAT / aSPU.

---

## Requirements

Before running this tutorial, please make sure that:

- the entire `NL-HSIM` repository has been downloaded,
- the working directory is set to the project root,
- and the following files are available in the project directory:

```r
source("environment.R")
source("nystrom_kpca_core.R")
source("preimage.R")
source("kpca.R")
source("lasso_model_omnibus.R")
source("sorted_Group8.RData")
source("aspu/aSPUd2.R")
source("aspu/Sum.R")
source("aspu/SumSqU.R")
source("aspu/UminPd.R")
source("aspu/PowerUniv.R")

required_pkgs <- c( "mvtnorm", "SKAT")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
load("sorted_Group8.RData")#Load real SNP groups

# 2) Simulate phenotype y
# -----------------------------
n_x <- 400
set.seed(gg)
idx <- sample(nrow(sorted_Group8[[1]]), n_x)

x11 <- as.matrix(sorted_Group8[[1]][idx, , drop = FALSE])
x22 <- as.matrix(sorted_Group8[[2]][idx, , drop = FALSE])

# signal strength
c1 <- 0
c2 <- 1
# ----- choose scenario -----
# nonlinear CWSM:
# hx1 <- rowSums(apply(x11[,1:15, drop=FALSE], 2, function(x) cos(pi*(x^2))))
# hx2 <- rowSums(apply(x22[,1:50, drop=FALSE], 2, function(x) cos(pi*(x^2))))
# nonlinear DSSM:
# hx1 <- cos(pi*(x11[,1]^2)); hx2 <- cos(pi*(x22[,1]^2))
# linear CWSM:
# hx1 <- rowSums(apply(x11[,1:15, drop=FALSE], 2, function(x) x))
# hx2 <- rowSums(apply(x22[,1:50, drop=FALSE], 2, function(x) x))
# linear DSSM:  
hx1 <- x11[, 1]
hx2 <- x22[, 1]

set.seed(gg * 100)
epsilon <- rnorm(n_x, 0, 1)
y <- epsilon + c1 * hx1 + c2 * hx2
# -----------------------------
# 3) Build full matrix & split train/test
# -----------------------------
all <- do.call(cbind, sorted_Group8)
all_idx <- all[idx, , drop = FALSE]

set.seed(gg)
SAM <- sample(x = nrow(all_idx), size = nrow(all_idx) / 2)

all_train <- as.matrix(all_idx[SAM,  , drop = FALSE])
all_test  <- as.matrix(all_idx[-SAM, , drop = FALSE])

y_train <- y[SAM]
y_test  <- y[-SAM]

# keep your original behavior: center test response
Y <- y_test - mean(y_test)

# -----------------------------
# 4) DCSIS screening on concatenated SNPs (train), apply to test
# -----------------------------
A1 <- MFSIS(all_train, y_train, method = "DCSIS")
sisdc <- all_test[, A1, drop = FALSE]
# -----------------------------
# 5) Map screened SNPs back to gene-wise list, drop empty groups
# -----------------------------
DC_test <- vector("list", length(sorted_Group8))
original_indices <- integer(0)

for (i in seq_along(sorted_Group8)) {
  keep_cols <- colnames(sisdc) %in% colnames(sorted_Group8[[i]])
  DC_test[[i]] <- sisdc[, keep_cols, drop = FALSE]
  if (ncol(DC_test[[i]]) > 0) original_indices <- c(original_indices, i)
}

empty <- which(sapply(DC_test, function(x) is.null(x) || ncol(x) == 0))
Group <- DC_test[-empty]

# where original gene1/gene2 land after dropping empties
idx_gene1 <- which(original_indices == 1)
idx_gene2 <- which(original_indices == 2)
gene1_exists <- length(idx_gene1) > 0
gene2_exists <- length(idx_gene2) > 0

cat("[INFO] #groups after DCSIS:", length(Group), "\n")
cat("[INFO] gene1_exists:", gene1_exists, " gene2_exists:", gene2_exists, "\n")

# -----------------------------
# 6) KPCA per group via modules -> PX_PA / PX_A
# -----------------------------
px_list_PA <- vector("list", length(Group))
px_list_A  <- vector("list", length(Group))
num_of_PC_PA <- integer(length(Group))
num_of_PC_A  <- integer(length(Group))
for (nn in seq_along(Group)) {
  cat("[KPCA] group:", nn, " n =", nrow(Group[[nn]]), " p =", ncol(Group[[nn]]), "\n")
  kp <- kpca_module_group(
    X_group = as.matrix(Group[[nn]]),
    n_total = nrow(Group[[nn]]),   # match your old script behavior
    seed    = gg * 1000 + nn,      # key fix: different folds per group
    nfold   = 10L
  )
  px_list_PA[[nn]] <- kp$PX_PA
  px_list_A[[nn]]  <- kp$PX_A
  num_of_PC_PA[nn] <- ncol(kp$PX_PA)
  num_of_PC_A[nn]  <- ncol(kp$PX_A)
}

X_PA <- do.call(cbind, px_list_PA)
X_A  <- do.call(cbind, px_list_A)
cat("[INFO] Total PCs: PA =", ncol(X_PA), " A =", ncol(X_A), "\n")# -----------------------------
# 7) de-sparsified lasso inference (PA/A)
# -----------------------------
fit_PA <- lasso.proj(
  x = as.matrix(X_PA),
  y = Y,
  multiplecorr.method = "WY",
  parallel = TRUE,
  ncores = getOption("mc.cores", 2L),
  robust = TRUE
)

fit_A <- lasso.proj(
  x = as.matrix(X_A),
  y = Y,
  multiplecorr.method = "WY",
  parallel = TRUE,
  ncores = getOption("mc.cores", 2L),
  robust = TRUE
)

pval_px_PA <- fit_PA$pval
cov_PA     <- fit_PA$beta.cov

pval_px_A  <- fit_A$pval
cov_A      <- fit_A$beta.cov

# -----------------------------
# 8) Multi-group omnibus (MinP(WY) + iART-A + ACATO) + FDR
# -----------------------------
split_blocks <- function(M, sizes) {
  out <- vector("list", length(sizes))
  st <- 1
  for (i in seq_along(sizes)) {
    sz <- sizes[i]
    if (sz <= 0) {
      out[[i]] <- NULL
    } else {
      out[[i]] <- M[st:(st+sz-1), st:(st+sz-1), drop = FALSE]
      st <- st + sz
    }
  }
  out
}

split_vec <- function(v, sizes) {
  out <- vector("list", length(sizes))
  st <- 1
  for (i in seq_along(sizes)) {
    sz <- sizes[i]
    if (sz <= 0) {
      out[[i]] <- numeric(0)
    } else {
      out[[i]] <- v[st:(st+sz-1)]
      st <- st + sz
    }
  }
  out
}

group_omnibus <- function(pvals, cov_mat, sizes) {
  cov_list <- split_blocks(cov_mat, sizes)
  p_list   <- split_vec(pvals, sizes)
  
  PG <- rep(1, length(sizes))
  P_arta <- rep(1, length(sizes))
  
  for (i in seq_along(sizes)) {
    L <- sizes[i]
    if (L <= 0) {
      PG[i] <- 1; P_arta[i] <- 1
      next
    }
    
    # MinP + WY
    PG[i] <- min(p.adjust.wy(cov = cov_list[[i]], pval = p_list[[i]]))
    
    # iART-A + ACATO
    if (L == 1) {
      P_arta[i] <- p_list[[i]][1]
    } else {
      P1 <- sort(p_list[[i]])
      tmp <- sapply(2:L, function(k) ART.A(P1, k, L)[1])
      P_arta[i] <- ACATO(tmp)
      if (P_arta[i] == 1) P_arta[i] <- 1 - 1/(L-1)  # keep your original guard
    }
  }
  
  Omnibus <- mapply(function(a, b) ACATO(c(a, b)), PG, P_arta)
  list(PG = PG, P_arta = P_arta, Omnibus = Omnibus)
}

res_PA <- group_omnibus(pval_px_PA, cov_PA, num_of_PC_PA)
res_A  <- group_omnibus(pval_px_A,  cov_A,  num_of_PC_A)

Omnibus_PA   <- res_PA$Omnibus
Omnibus_A    <- res_A$Omnibus
Omnibus_A_PA <- mapply(function(a, b) ACATO(c(a, b)), Omnibus_PA, Omnibus_A)

# FDR across gene-groups
Omnibus_PA_fdr   <- p.adjust(Omnibus_PA,   method = "fdr")
Omnibus_A_fdr    <- p.adjust(Omnibus_A,    method = "fdr")
Omnibus_A_PA_fdr <- p.adjust(Omnibus_A_PA, method = "fdr")

# gene1/gene2 extraction (if exist after dropping empties)
Omnibus_1_PA <- ifelse(gene1_exists, Omnibus_PA_fdr[idx_gene1], 1)
Omnibus_2_PA <- ifelse(gene2_exists, Omnibus_PA_fdr[idx_gene2], 1)

Omnibus_1_A  <- ifelse(gene1_exists, Omnibus_A_fdr[idx_gene1], 1)
Omnibus_2_A  <- ifelse(gene2_exists, Omnibus_A_fdr[idx_gene2], 1)

Omnibus_1_A_PA <- ifelse(gene1_exists, Omnibus_A_PA_fdr[idx_gene1], 1)
Omnibus_2_A_PA <- ifelse(gene2_exists, Omnibus_A_PA_fdr[idx_gene2], 1)

cat("\n==== Results (FDR) ====\n")
cat("Gene1: PA =", Omnibus_1_PA, " A =", Omnibus_1_A, " A+PA =", Omnibus_1_A_PA, "\n")
cat("Gene2: PA =", Omnibus_2_PA, " A =", Omnibus_2_A, " A+PA =", Omnibus_2_A_PA, "\n")

####################sakt and aspu##############

source('aSPUd2.R')
source('Sum.R')
source('SumSqU.R')
source('UminPd.R')
source('PowerUniv.R')
library(mvtnorm)
library(SKAT)

p_SKAT_Gene<-c()
p_SKAT_Gene1<-c()
p_aSPU_Gene<-c()
for (bb in 1:length(sorted_Group8)) {
  print(bb)
  
  X_test<-sorted_Group8[[bb]]
  x_test<-as.matrix(X_test[idx,])
  out1 <- SKAT_Null_Model(y ~ 1, out_type="C")
  p_SKAT_Gene[bb]<-SKAT(x_test, out1)$p.value
  p_SKAT_Gene1[bb]<-SKAT(x_test, out1,kernel="IBS")$p.value
  
  out2 <- aSPUd2(as.vector(y), as.matrix(x_test), cov = NULL,model = "gaussian") 
  p_aSPU_Gene[bb]<-out2[4]  
}
p_SKAT_Gene_fdr<-p.adjust( p_SKAT_Gene,method='fdr')
p_SKAT_Gene1_fdr<-p.adjust( p_SKAT_Gene1,method='fdr')
p_aSPU_Gene_fdr<-p.adjust( p_aSPU_Gene,method='fdr')


p_value_raw <- as.numeric(c(
  Omnibus_1_PA, Omnibus_2_PA,
  Omnibus_1_A, Omnibus_2_A,
  Omnibus_1_A_PA, Omnibus_2_A_PA,
  p_SKAT_Gene_fdr[1], p_SKAT_Gene_fdr[2],
  p_SKAT_Gene1_fdr[1], p_SKAT_Gene1_fdr[2],
  p_aSPU_Gene_fdr[1], p_aSPU_Gene_fdr[2]
))  # 保留小数点后6位，防止变成科学计数法

p_value_clean <- ifelse(p_value_raw <= 9e-4, 0, p_value_raw)
p_value_formatted <- formatC(p_value_clean, format = "f", digits = 5)
p_value <- as.numeric(p_value_formatted)







