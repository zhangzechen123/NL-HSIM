load('C:/research/paper/NL-HSIM/sorted_Group8.RData')
library(energy)
library(pheatmap)
source("C:/research/paper/NL-HSIM/environment.R")
source("C:/research/paper/NL-HSIM/nystrom_kpca_core.R")
source("C:/research/paper/NL-HSIM/preimage.R")
source("C:/research/paper/NL-HSIM/kpca.R")
source("C:/research/paper/NL-HSIM/lasso_model_omnibus.R")
set.seed(1)
gene_idx <- sample(seq_along(sorted_Group8), 10, replace = FALSE)
Group<-sorted_Group8[gene_idx]
cat("Sampled gene indices:\n")
print(gene_idx)
cat("Number of SNPs in sampled genes:\n")
print(sapply(Group, ncol))

# -----------------------------
# 2. Build real gene-wise KPC blocks
# -----------------------------
px_list_A <- vector("list", length(Group))

for (nn in seq_along(Group)) {
  cat("Processing gene block:", nn, "\n")
  
  b1 <- stdv(Group[[nn]])
  mm <- max(2, round(0.05 * nrow(Group[[nn]])))   # 更稳一点
  count <- 0
  max_try <- 100
  
  repeat {
    count <- count + 1
    tryCatch({
      center1 <- eff_kmeans(as.matrix(Group[[nn]]), mm, 5)
      if (!is.matrix(center1) || nrow(center1) < 2 || anyNA(center1)) {
        stop("Invalid center matrix")
      }
      break
    }, error = function(e) {
      message("Attempt ", count, ": ", e$message)
    })
    
    if (count >= max_try) {
      stop("eff_kmeans failed after 100 attempts.")
    }
  }
  
  W1 <- exp(-sqdist(t(as.matrix(center1)), t(as.matrix(center1))) / b1)
  E1 <- exp(-sqdist(t(as.matrix(Group[[nn]])), t(as.matrix(center1))) / b1)
  
  eig1 <- eigen(W1)
  vec1 <- eig1$vectors
  val1 <- eig1$values
  pidx1 <- which(val1 > 1e-6)
  
  inva11 <- diag(val1[pidx1]^(-0.5))
  G1 <- E1 %*% vec1[, pidx1] %*% inva11
  
  H <- diag(nrow(G1)) - matrix(1 / nrow(G1), nrow(G1), nrow(G1))
  G1c <- H %*% G1
  cen_K1 <- G1c %*% t(G1c)
  
  K_SVD1 <- svd(cen_K1)
  W11 <- K_SVD1$v
  D11 <- K_SVD1$d
  
  selected_pc1 <- which(D11 > mean(D11, na.rm = TRUE))
  num_of_PC_A1 <- max(length(selected_pc1), 1)
  
  px_list_A[[nn]] <- matrix(NA, nrow(Group[[nn]]), num_of_PC_A1)
  for (i in 1:num_of_PC_A1) {
    px_list_A[[nn]][, i] <- W11[, i] * sqrt(D11[i])
  }
}

cat("Retained KPC dimensions for real blocks:\n")
print(sapply(px_list_A, ncol))

# -----------------------------
# 3. Independent Gaussian baseline
#    Match final KPC block dimensions
# -----------------------------
null_dims <- sapply(px_list_A, ncol)
n <- nrow(px_list_A[[1]])

set.seed(1)
null_list <- lapply(null_dims, function(p) {
  matrix(rnorm(n * p), nrow = n, ncol = p)
})

cat("Retained KPC dimensions for null baseline:\n")
print(sapply(null_list, ncol))

# -----------------------------
# 4. Function: calculate three block-level matrices
# -----------------------------
calc_block_metrics <- function(block_list) {
  G <- length(block_list)
  
  block_mean_abs_cor <- matrix(NA, G, G)
  block_max_abs_cor  <- matrix(NA, G, G)
  block_dcor         <- matrix(NA, G, G)
  
  for (i in 1:G) {
    for (j in i:G) {
      if (i == j) {
        block_mean_abs_cor[i, j] <- 1
        block_max_abs_cor[i, j]  <- 1
        block_dcor[i, j]         <- 1
      } else {
        Xi <- as.matrix(block_list[[i]])
        Xj <- as.matrix(block_list[[j]])
        
        # pairwise column correlations
        Cij <- cor(Xi, Xj, use = "pairwise.complete.obs")
        vals <- abs(as.vector(Cij))
        vals <- vals[is.finite(vals)]
        
        block_mean_abs_cor[i, j] <- mean(vals)
        block_max_abs_cor[i, j]  <- max(vals)
        block_dcor[i, j]         <- dcor(Xi, Xj)
        
        block_mean_abs_cor[j, i] <- block_mean_abs_cor[i, j]
        block_max_abs_cor[j, i]  <- block_max_abs_cor[i, j]
        block_dcor[j, i]         <- block_dcor[i, j]
      }
    }
  }
  
  list(
    block_mean_abs_cor = block_mean_abs_cor,
    block_max_abs_cor  = block_max_abs_cor,
    block_dcor         = block_dcor
  )
}


# -----------------------------
# 5. Compute real and null matrices
# -----------------------------
real_metrics <- calc_block_metrics(px_list_A)
null_metrics <- calc_block_metrics(null_list)

block_mean_abs_cor      <- real_metrics$block_mean_abs_cor
block_max_abs_cor       <- real_metrics$block_max_abs_cor
block_dcor              <- real_metrics$block_dcor

null_block_mean_abs_cor <- null_metrics$block_mean_abs_cor
null_block_max_abs_cor  <- null_metrics$block_max_abs_cor
null_block_dcor         <- null_metrics$block_dcor

# -----------------------------
# 6. Comparison table (three metrics only)
# -----------------------------
compare_table <- data.frame(
  Metric = c(
    "Mean block mean |cor|",
    "Mean block max |cor|",
    "Mean block dCor"
  ),
  Real = c(
    mean(block_mean_abs_cor[upper.tri(block_mean_abs_cor)]),
    mean(block_max_abs_cor[upper.tri(block_max_abs_cor)]),
    mean(block_dcor[upper.tri(block_dcor)])
  ),
  Independent_Gaussian = c(
    mean(null_block_mean_abs_cor[upper.tri(null_block_mean_abs_cor)]),
    mean(null_block_max_abs_cor[upper.tri(null_block_max_abs_cor)]),
    mean(null_block_dcor[upper.tri(null_block_dcor)])
  )
)

print(compare_table)
# -----------------------------
# 7. Add row/col names
# -----------------------------
gene_lab <- paste0("Gene", 1:10)

rownames(block_mean_abs_cor)      <- gene_lab
colnames(block_mean_abs_cor)      <- gene_lab
rownames(null_block_mean_abs_cor) <- gene_lab
colnames(null_block_mean_abs_cor) <- gene_lab

rownames(block_max_abs_cor)       <- gene_lab
colnames(block_max_abs_cor)       <- gene_lab
rownames(null_block_max_abs_cor)  <- gene_lab
colnames(null_block_max_abs_cor)  <- gene_lab

rownames(block_dcor)              <- gene_lab
colnames(block_dcor)              <- gene_lab
rownames(null_block_dcor)         <- gene_lab
colnames(null_block_dcor)         <- gene_lab

# -----------------------------
# 8. Heatmap settings
# -----------------------------
bk <- seq(0, 1, length.out = 101)
hm_colors <- colorRampPalette(c("white", "dodgerblue3"))(100)

plot_hm_to_file <- function(mat, title_text, file_name) {
  png(file_name, width = 6, height = 5, units = "in", res = 300)
  
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = hm_colors,
    breaks = bk,
    border_color = "grey90",
    main = title_text,
    fontsize_number = 8,
    fontsize = 10
  )
  
  dev.off()
  cat("Saved:", file_name, "\n")
}

# -----------------------------
# 9. Output directory
# -----------------------------
outdir <- "C:/research/paper/review"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

write.csv(compare_table, "real_vs_independent_gaussian_comparison.csv", row.names = FALSE)

# -----------------------------
# 10. Save six heatmaps
# -----------------------------

# Pair 1: Mean absolute correlation
plot_hm_to_file(
  block_mean_abs_cor,
  "Real Data: Mean Absolute Cross-block Correlation",
  "FigS12_Real_MeanCor.png"
)
plot_hm_to_file(
  null_block_mean_abs_cor,
  "Null Baseline: Mean Absolute Cross-block Correlation",
  "FigS13_Null_MeanCor.png"
)

# Pair 2: Maximum absolute correlation
plot_hm_to_file(
  block_max_abs_cor,
  "Real Data: Maximum Absolute Cross-block Correlation",
  "FigS14_Real_MaxCor.png"
)
plot_hm_to_file(
  null_block_max_abs_cor,
  "Null Baseline: Maximum Absolute Cross-block Correlation",
  "FigS15_Null_MaxCor.png"
)

# Pair 3: Distance correlation
plot_hm_to_file(
  block_dcor,
  "Real Data: Block Distance Correlation (dCor)",
  "FigS16_Real_dCor.png"
)
plot_hm_to_file(
  null_block_dcor,
  "Null Baseline: Block Distance Correlation (dCor)",
  "FigS17_Null_dCor.png"
)

# -----------------------------
# 11. Save workspace (optional but recommended)
# -----------------------------
save.image(file = "matrix_instability_workspace.RData")