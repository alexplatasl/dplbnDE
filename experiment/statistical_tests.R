#!/usr/bin/env Rscript
# =============================================================================
# Statistical Analysis of Experiment Results
# Friedman test, Holm's post-hoc, Wilcoxon, Critical Difference diagrams
# =============================================================================

suppressPackageStartupMessages({
  library(dplbnDE)
})

OUT_DIR <- Sys.getenv("OUT_DIR", "experiment/results")
results <- read.csv(file.path(OUT_DIR, "full_results.csv"), stringsAsFactors = FALSE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# Remove WANBIA from DE comparisons (add back for accuracy comparison)
de_results <- results[results$algorithm != "wanbia", ]
algorithms <- unique(de_results$algorithm)
datasets <- unique(de_results$dataset)
structures <- unique(de_results$structure)

# =============================================================================
# Helper: Compute summary per dataset (worst CLL of 31 runs = robustness)
# =============================================================================
compute_summary <- function(df, metric, fun = min) {
  agg <- aggregate(as.formula(paste(metric, "~ dataset + structure + algorithm")),
                   data = df, FUN = fun, na.rm = TRUE)
  agg
}

# =============================================================================
# Friedman Test + Holm's Post-Hoc
# =============================================================================
friedman_analysis <- function(summary_df, metric_col, struct_name, metric_name) {
  cat("\n=== Friedman Test:", metric_name, "- Structure:", struct_name, "===\n")

  sub <- summary_df[summary_df$structure == struct_name, ]
  wide <- reshape(sub[, c("dataset", "algorithm", metric_col)],
                  idvar = "dataset", timevar = "algorithm",
                  direction = "wide")
  mat <- as.matrix(wide[, -1])
  colnames(mat) <- gsub(paste0(metric_col, "."), "", colnames(mat))

  # Rank per dataset (higher = better for CLL and accuracy)
  ranks <- t(apply(mat, 1, function(x) rank(-x)))
  avg_ranks <- colMeans(ranks, na.rm = TRUE)
  n <- nrow(mat)
  k <- ncol(mat)

  # Friedman statistic
  chi2 <- (12 * n / (k * (k + 1))) * (sum(avg_ranks^2) - k * (k + 1)^2 / 4)
  p_friedman <- pchisq(chi2, df = k - 1, lower.tail = FALSE)

  cat("Chi-squared:", round(chi2, 4), " p-value:", format.pval(p_friedman), "\n")
  cat("\nAverage Ranks:\n")
  sorted <- sort(avg_ranks)
  for (nm in names(sorted)) {
    cat(sprintf("  %-15s %.2f\n", nm, sorted[nm]))
  }

  # Holm's post-hoc (control = best ranked)
  control <- names(sorted)[1]
  others <- names(sorted)[-1]
  z_values <- numeric(length(others))
  p_raw <- numeric(length(others))

  for (j in seq_along(others)) {
    z_values[j] <- (avg_ranks[others[j]] - avg_ranks[control]) /
      sqrt(k * (k + 1) / (6 * n))
    p_raw[j] <- 2 * pnorm(-abs(z_values[j]))
  }

  # Sort by p-value for Holm correction
  ord <- order(p_raw)
  p_holm <- p.adjust(p_raw[ord], method = "holm")
  p_bonf <- p.adjust(p_raw[ord], method = "bonferroni")
  p_hoch <- p.adjust(p_raw[ord], method = "hochberg")

  cat("\nPost-hoc (control:", control, "):\n")
  cat(sprintf("%-20s %10s %10s %10s %10s\n", "Algorithm", "p-unadj", "p-Bonf", "p-Holm", "p-Hochberg"))
  for (j in seq_along(ord)) {
    sig <- ifelse(p_holm[j] < 0.05, "*", "")
    cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f %s\n",
                others[ord[j]], p_raw[ord[j]], p_bonf[j], p_holm[j], p_hoch[j], sig))
  }

  list(avg_ranks = avg_ranks, p_friedman = p_friedman,
       control = control, sorted_ranks = sorted)
}

# =============================================================================
# Run analysis for each metric × structure
# =============================================================================

# CLL analysis (worst case = min over 31 runs)
cll_worst <- compute_summary(de_results, "CLL", min)

# Training accuracy (mean over 31 runs)
acc_mean <- compute_summary(de_results, "train_accuracy", mean)

# CV accuracy (mean over 31 runs)
cv_mean <- compute_summary(de_results, "cv_accuracy", mean)

cat("\n######################################################################\n")
cat("# CLL ANALYSIS (worst case of 31 runs — robustness measure)\n")
cat("######################################################################\n")
for (s in structures) {
  friedman_analysis(cll_worst, "CLL", s, "CLL")
}

cat("\n######################################################################\n")
cat("# TRAINING ACCURACY ANALYSIS (mean of 31 runs)\n")
cat("######################################################################\n")
for (s in structures) {
  friedman_analysis(acc_mean, "train_accuracy", s, "Training Accuracy")
}

cat("\n######################################################################\n")
cat("# 5-FOLD CV ACCURACY ANALYSIS (mean of 31 runs)\n")
cat("######################################################################\n")
for (s in structures) {
  friedman_analysis(cv_mean, "cv_accuracy", s, "CV Accuracy")
}

# =============================================================================
# Summary table: Average rank across all structures
# =============================================================================
cat("\n######################################################################\n")
cat("# OVERALL AVERAGE RANKS (across all structures)\n")
cat("######################################################################\n")

overall_cll <- compute_summary(de_results, "CLL", min)
wide_all <- reshape(overall_cll[, c("dataset", "structure", "algorithm", "CLL")],
                    idvar = c("dataset", "structure"), timevar = "algorithm",
                    direction = "wide")
mat_all <- as.matrix(wide_all[, -(1:2)])
colnames(mat_all) <- gsub("CLL.", "", colnames(mat_all))
ranks_all <- t(apply(mat_all, 1, function(x) rank(-x)))
overall_avg <- sort(colMeans(ranks_all, na.rm = TRUE))

cat("\nCLL Overall Ranks:\n")
for (nm in names(overall_avg)) {
  cat(sprintf("  %d. %-15s %.2f\n", which(names(overall_avg) == nm), nm, overall_avg[nm]))
}

# =============================================================================
# Wilcoxon signed-rank tests (pairwise vs JADE with Archive)
# =============================================================================
cat("\n######################################################################\n")
cat("# WILCOXON SIGNED-RANK: each algorithm vs jade_archive\n")
cat("######################################################################\n")

for (s in structures) {
  cat("\nStructure:", s, "\n")
  sub <- cll_worst[cll_worst$structure == s, ]
  jade_vals <- sub$CLL[sub$algorithm == "jade_archive"]
  jade_ds <- sub$dataset[sub$algorithm == "jade_archive"]

  for (alg in setdiff(algorithms, "jade_archive")) {
    alg_vals <- sub$CLL[sub$algorithm == alg]
    alg_ds <- sub$dataset[sub$algorithm == alg]
    # Align by dataset
    common <- intersect(jade_ds, alg_ds)
    j <- jade_vals[match(common, jade_ds)]
    a <- alg_vals[match(common, alg_ds)]

    wt <- tryCatch(wilcox.test(a, j, paired = TRUE), error = function(e) list(p.value = NA))
    wins <- sum(a > j, na.rm = TRUE)
    losses <- sum(a < j, na.rm = TRUE)
    ties <- sum(a == j, na.rm = TRUE)
    cat(sprintf("  %-15s p=%.4f  W/L/T: %d/%d/%d %s\n",
                alg, wt$p.value, wins, losses, ties,
                ifelse(!is.na(wt$p.value) && wt$p.value < 0.05, "*", "")))
  }
}

# =============================================================================
# Save detailed summary table
# =============================================================================
cat("\nSaving summary tables...\n")

# Per-structure average ranks
for (s in structures) {
  sub <- cll_worst[cll_worst$structure == s, ]
  wide <- reshape(sub[, c("dataset", "algorithm", "CLL")],
                  idvar = "dataset", timevar = "algorithm", direction = "wide")
  write.csv(wide, file.path(OUT_DIR, "tables", paste0("cll_worst_", s, ".csv")),
            row.names = FALSE)
}

# Overall summary
all_results_summary <- results[results$algorithm != "wanbia", ]
overall <- aggregate(cbind(CLL, train_accuracy, cv_accuracy, n_evals, time) ~
                       algorithm + structure,
                     data = all_results_summary, FUN = mean, na.rm = TRUE)
write.csv(overall, file.path(OUT_DIR, "tables", "overall_means.csv"), row.names = FALSE)

cat("Done! Check", OUT_DIR, "\n")
