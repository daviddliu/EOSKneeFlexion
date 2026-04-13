#!/usr/bin/env Rscript
# Analysis 15: PCA + DAG-style adjustment for subgroup preop S1PT > 20 (elevated pelvic tilt)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(boot)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
EXCLUDE_LOW_VOLUME_SURGEONS <- FALSE

PT_THRESHOLD <- 20

db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

if (EXCLUDE_LOW_VOLUME_SURGEONS) {
  df <- filter_low_volume_surgeons(df, min_surgeon_cases = 10)
}

df$change_lordosis <- df$LAT6W_L1_S1 - df$LATpre_L1_S1

if ("LATpre_L1PA" %in% names(df) && "LATpre_T4PA" %in% names(df)) {
  df$LATpre_T4_L1_PA <- df$LATpre_L1PA - df$LATpre_T4PA
  cat("Calculated LATpre_T4_L1_PA = LATpre_L1PA - LATpre_T4PA\n")
} else {
  cat("WARNING: Could not find LATpre_L1PA or LATpre_T4PA columns.\n")
  df$LATpre_T4_L1_PA <- NA
}

age_var <- NULL
if ("demo_Age" %in% names(df)) {
  age_var <- "demo_Age"
  cat("Using demo_Age as age variable\n")
} else if ("demo_age" %in% names(df)) {
  age_var <- "demo_age"
  cat("Using demo_age as age variable\n")
} else if ("age" %in% names(df)) {
  age_var <- "age"
  cat("Using age as age variable\n")
} else if ("Age" %in% names(df)) {
  age_var <- "Age"
  cat("Using Age as age variable\n")
} else {
  cat("WARNING: Age variable not found. Will exclude from PCA / DAG adjustment.\n")
}

if ("LATpre_PI_LL" %in% names(df)) {
  df$preop_PI_LL <- df$LATpre_PI_LL
  cat("Using LATpre_PI_LL as preoperative PI-LL mismatch\n")
} else if ("PI" %in% names(df)) {
  df$preop_PI_LL <- df$PI - df$LATpre_L1_S1
  cat("Calculating preoperative PI-LL mismatch as PI - LATpre_L1_S1\n")
} else {
  cat("WARNING: Could not find PI or LATpre_PI_LL; using LAT6W_PI_LL as preop PI-LL proxy.\n")
  df$preop_PI_LL <- df$LAT6W_PI_LL
}

n_before <- nrow(df)
df <- df %>%
  filter(!is.na(LATpre_S1PT), LATpre_S1PT > PT_THRESHOLD)
n_after <- nrow(df)
cat(sprintf("\nSubgroup: preop S1PT > %g  →  %d rows (excluded %d without PT or PT ≤ %g)\n",
            PT_THRESHOLD, n_after, n_before - n_after, PT_THRESHOLD))

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANALYSIS 15: PCA (analysis5-style) + DAG adjustment (analysis8-style)\n")
cat("Subgroup: LATpre_S1PT > ", PT_THRESHOLD, "\n", sep = "")
cat(rep("=", 80), "\n\n", sep = "")

# ==============================================================================
# Part A: PCA confounder analysis (mirrors analysis5_PCA_confounder_analysis.R)
# ==============================================================================

if (!is.null(age_var)) {
  df_clean <- df %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
           !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) &
           !is.na(LATpre_S1PT) & !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA) &
           !is.na(.data[[age_var]]))
} else {
  df_clean <- df %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
           !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) &
           !is.na(LATpre_S1PT) & !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA))
}

cat("\n=== Part A: PCA-Based Confounder Analysis ===\n")
cat("Complete cases for PCA + KF + outcome:", nrow(df_clean), "\n\n")

if (nrow(df_clean) < 15) {
  stop("Too few patients for PCA/regression in S1PT > ", PT_THRESHOLD, " subgroup (n = ", nrow(df_clean), ").")
}

confounder_vars <- c(
  "LATpre_S1PI",
  "LATpre_L1_S1",
  "LATpre_L4_S1",
  "LATpre_T2_T12",
  "LATpre_S1PT",
  "LATpre_SVA_C2_S1",
  "LATpre_T4_L1_PA"
)
confounder_labels <- c(
  "S1PI (preop)",
  "Lordosis L1-S1",
  "Lordosis L4-S1",
  "Thoracic Kyphosis T2-T12",
  "Pelvic Tilt S1PT",
  "SVA",
  "T4-L1 PA"
)
if (!is.null(age_var)) {
  confounder_vars <- c(confounder_vars, age_var)
  confounder_labels <- c(confounder_labels, "Age")
}

X_confounders <- df_clean[, confounder_vars]
colnames(X_confounders) <- confounder_labels

cat("Correlation matrix of confounders (before PCA):\n")
print(round(cor(X_confounders), 3))
cat("\n")

X_scaled <- scale(X_confounders)
pca_result <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
pca_components <- pca_result$x
pca_loadings <- pca_result$rotation
pca_summary <- summary(pca_result)

variance_explained <- pca_summary$importance[2, ]
cumulative_variance <- pca_summary$importance[3, ]
eigenvalues <- pca_result$sdev^2
n_vars <- length(confounder_vars)

broken_stick <- function(n_vars0, n_components) {
  expected <- numeric(n_components)
  for (i in 1:n_components) {
    expected[i] <- sum(1 / (i:n_vars0))
  }
  expected / n_vars0
}
broken_stick_expected <- broken_stick(n_vars, length(eigenvalues))
n_components_kaiser <- sum(eigenvalues > 1)
n_components_broken_stick <- sum(eigenvalues > broken_stick_expected)

n_components_80 <- sum(cumulative_variance < 0.80) + 1
if (n_components_80 > length(cumulative_variance)) n_components_80 <- length(cumulative_variance)
n_components_90 <- sum(cumulative_variance < 0.90) + 1
if (n_components_90 > length(cumulative_variance)) n_components_90 <- length(cumulative_variance)

set.seed(123)
n_iterations <- 100
n_obs <- nrow(X_scaled)
random_eigenvalues <- matrix(0, nrow = n_iterations, ncol = n_vars)
for (iter in 1:n_iterations) {
  random_data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  random_pca <- prcomp(random_data, center = TRUE, scale. = TRUE)
  random_eigenvalues[iter, ] <- random_pca$sdev^2
}
parallel_threshold <- apply(random_eigenvalues, 2, quantile, probs = 0.95)
n_components_parallel <- sum(eigenvalues > parallel_threshold)

n_ic_max <- min(7, length(cumulative_variance))
ic_results <- data.frame(
  N_Components = seq_len(n_ic_max),
  Variance_Explained = cumulative_variance[seq_len(n_ic_max)] * 100,
  AIC = NA_real_,
  BIC = NA_real_,
  Knee_Flexion_P = NA_real_
)
for (n_comp in seq_len(n_ic_max)) {
  df_temp <- df_clean
  for (i in 1:n_comp) df_temp[[paste0("PC", i)]] <- pca_components[, i]
  pca_formula_temp <- as.formula(paste(
    "change_lordosis ~ LATpre_LL_KneeAngle +",
    paste(paste0("PC", 1:n_comp), collapse = " + ")
  ))
  model_temp <- lm(pca_formula_temp, data = df_temp)
  ic_results$AIC[n_comp] <- AIC(model_temp)
  ic_results$BIC[n_comp] <- BIC(model_temp)
  ic_results$Knee_Flexion_P[n_comp] <- summary(model_temp)$coefficients[2, 4]
}
optimal_aic <- which.min(ic_results$AIC)
optimal_bic <- which.min(ic_results$BIC)
n_components_use <- optimal_bic

cat(sprintf("PCA: BIC-optimal number of components = %d\n", n_components_use))
print(ic_results, row.names = FALSE)
cat("\n")

loadings_df <- as.data.frame(pca_loadings[, 1:n_components_use, drop = FALSE])
colnames(loadings_df) <- paste0("PC", 1:n_components_use)

df_pca <- df_clean
for (i in 1:n_components_use) df_pca[[paste0("PC", i)]] <- pca_components[, i]

model1 <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = df_pca)
if (!is.null(age_var)) {
  model2 <- lm(
    as.formula(paste(
      "change_lordosis ~ LATpre_LL_KneeAngle + LATpre_S1PI + LATpre_L1_S1 +",
      "LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 +",
      "LATpre_T4_L1_PA +", age_var
    )),
    data = df_pca
  )
} else {
  model2 <- lm(
    change_lordosis ~ LATpre_LL_KneeAngle + LATpre_S1PI + LATpre_L1_S1 +
      LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA,
    data = df_pca
  )
}
pca_formula <- as.formula(paste(
  "change_lordosis ~ LATpre_LL_KneeAngle +",
  paste(paste0("PC", 1:n_components_use), collapse = " + ")
))
model3 <- lm(pca_formula, data = df_pca)

summary1 <- summary(model1)
summary2 <- summary(model2)
summary3 <- summary(model3)

n_pred_m2 <- length(coef(model2)) - 1
comparison_df <- data.frame(
  Model = c("Simple (KF only)", "Multiple (original confounders)",
            paste0("Multiple (PCA, ", n_components_use, " PCs)")),
  Knee_Flexion_Coef = c(
    summary1$coefficients[2, 1],
    summary2$coefficients[2, 1],
    summary3$coefficients[2, 1]
  ),
  Knee_Flexion_P = c(
    summary1$coefficients[2, 4],
    summary2$coefficients[2, 4],
    summary3$coefficients[2, 4]
  ),
  R_squared = c(summary1$r.squared, summary2$r.squared, summary3$r.squared),
  Adj_R_squared = c(summary1$adj.r.squared, summary2$adj.r.squared, summary3$adj.r.squared),
  N_Predictors = c(1, n_pred_m2, 1 + n_components_use)
)

calculate_vif <- function(model) {
  X <- model.matrix(model)[, -1, drop = FALSE]
  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)
  for (i in seq_len(ncol(X))) {
    other_preds <- X[, -i, drop = FALSE]
    vif_model <- lm(X[, i] ~ other_preds)
    r_squared <- summary(vif_model)$r.squared
    vif_values[i] <- if (r_squared >= 0.9999) Inf else 1 / (1 - r_squared)
  }
  vif_values
}
vif_pca <- calculate_vif(model3)
cat("VIF (KF + PCA model):\n")
print(round(vif_pca, 3))
cat("\n")

cat("=== PCA model summaries (KF coefficient) ===\n")
cat(sprintf("Model 1 KF: beta = %.4f, p = %.4e, R2 = %.4f\n",
            summary1$coefficients[2, 1], summary1$coefficients[2, 4], summary1$r.squared))
cat(sprintf("Model 2 KF: beta = %.4f, p = %.4e, R2 = %.4f\n",
            summary2$coefficients[2, 1], summary2$coefficients[2, 4], summary2$r.squared))
cat(sprintf("Model 3 KF: beta = %.4f, p = %.4e, R2 = %.4f\n\n",
            summary3$coefficients[2, 1], summary3$coefficients[2, 4], summary3$r.squared))

print(comparison_df, row.names = FALSE)
cat("\n")

if (!dir.exists("results")) dir.create("results")

png("results/analysis15_scree_plot.png", width = 10, height = 6, units = "in", res = 300)
scree_df <- data.frame(
  Component = seq_along(eigenvalues),
  Eigenvalue = eigenvalues,
  Variance_Explained = variance_explained * 100,
  Cumulative_Variance = cumulative_variance * 100
)
print(
  ggplot(scree_df, aes(x = Component, y = Eigenvalue)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "blue", size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(x = "Principal Component", y = "Eigenvalue",
         title = sprintf("Analysis 15: Scree (S1PT > %g)", PT_THRESHOLD)) +
    theme_minimal()
)
dev.off()

png("results/analysis15_variance_explained.png", width = 10, height = 6, units = "in", res = 300)
variance_df <- data.frame(
  Component = seq_along(variance_explained),
  Variance = variance_explained * 100,
  Cumulative = cumulative_variance * 100
)
print(
  ggplot(variance_df, aes(x = Component)) +
    geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    labs(x = "Principal Component", y = "Variance explained (%)",
         title = sprintf("Analysis 15: PCA variance (S1PT > %g)", PT_THRESHOLD)) +
    theme_minimal()
)
dev.off()

png("results/analysis15_loadings_heatmap.png", width = 10, height = 8, units = "in", res = 300)
loadings_long <- loadings_df %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Component", values_to = "Loading")
print(
  ggplot(loadings_long, aes(x = Component, y = Variable, fill = Loading)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    labs(title = sprintf("Analysis 15: PCA loadings (S1PT > %g)", PT_THRESHOLD), fill = "Loading") +
    theme_minimal()
)
dev.off()

coef_comparison <- data.frame(
  Model = comparison_df$Model,
  Coefficient = comparison_df$Knee_Flexion_Coef,
  Lower_CI = c(confint(model1)[2, 1], confint(model2)[2, 1], confint(model3)[2, 1]),
  Upper_CI = c(confint(model1)[2, 2], confint(model2)[2, 2], confint(model3)[2, 2]),
  P_value = comparison_df$Knee_Flexion_P
)
png("results/analysis15_model_comparison.png", width = 10, height = 6, units = "in", res = 300)
print(
  ggplot(coef_comparison, aes(x = Model, y = Coefficient)) +
    geom_point(size = 3, color = "red") +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(y = "Knee flexion coefficient", x = NULL,
         title = sprintf("Analysis 15: KF coefficient by model (S1PT > %g)", PT_THRESHOLD)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
)
dev.off()

pca_formula_no_kf <- as.formula(paste("change_lordosis ~", paste(paste0("PC", 1:n_components_use), collapse = " + ")))
covariates_model_pca <- lm(pca_formula_no_kf, data = df_pca)
df_pca$resid_from_pca_covariates <- residuals(covariates_model_pca)
t_stat_kf_pca <- summary3$coefficients[2, 3]
df_residual_pca <- summary3$df[2]
r_partial_pca <- t_stat_kf_pca / sqrt(t_stat_kf_pca^2 + df_residual_pca)
r2_partial_pca <- summary3$r.squared - summary(covariates_model_pca)$r.squared

png("results/analysis15_residual_plot.png", width = 10, height = 8, units = "in", res = 300)
print(
  ggplot(df_pca, aes(x = LATpre_LL_KneeAngle, y = resid_from_pca_covariates)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Preoperative knee flexion",
      y = paste0("Residual ΔLL after PC1–PC", n_components_use),
      title = sprintf("Analysis 15: KF vs residual ΔLL | PCA (S1PT > %g)", PT_THRESHOLD)
    ) +
    theme_minimal() +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0(
        "Partial r = ", round(r_partial_pca, 3),
        "\nPartial R² = ", round(r2_partial_pca, 3),
        "\np (KF in full PCA model) = ", formatC(summary3$coefficients[2, 4], format = "e", digits = 2)
      ),
      hjust = 1.05, vjust = 1.2, size = 3.8
    )
)
dev.off()

cat("Saved PCA plots under results/analysis15_*.png\n\n")

# ==============================================================================
# Part B: DAG adjustment set (mirrors analysis8_DAG.R core regression)
# ==============================================================================

adjustment_set <- c(
  "preop_PI_LL",
  "LATpre_L1_S1",
  "LATpre_L4_S1",
  "LATpre_T2_T12",
  "LATpre_S1PI",
  "LATpre_S1PT",
  "LATpre_SVA_C2_S1",
  "LATpre_T4_L1_PA"
)
if (!is.null(age_var)) adjustment_set <- c(adjustment_set, age_var)

cat("=== Part B: DAG-based adjustment set ===\n")
cat("Covariates:", paste(adjustment_set, collapse = ", "), "\n\n")

df_dag <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis))
complete_vars <- c("LATpre_LL_KneeAngle", "change_lordosis", adjustment_set)
df_complete <- df_dag %>%
  filter(complete.cases(select(., all_of(complete_vars))))

cat("Complete cases for DAG models:", nrow(df_complete), "\n\n")
if (nrow(df_complete) < 15) {
  stop("Too few complete cases for DAG models in subgroup.")
}

formula_adj <- as.formula(paste(
  "change_lordosis ~ LATpre_LL_KneeAngle +",
  paste(adjustment_set, collapse = " + ")
))
model_unadj <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = df_complete)
model_adj <- lm(formula_adj, data = df_complete)
formula_no_kf <- as.formula(paste("change_lordosis ~", paste(adjustment_set, collapse = " + ")))
model_no_kf <- lm(formula_no_kf, data = df_complete)

su_u <- summary(model_unadj)
su_a <- summary(model_adj)
su_nk <- summary(model_no_kf)

cat("--- Unadjusted: change_lordosis ~ KF ---\n")
cat(sprintf("  KF beta = %.4f, 95%% CI [%.4f, %.4f], p = %.4e, R² = %.4f\n",
            coef(model_unadj)[2], confint(model_unadj)[2, 1], confint(model_unadj)[2, 2],
            su_u$coefficients[2, 4], su_u$r.squared))

cat("\n--- Adjusted (DAG set): change_lordosis ~ KF + ... ---\n")
cat(sprintf("  KF beta = %.4f, 95%% CI [%.4f, %.4f], p = %.4e, R² = %.4f\n",
            coef(model_adj)[2], confint(model_adj)[2, 1], confint(model_adj)[2, 2],
            su_a$coefficients[2, 4], su_a$r.squared))

partial_r2_kf <- su_a$r.squared - su_nk$r.squared
cat(sprintf("\nPartial R² increment for KF (full vs without KF): %.4f\n", partial_r2_kf))

vif_dag <- calculate_vif(model_adj)
cat("\nVIF (DAG adjusted model):\n")
print(round(vif_dag, 3))
cat("\n")

set.seed(123)
boot_kf <- function(data, indices) {
  d <- data[indices, , drop = FALSE]
  coef(lm(formula_adj, data = d))[2]
}
boot_out <- boot::boot(data = df_complete, statistic = boot_kf, R = 1000)
boot_ci <- boot::boot.ci(boot_out, type = c("perc", "bca"))
cat("Bootstrap 95% CI for KF coefficient (percentile):\n")
if (!is.null(boot_ci$percent)) {
  cat(sprintf("  [%.4f, %.4f]\n", boot_ci$percent[4], boot_ci$percent[5]))
}

cat("\n=== Analysis 15 complete ===\n")
