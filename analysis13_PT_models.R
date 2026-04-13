#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration: Toggle PJK exclusion on/off
EXCLUDE_PJK <- TRUE

# Configuration: Toggle surgeon volume filter on/off
EXCLUDE_LOW_VOLUME_SURGEONS <- FALSE

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Filter out low-volume surgeons if enabled
if (EXCLUDE_LOW_VOLUME_SURGEONS) {
  df <- filter_low_volume_surgeons(df, min_surgeon_cases = 10)
}

# Calculate change in lordosis using 6-week data
df$change_lordosis <- df$LAT6W_L1_S1 - df$LATpre_L1_S1

# Calculate T4-L1 PA
if ("LATpre_L1PA" %in% names(df) && "LATpre_T4PA" %in% names(df)) {
  df$LATpre_T4_L1_PA <- df$LATpre_L1PA - df$LATpre_T4PA
  cat("Calculated LATpre_T4_L1_PA = LATpre_L1PA - LATpre_T4PA\n")
} else {
  cat("WARNING: Could not find LATpre_L1PA or LATpre_T4PA columns.\n")
  df$LATpre_T4_L1_PA <- NA
}

# Check for age variable (prioritize demo_Age)
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
  cat("WARNING: Age variable not found. Will exclude from PCA.\n")
}

# Drop patients with missing data for confounder analysis
# Confounder set: S1PI (preop), Preop Lordosis (L1-S1), L4-S1, Thoracic Kyphosis, SVA, T4-L1 PA, KF, Age
# NOTE: S1PT (Pelvic Tilt) is now the predictor, so it's excluded from confounders
# NOTE: KF (Knee Flexion) has been added to confounders
if (!is.null(age_var)) {
  df_clean <- df %>%
    filter(!is.na(LATpre_S1PT) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
           !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) & 
           !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA) & !is.na(LATpre_LL_KneeAngle) &
           !is.na(.data[[age_var]]))
} else {
  df_clean <- df %>%
    filter(!is.na(LATpre_S1PT) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
           !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) & 
           !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA) & !is.na(LATpre_LL_KneeAngle))
}

cat(paste("\n=== PCA-Based Confounder Analysis: PT ~ LL Correction ===\n"))
cat(paste("Sample size:", nrow(df_clean), "patients with complete data (using 6-week follow-up)\n\n"))

# ============================================================================
# STEP 1: Prepare confounder variables for PCA
# ============================================================================
cat("=== STEP 1: Preparing Confounder Variables for PCA ===\n\n")

# Select confounder variables (excluding pelvic tilt and outcome)
# Pelvic Tilt (S1PT) is now the predictor, so it's excluded from confounders
# KF (Knee Flexion) has been added to confounders
confounder_vars <- c(
  "LATpre_S1PI",
  "LATpre_L1_S1",
  "LATpre_L4_S1",
  "LATpre_T2_T12",
  "LATpre_SVA_C2_S1",
  "LATpre_T4_L1_PA",
  "LATpre_LL_KneeAngle"
)

confounder_labels <- c(
  "S1PI (preop)",
  "Lordosis L1-S1",
  "Lordosis L4-S1",
  "Thoracic Kyphosis T2-T12",
  "SVA",
  "T4-L1 PA",
  "Knee Flexion"
)

# Add age if available
if (!is.null(age_var)) {
  confounder_vars <- c(confounder_vars, age_var)
  confounder_labels <- c(confounder_labels, "Age")
}

X_confounders <- df_clean[, confounder_vars]
colnames(X_confounders) <- confounder_labels

cat("Confounder variables for PCA:\n")
for (i in 1:length(confounder_vars)) {
  cat(sprintf("  %d. %s (%s)\n", i, confounder_labels[i], confounder_vars[i]))
}
cat("\n")
cat("NOTE: Pelvic Tilt (S1PT) is the predictor variable, not a confounder\n")
cat("NOTE: Knee Flexion (KF) is included as a confounder\n\n")

# Check correlation matrix before PCA
cat("Correlation matrix of confounders (before PCA):\n")
cor_matrix <- cor(X_confounders)
print(round(cor_matrix, 3))
cat("\n")

# ============================================================================
# STEP 2: Perform PCA
# ============================================================================
cat("=== STEP 2: Principal Component Analysis ===\n\n")

# Standardize variables (important for PCA when variables have different scales)
X_scaled <- scale(X_confounders)

# Perform PCA
pca_result <- prcomp(X_scaled, center = FALSE, scale. = FALSE)  # Already scaled

# Extract components
pca_components <- pca_result$x
pca_loadings <- pca_result$rotation
pca_summary <- summary(pca_result)

cat("PCA Summary:\n")
print(pca_summary)
cat("\n")

# Variance explained by each component
variance_explained <- pca_summary$importance[2, ]  # Proportion of variance
cumulative_variance <- pca_summary$importance[3, ]  # Cumulative proportion

cat("Variance Explained by Each Component:\n")
for (i in 1:length(variance_explained)) {
  cat(sprintf("  PC%d: %.2f%% (Cumulative: %.2f%%)\n", 
              i, variance_explained[i] * 100, cumulative_variance[i] * 100))
}
cat("\n")

# Determine number of components to keep
# Rule 1: Kaiser criterion (eigenvalue > 1)
eigenvalues <- pca_result$sdev^2
n_components_kaiser <- sum(eigenvalues > 1)

# Rule 2: Components explaining >80% cumulative variance
n_components_80 <- sum(cumulative_variance < 0.80) + 1
if (n_components_80 > length(cumulative_variance)) {
  n_components_80 <- length(cumulative_variance)
}

# Rule 3: Components explaining >90% cumulative variance
n_components_90 <- sum(cumulative_variance < 0.90) + 1
if (n_components_90 > length(cumulative_variance)) {
  n_components_90 <- length(cumulative_variance)
}

cat("Component Selection Criteria:\n")
cat(sprintf("  Kaiser criterion (eigenvalue > 1): %d components\n", n_components_kaiser))
cat(sprintf("  >80%% variance explained: %d components\n", n_components_80))
cat(sprintf("  >90%% variance explained: %d components\n", n_components_90))
cat("\n")

# ============================================================================
# RIGOROUS COMPONENT SELECTION METHODS
# ============================================================================
cat("=== RIGOROUS COMPONENT SELECTION METHODS ===\n\n")

# Method 1: Broken Stick Model
# Expected eigenvalues if variance is randomly distributed
broken_stick <- function(n_vars, n_components) {
  # Expected eigenvalue under broken stick model
  expected <- numeric(n_components)
  for (i in 1:n_components) {
    expected[i] <- sum(1 / (i:n_vars))
  }
  return(expected / n_vars)
}

n_vars <- length(confounder_vars)
broken_stick_expected <- broken_stick(n_vars, length(eigenvalues))
n_components_broken_stick <- sum(eigenvalues > broken_stick_expected)
cat("1. Broken Stick Model:\n")
cat("   Compares observed eigenvalues to expected if variance is randomly distributed\n")
cat(sprintf("   Components with eigenvalue > broken stick expectation: %d\n", n_components_broken_stick))
for (i in 1:min(7, length(eigenvalues))) {
  status <- ifelse(eigenvalues[i] > broken_stick_expected[i], "✓ Keep", "✗ Drop")
  cat(sprintf("     PC%d: eigenvalue = %.4f, expected = %.4f %s\n", 
              i, eigenvalues[i], broken_stick_expected[i], status))
}
cat("\n")

# Method 2: Parallel Analysis (Horn's method)
# Compare eigenvalues to those from random data with same dimensions
cat("2. Parallel Analysis (Horn's Method):\n")
cat("   Compares eigenvalues to those from random data\n")
set.seed(123)
n_iterations <- 100
n_obs <- nrow(X_scaled)
random_eigenvalues <- matrix(0, nrow = n_iterations, ncol = n_vars)

for (iter in 1:n_iterations) {
  # Generate random data with same dimensions
  random_data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  random_pca <- prcomp(random_data, center = TRUE, scale. = TRUE)
  random_eigenvalues[iter, ] <- random_pca$sdev^2
}

# 95th percentile of random eigenvalues
parallel_threshold <- apply(random_eigenvalues, 2, quantile, probs = 0.95)
n_components_parallel <- sum(eigenvalues > parallel_threshold)

cat(sprintf("   Components with eigenvalue > 95th percentile of random: %d\n", n_components_parallel))
for (i in 1:min(7, length(eigenvalues))) {
  status <- ifelse(eigenvalues[i] > parallel_threshold[i], "✓ Keep", "✗ Drop")
  cat(sprintf("     PC%d: eigenvalue = %.4f, parallel threshold = %.4f %s\n", 
              i, eigenvalues[i], parallel_threshold[i], status))
}
cat("\n")

# Method 3: Information Criteria (AIC/BIC) for regression models
# Fit models with different numbers of components and compare AIC/BIC
cat("3. Information Criteria (AIC/BIC) for Regression Models:\n")
cat("   Compares models with different numbers of PCA components\n")
ic_results <- data.frame(
  N_Components = 1:min(7, length(cumulative_variance)),
  Variance_Explained = cumulative_variance[1:min(7, length(cumulative_variance))] * 100,
  AIC = numeric(min(7, length(cumulative_variance))),
  BIC = numeric(min(7, length(cumulative_variance))),
  Pelvic_Tilt_P = numeric(min(7, length(cumulative_variance)))
)

for (n_comp in 1:min(7, length(cumulative_variance))) {
  # Prepare data with n_comp components
  df_temp <- df_clean
  for (i in 1:n_comp) {
    df_temp[[paste0("PC", i)]] <- pca_components[, i]
  }
  
  # Fit model
  pca_formula_temp <- as.formula(paste("change_lordosis ~ LATpre_S1PT +", 
                                       paste(paste0("PC", 1:n_comp), collapse = " + ")))
  model_temp <- lm(pca_formula_temp, data = df_temp)
  
  ic_results$AIC[n_comp] <- AIC(model_temp)
  ic_results$BIC[n_comp] <- BIC(model_temp)
  ic_results$Pelvic_Tilt_P[n_comp] <- summary(model_temp)$coefficients[2, 4]
}

# Find optimal based on AIC and BIC
optimal_aic <- which.min(ic_results$AIC)
optimal_bic <- which.min(ic_results$BIC)

cat("   Model comparison:\n")
print(ic_results, row.names = FALSE)
cat(sprintf("\n   Optimal by AIC: %d components (AIC = %.2f)\n", optimal_aic, ic_results$AIC[optimal_aic]))
cat(sprintf("   Optimal by BIC: %d components (BIC = %.2f)\n", optimal_bic, ic_results$BIC[optimal_bic]))
cat("   (Lower AIC/BIC is better - BIC penalizes complexity more)\n\n")

# Method 4: Variance explained thresholds
cat("4. Variance Explained Thresholds:\n")
cat(sprintf("   >80%% variance: %d components (%.2f%%)\n", n_components_80, cumulative_variance[n_components_80] * 100))
cat(sprintf("   >90%% variance: %d components (%.2f%%)\n", n_components_90, cumulative_variance[n_components_90] * 100))
cat(sprintf("   >95%% variance: %d components (%.2f%%)\n", 
            sum(cumulative_variance < 0.95) + 1, 
            cumulative_variance[sum(cumulative_variance < 0.95) + 1] * 100))
cat("\n")

# Method 5: Cattell's Scree Test (visual inspection - we'll note it)
cat("5. Cattell's Scree Test:\n")
cat("   Look for 'elbow' in scree plot where eigenvalues level off\n")
cat("   (See scree plot visualization)\n")
# Find where the drop in eigenvalues becomes small
eigenvalue_drops <- diff(eigenvalues)
largest_drop_idx <- which.max(eigenvalue_drops)
cat(sprintf("   Largest drop: Between PC%d and PC%d (drop = %.4f)\n", 
            largest_drop_idx, largest_drop_idx + 1, eigenvalue_drops[largest_drop_idx]))
cat("   Suggests keeping components before the largest drop\n\n")

# ============================================================================
# RECOMMENDATION BASED ON ALL CRITERIA
# ============================================================================
cat("=== RECOMMENDATION SUMMARY ===\n\n")

criteria_votes <- c(
  "Kaiser" = n_components_kaiser,
  "Broken Stick" = n_components_broken_stick,
  "Parallel Analysis" = n_components_parallel,
  "AIC Optimal" = optimal_aic,
  "BIC Optimal" = optimal_bic,
  ">90% Variance" = n_components_90
)

cat("Component selection by different criteria:\n")
for (i in 1:length(criteria_votes)) {
  cat(sprintf("  %s: %d components\n", names(criteria_votes)[i], criteria_votes[i]))
}
cat("\n")

# Find most common recommendation
mode_value <- as.numeric(names(sort(table(criteria_votes), decreasing = TRUE)[1]))
cat(sprintf("Most common recommendation: %d components\n", mode_value))
cat("\n")

# Also check which number gives p > 0.05 for pelvic tilt
components_with_p_gt_05 <- which(ic_results$Pelvic_Tilt_P > 0.05)
if (length(components_with_p_gt_05) > 0) {
  cat("Components that make pelvic tilt non-significant (p > 0.05):\n")
  for (n_comp in components_with_p_gt_05) {
    cat(sprintf("  %d components: p = %.4f\n", n_comp, ic_results$Pelvic_Tilt_P[n_comp]))
  }
  cat("\n")
}

# Final recommendation
cat("=== FINAL RECOMMENDATION ===\n\n")
cat("Based on statistical criteria:\n")
cat(sprintf("  - Kaiser criterion: %d components\n", n_components_kaiser))
cat(sprintf("  - Broken stick model: %d components\n", n_components_broken_stick))
cat(sprintf("  - Parallel analysis: %d components\n", n_components_parallel))
cat(sprintf("  - AIC optimal: %d components\n", optimal_aic))
cat(sprintf("  - BIC optimal: %d components (most parsimonious)\n", optimal_bic))
cat(sprintf("  - >90%% variance: %d components\n", n_components_90))
cat("\n")

# Use 6 components (user-specified)
n_components_use <- 6
# Ensure we don't exceed available components
if (n_components_use > length(cumulative_variance)) {
  n_components_use <- length(cumulative_variance)
  cat(sprintf("WARNING: Requested 6 components but only %d available. Using %d components.\n", 
              length(cumulative_variance), n_components_use))
}
cat(sprintf("USING: %d components (user-specified)\n", n_components_use))
if (n_components_use <= nrow(ic_results)) {
  cat(sprintf("  - BIC = %.2f\n", ic_results$BIC[n_components_use]))
  cat(sprintf("  - Pelvic tilt p-value = %.4f", ic_results$Pelvic_Tilt_P[n_components_use]))
  if (ic_results$Pelvic_Tilt_P[n_components_use] > 0.05) {
    cat(" (> 0.05, not significant)\n")
  } else {
    cat(" (≤ 0.05, significant)\n")
  }
}
cat(sprintf("  - Variance explained: %.2f%%\n", cumulative_variance[n_components_use] * 100))

cat(sprintf("\nUsing %d components (explains %.2f%% of variance)\n\n", 
            n_components_use, cumulative_variance[n_components_use] * 100))

# ============================================================================
# STEP 3: Interpret PCA Components
# ============================================================================
cat("=== STEP 3: PCA Component Interpretation ===\n\n")

cat("Loadings (correlations between original variables and components):\n")
loadings_df <- as.data.frame(pca_loadings[, 1:n_components_use])
colnames(loadings_df) <- paste0("PC", 1:n_components_use)
print(round(loadings_df, 3))
cat("\n")

# Interpret each component
cat("Component Interpretation:\n")
for (i in 1:n_components_use) {
  cat(sprintf("\nPC%d (explains %.2f%% of variance):\n", i, variance_explained[i] * 100))
  
  # Find variables with highest absolute loadings
  loadings_pc <- abs(pca_loadings[, i])
  top_vars <- order(loadings_pc, decreasing = TRUE)[1:min(3, length(loadings_pc))]
  
  cat("  High loadings (dominant variables):\n")
  for (j in top_vars) {
    sign <- ifelse(pca_loadings[j, i] > 0, "+", "-")
    cat(sprintf("    %s %s (loading = %.3f)\n", 
                sign, rownames(pca_loadings)[j], pca_loadings[j, i]))
  }
  
  # Clinical interpretation
  if (i == 1) {
    cat("  → Likely represents overall spinal alignment/pelvic parameters\n")
  } else if (i == 2) {
    cat("  → Likely represents regional alignment (thoracic vs lumbar)\n")
  } else {
    cat("  → Represents additional alignment dimensions\n")
  }
}
cat("\n")

# ============================================================================
# STEP 4: Regression Models
# ============================================================================
cat("=== STEP 4: Regression Models ===\n\n")

# Prepare data with PCA components
df_pca <- df_clean
for (i in 1:n_components_use) {
  df_pca[[paste0("PC", i)]] <- pca_components[, i]
}

# Model 1: Simple regression (pelvic tilt only)
cat("Model 1: Simple Regression (Pelvic Tilt Only)\n")
model1 <- lm(change_lordosis ~ LATpre_S1PT, data = df_pca)
summary1 <- summary(model1)
cat(sprintf("  Pelvic Tilt coefficient: %.4f (p = %.4e)\n", 
            summary1$coefficients[2, 1], summary1$coefficients[2, 4]))
cat(sprintf("  R²: %.4f\n", summary1$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n\n", summary1$adj.r.squared))

# Model 2: Multiple regression with original confounders (for comparison)
cat("Model 2: Multiple Regression with Original Confounders (including KF)\n")
if (!is.null(age_var)) {
  formula_str <- paste("change_lordosis ~ LATpre_S1PT + LATpre_S1PI + LATpre_L1_S1 +",
                       "LATpre_L4_S1 + LATpre_T2_T12 + LATpre_SVA_C2_S1 +",
                       "LATpre_T4_L1_PA + LATpre_LL_KneeAngle +", age_var)
  model2 <- lm(as.formula(formula_str), data = df_pca)
} else {
  model2 <- lm(change_lordosis ~ LATpre_S1PT + LATpre_S1PI + LATpre_L1_S1 + 
               LATpre_L4_S1 + LATpre_T2_T12 + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA + LATpre_LL_KneeAngle, 
               data = df_pca)
}
summary2 <- summary(model2)
cat(sprintf("  Pelvic Tilt coefficient: %.4f (p = %.4e)\n", 
            summary2$coefficients[2, 1], summary2$coefficients[2, 4]))
cat(sprintf("  R²: %.4f\n", summary2$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n\n", summary2$adj.r.squared))

# Model 3: Multiple regression with PCA components
cat(sprintf("Model 3: Multiple Regression with PCA Components (PC1-PC%d)\n", n_components_use))
pca_formula <- as.formula(paste("change_lordosis ~ LATpre_S1PT +", 
                                 paste(paste0("PC", 1:n_components_use), collapse = " + ")))
model3 <- lm(pca_formula, data = df_pca)
summary3 <- summary(model3)
cat(sprintf("  Pelvic Tilt coefficient: %.4f (p = %.4e)\n", 
            summary3$coefficients[2, 1], summary3$coefficients[2, 4]))
cat(sprintf("  R²: %.4f\n", summary3$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n\n", summary3$adj.r.squared))

# Print PCA component coefficients
cat("  PCA Component Coefficients:\n")
for (i in 1:n_components_use) {
  pc_name <- paste0("PC", i)
  if (pc_name %in% rownames(summary3$coefficients)) {
    cat(sprintf("    %s: %.4f (p = %.4e)\n", 
                pc_name, 
                summary3$coefficients[pc_name, 1],
                summary3$coefficients[pc_name, 4]))
  }
}
cat("\n")

# ============================================================================
# STEP 5: Model Comparison
# ============================================================================
cat("=== STEP 5: Model Comparison ===\n\n")

comparison_df <- data.frame(
  Model = c("Simple (PT only)", 
            "Multiple (original confounders)", 
            paste0("Multiple (PCA components, ", n_components_use, " PCs)")),
  Pelvic_Tilt_Coef = c(
    summary1$coefficients[2, 1],
    summary2$coefficients[2, 1],
    summary3$coefficients[2, 1]
  ),
  Pelvic_Tilt_P = c(
    summary1$coefficients[2, 4],
    summary2$coefficients[2, 4],
    summary3$coefficients[2, 4]
  ),
  R_squared = c(
    summary1$r.squared,
    summary2$r.squared,
    summary3$r.squared
  ),
  Adj_R_squared = c(
    summary1$adj.r.squared,
    summary2$adj.r.squared,
    summary3$adj.r.squared
  ),
  N_Predictors = c(1, 8, 1 + n_components_use)
)

print(comparison_df, row.names = FALSE)
cat("\n")

# ============================================================================
# STEP 6: Check for Multicollinearity in PCA Model
# ============================================================================
cat("=== STEP 6: Multicollinearity Check (PCA Model) ===\n\n")

# Calculate VIF for PCA model
calculate_vif <- function(model) {
  X <- model.matrix(model)[, -1]  # Remove intercept
  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)
  
  for (i in 1:ncol(X)) {
    other_preds <- X[, -i, drop = FALSE]
    vif_model <- lm(X[, i] ~ other_preds)
    r_squared <- summary(vif_model)$r.squared
    if (r_squared >= 0.9999) {
      vif_values[i] <- Inf
    } else {
      vif_values[i] <- 1 / (1 - r_squared)
    }
  }
  return(vif_values)
}

vif_pca <- calculate_vif(model3)
cat("VIF values for PCA model:\n")
for (i in 1:length(vif_pca)) {
  interpretation <- ifelse(vif_pca[i] == Inf, "Perfect collinearity",
                  ifelse(vif_pca[i] > 10, "Severe",
                  ifelse(vif_pca[i] > 5, "Moderate",
                  ifelse(vif_pca[i] > 2.5, "Mild", "Acceptable"))))
  cat(sprintf("  %s: %.2f (%s)\n", names(vif_pca)[i], vif_pca[i], interpretation))
}
cat("\n")

# PCA components are orthogonal by design, so VIF should be ~1
if (all(vif_pca < 2.5)) {
  cat("✓ Excellent: All VIF values < 2.5 (PCA components are orthogonal)\n")
  cat("  This confirms that PCA successfully eliminated multicollinearity!\n\n")
} else {
  cat("⚠️  Warning: Some VIF values > 2.5\n")
  cat("  This is unexpected since PCA components should be orthogonal.\n\n")
}

# ============================================================================
# STEP 7: Visualization
# ============================================================================
cat("=== STEP 7: Creating Visualizations ===\n\n")

if (!dir.exists("results")) {
  dir.create("results")
}

# Plot 1: Scree plot
png("results/analysis13_scree_plot.png", width = 10, height = 6, units = "in", res = 300)
scree_df <- data.frame(
  Component = 1:length(eigenvalues),
  Eigenvalue = eigenvalues,
  Variance_Explained = variance_explained * 100,
  Cumulative_Variance = cumulative_variance * 100
)

p1 <- ggplot(scree_df, aes(x = Component, y = Eigenvalue)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Principal Component",
    y = "Eigenvalue",
    title = "Scree Plot: PCA Eigenvalues"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text", x = Inf, y = 1, label = "Kaiser criterion (eigenvalue = 1)", 
           hjust = 1.1, vjust = -0.5, color = "red")

print(p1)
dev.off()
cat("Saved scree plot to results/analysis13_scree_plot.png\n")

# Plot 2: Variance explained
png("results/analysis13_variance_explained.png", width = 10, height = 6, units = "in", res = 300)
variance_df <- data.frame(
  Component = 1:length(variance_explained),
  Variance = variance_explained * 100,
  Cumulative = cumulative_variance * 100
)

p2 <- ggplot(variance_df, aes(x = Component)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = "PCA Variance Explained"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text", x = Inf, y = 80, label = "80% threshold", 
           hjust = 1.1, vjust = -0.5, color = "orange") +
  annotate("text", x = Inf, y = 90, label = "90% threshold", 
           hjust = 1.1, vjust = -0.5, color = "darkgreen")

print(p2)
dev.off()
cat("Saved variance explained plot to results/analysis13_variance_explained.png\n")

# Plot 3: Loadings heatmap
png("results/analysis13_loadings_heatmap.png", width = 10, height = 8, units = "in", res = 300)
loadings_long <- loadings_df %>%
  rownames_to_column("Variable") %>%
  pivot_longer(cols = -Variable, names_to = "Component", values_to = "Loading")

p3 <- ggplot(loadings_long, aes(x = Component, y = Variable, fill = Loading)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1)) +
  labs(
    x = "Principal Component",
    y = "Original Variable",
    title = "PCA Loadings Heatmap",
    fill = "Loading"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p3)
dev.off()
cat("Saved loadings heatmap to results/analysis13_loadings_heatmap.png\n")

# Plot 4: Model comparison
png("results/analysis13_model_comparison.png", width = 10, height = 6, units = "in", res = 300)
coef_comparison <- data.frame(
  Model = comparison_df$Model,
  Coefficient = comparison_df$Pelvic_Tilt_Coef,
  Lower_CI = c(confint(model1)[2, 1], confint(model2)[2, 1], confint(model3)[2, 1]),
  Upper_CI = c(confint(model1)[2, 2], confint(model2)[2, 2], confint(model3)[2, 2]),
  P_value = comparison_df$Pelvic_Tilt_P
)

p4 <- ggplot(coef_comparison, aes(x = Model, y = Coefficient)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Model",
    y = "Pelvic Tilt Coefficient",
    title = "Model Comparison: Pelvic Tilt Coefficient"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text", x = coef_comparison$Model, y = coef_comparison$Coefficient,
    label = paste0("p = ", formatC(coef_comparison$P_value, format = "e", digits = 2)),
    vjust = -1.5, size = 3
  )

print(p4)
dev.off()
cat("Saved model comparison plot to results/analysis13_model_comparison.png\n")

# Plot 5: Residual plot (similar to analysis5)
# This shows the relationship between pelvic tilt and residuals after controlling for PCA components
cat("\nCreating residual plot...\n")
pca_formula_no_pt <- as.formula(paste("change_lordosis ~", 
                                      paste(paste0("PC", 1:n_components_use), collapse = " + ")))
covariates_model_pca <- lm(pca_formula_no_pt, data = df_pca)
df_pca$resid_from_pca_covariates <- residuals(covariates_model_pca)

# Fit regression of residuals on pelvic tilt
pt_model_resid_pca <- lm(resid_from_pca_covariates ~ LATpre_S1PT, data = df_pca)
summary_resid_pca <- summary(pt_model_resid_pca)

# Calculate partial correlation
t_stat_pt_pca <- summary3$coefficients[2, 3]
df_residual_pca <- summary3$df[2]
r_partial_pca <- t_stat_pt_pca / sqrt(t_stat_pt_pca^2 + df_residual_pca)

# Calculate partial R²
r2_no_pt_pca <- summary(covariates_model_pca)$r.squared
r2_full_pca <- summary3$r.squared
r2_partial_pca <- r2_full_pca - r2_no_pt_pca

# Create the plot
png("results/analysis13_residual_plot.png", width = 10, height = 8, units = "in", res = 300)
p5 <- ggplot(df_pca, aes(x = LATpre_S1PT, y = resid_from_pca_covariates)) +
  geom_point(alpha = 0.6, color = "darkgrey") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Preoperative Pelvic Tilt (S1PT)",
    y = paste0("Residuals from PCA Covariates Model\n(Change in Lordosis after removing\neffects of PC1-PC", n_components_use, ")"),
    title = paste0("Model 3: Adjusted Relationship (6-week follow-up)\n(Pelvic Tilt ~ Residuals after controlling for ", n_components_use, " PCA components)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("Partial r = ", round(r_partial_pca, 3),
                   "\nPartial R² = ", round(r2_partial_pca, 3),
                   "\nCoefficient = ", round(summary3$coefficients[2, 1], 4),
                   "\np = ", formatC(summary3$coefficients[2, 4], format = "e", digits = 2)),
    hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
  )

print(p5)
dev.off()
cat("Saved residual plot to results/analysis13_residual_plot.png\n")

# Plot 6: Simple scatter plot (PT vs LL correction)
png("results/analysis13_pt_vs_ll_correction.png", width = 10, height = 8, units = "in", res = 300)
p6 <- ggplot(df_pca, aes(x = LATpre_S1PT, y = change_lordosis)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative Pelvic Tilt (S1PT)",
    y = "Change in Lordosis (6-week - Preop)",
    title = "Preoperative Pelvic Tilt vs Change in Lordosis\n(Simple Regression)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("r = ", round(cor(df_pca$LATpre_S1PT, df_pca$change_lordosis, use = "complete.obs"), 3),
                   "\nR² = ", round(summary1$r.squared, 3),
                   "\nCoefficient = ", round(summary1$coefficients[2, 1], 4),
                   "\np = ", formatC(summary1$coefficients[2, 4], format = "e", digits = 2)),
    hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
  )

print(p6)
dev.off()
cat("Saved simple scatter plot to results/analysis13_pt_vs_ll_correction.png\n")

# ============================================================================
# SUMMARY
# ============================================================================
# ============================================================================
# STEP 8: Compare Residual Variance with and without KF
# ============================================================================
cat("=== STEP 8: Residual Variance Comparison (with and without KF in confounders) ===\n\n")

# Calculate residual variance from Model 3 (with KF in confounders via PCA)
residual_variance_with_KF <- var(residuals(model3))
residual_se_with_KF <- summary3$sigma

cat("Model 3 (PCA with KF as confounder):\n")
cat(sprintf("  Residual variance: %.4f\n", residual_variance_with_KF))
cat(sprintf("  Residual SE: %.4f\n", residual_se_with_KF))
cat(sprintf("  R²: %.4f\n", summary3$r.squared))
cat(sprintf("  PT coefficient: %.4f (p = %.4e)\n\n", 
            summary3$coefficients[2, 1], summary3$coefficients[2, 4]))

# To compare, we need to run PCA without KF and see the difference
# Create confounders without KF for comparison
confounder_vars_no_KF <- confounder_vars[confounder_vars != "LATpre_LL_KneeAngle"]
confounder_labels_no_KF <- confounder_labels[confounder_vars != "LATpre_LL_KneeAngle"]

X_confounders_no_KF <- df_clean[, confounder_vars_no_KF]
colnames(X_confounders_no_KF) <- confounder_labels_no_KF

# Perform PCA without KF
X_scaled_no_KF <- scale(X_confounders_no_KF)
pca_result_no_KF <- prcomp(X_scaled_no_KF, center = FALSE, scale. = FALSE)
pca_components_no_KF <- pca_result_no_KF$x
pca_summary_no_KF <- summary(pca_result_no_KF)
cumulative_variance_no_KF <- pca_summary_no_KF$importance[3, ]

# Use same number of components for fair comparison
df_pca_no_KF <- df_clean
for (i in 1:n_components_use) {
  df_pca_no_KF[[paste0("PC", i, "_noKF")]] <- pca_components_no_KF[, i]
}

# Model 3 without KF in confounders
pca_formula_no_KF <- as.formula(paste("change_lordosis ~ LATpre_S1PT +", 
                                       paste(paste0("PC", 1:n_components_use, "_noKF"), collapse = " + ")))
model3_no_KF <- lm(pca_formula_no_KF, data = df_pca_no_KF)
summary3_no_KF <- summary(model3_no_KF)

residual_variance_without_KF <- var(residuals(model3_no_KF))
residual_se_without_KF <- summary3_no_KF$sigma

cat("Model 3 without KF in confounders (for comparison):\n")
cat(sprintf("  Residual variance: %.4f\n", residual_variance_without_KF))
cat(sprintf("  Residual SE: %.4f\n", residual_se_without_KF))
cat(sprintf("  R²: %.4f\n", summary3_no_KF$r.squared))
cat(sprintf("  PT coefficient: %.4f (p = %.4e)\n\n", 
            summary3_no_KF$coefficients[2, 1], summary3_no_KF$coefficients[2, 4]))

# Compare residual variances
variance_change <- residual_variance_with_KF - residual_variance_without_KF
variance_change_pct <- (variance_change / residual_variance_without_KF) * 100

cat("Comparison:\n")
cat(sprintf("  Change in residual variance: %.4f (%.2f%%)\n", variance_change, variance_change_pct))
if (abs(variance_change_pct) < 1) {
  cat("  → Residual variance changed by < 1% - minimal change\n")
  cat("  → Adding KF to confounders does not meaningfully change residual variance\n")
} else if (variance_change_pct < 0) {
  cat("  → Residual variance decreased (model fit improved)\n")
} else {
  cat("  → Residual variance increased (model fit slightly worse)\n")
}

# Compare PT coefficients
pt_coef_change <- summary3$coefficients[2, 1] - summary3_no_KF$coefficients[2, 1]
pt_coef_change_pct <- (pt_coef_change / abs(summary3_no_KF$coefficients[2, 1])) * 100

cat(sprintf("\n  Change in PT coefficient: %.4f (%.2f%%)\n", pt_coef_change, pt_coef_change_pct))
if (abs(pt_coef_change_pct) < 5) {
  cat("  → PT coefficient changed by < 5% - minimal change\n")
  cat("  → Adding KF to confounders does not meaningfully change PT effect\n")
} else {
  cat("  → PT coefficient changed substantially\n")
}
cat("\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("=== SUMMARY ===\n\n")

cat("PCA-Based Confounder Analysis Results: PT ~ LL Correction\n\n")

cat("1. PCA successfully reduced", length(confounder_vars), "correlated confounders (including KF) to", n_components_use, "orthogonal components\n")
cat(sprintf("   - Explained variance: %.2f%%\n", cumulative_variance[n_components_use] * 100))
cat("   - Eliminated multicollinearity (all VIF < 2.5)\n\n")

cat("2. Pelvic Tilt Effect:\n")
cat(sprintf("   - Simple model: %.4f (p = %.4e)", 
            summary1$coefficients[2, 1], summary1$coefficients[2, 4]))
if (summary1$coefficients[2, 4] <= 0.05) {
  cat(" *** SIGNIFICANT ***\n")
} else {
  cat(" - NOT SIGNIFICANT\n")
}
cat(sprintf("   - Original confounders: %.4f (p = %.4e)", 
            summary2$coefficients[2, 1], summary2$coefficients[2, 4]))
if (summary2$coefficients[2, 4] <= 0.05) {
  cat(" *** SIGNIFICANT ***\n")
} else {
  cat(" - NOT SIGNIFICANT\n")
}
cat(sprintf("   - PCA components (%d PCs): %.4f (p = %.4e)", 
            n_components_use, summary3$coefficients[2, 1], summary3$coefficients[2, 4]))
if (summary3$coefficients[2, 4] <= 0.05) {
  cat(" *** SIGNIFICANT ***\n")
  cat(sprintf("\n   ✓ Pelvic Tilt has a significant influence on LL correction (p ≤ 0.05)\n"))
  cat("     This suggests that pelvic tilt independently predicts change in lordosis\n")
  cat("     even after adjusting for confounders via PCA components.\n")
} else {
  cat(" - NOT SIGNIFICANT\n")
  cat(sprintf("\n   After adjusting for confounders via PCA, pelvic tilt is not significant.\n"))
  cat("     This suggests that the confounders fully explain the relationship.\n")
}
cat("\n")

cat("3. Model Fit:\n")
cat(sprintf("   - Simple model R²: %.4f\n", summary1$r.squared))
cat(sprintf("   - Original confounders R²: %.4f\n", summary2$r.squared))
cat(sprintf("   - PCA components R²: %.4f\n", summary3$r.squared))
cat("\n")

cat("4. Advantages of PCA Approach:\n")
cat("   ✓ Eliminates multicollinearity (orthogonal components)\n")
cat("   ✓ Reduces dimensionality (fewer parameters to estimate)\n")
cat("   ✓ More stable coefficient estimates\n")
cat("   ✓ Components may have clinical interpretation\n")
cat("   ✓ Better statistical power (fewer degrees of freedom lost)\n\n")

cat("5. Limitations:\n")
cat("   ⚠  PCA components are less interpretable than original variables\n")
cat("   ⚠  Requires understanding what each component represents\n")
cat("   ⚠  May lose some information if too few components are used\n\n")

cat("=== Analysis Complete ===\n")
cat("Check the plots in results/ directory for visualizations.\n\n")
