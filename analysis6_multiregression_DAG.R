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
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Filter out low-volume surgeons if enabled
if (EXCLUDE_LOW_VOLUME_SURGEONS) {
  df <- filter_low_volume_surgeons(df, min_surgeon_cases = 10)
}

# Calculate change in lordosis using 6-week data
df$change_lordosis <- df$LAT6W_L1_S1 - df$LATpre_L1_S1

# Drop patients with missing data for DAG-based confounder analysis
# Covariates based on DAG: S1PI (preop), TK (T2-T12), S1PT, SVA, Age (PI-LL, L1-S1 Lordosis, L1PA removed)
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
  cat("WARNING: Age variable not found. Will exclude from model.\n")
}

if (!is.null(age_var)) {
  df_clean <- df %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
           !is.na(LATpre_T2_T12) & !is.na(LATpre_S1PT) & 
           !is.na(LATpre_SVA_C2_S1) & !is.na(.data[[age_var]]))
} else {
  df_clean <- df %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
           !is.na(LATpre_T2_T12) & !is.na(LATpre_S1PT) & 
           !is.na(LATpre_SVA_C2_S1))
}

cat(paste("\n=== DAG-Based Multiple Regression Analysis ===\n"))
cat(paste("Sample size:", nrow(df_clean), "patients with complete data (using 6-week follow-up)\n\n"))

cat("Variables included in model (based on DAG):\n")
cat("  - Preop Knee Flexion (LATpre_LL_KneeAngle) - Main predictor\n")
cat("  - Preop S1PI (LATpre_S1PI)\n")
cat("  - Preop Thoracic Kyphosis T2-T12 (LATpre_T2_T12)\n")
cat("  - Preop S1PT - Pelvic Tilt (LATpre_S1PT)\n")
cat("  - Preop SVA (LATpre_SVA_C2_S1)\n")
if (!is.null(age_var)) {
  cat(sprintf("  - Age (%s)\n", age_var))
}
cat("  - Change in Lordosis (LAT6W_L1_S1 - LATpre_L1_S1, 6-week follow-up) - Outcome\n")
cat("  (Note: PI-LL, L1-S1 Lordosis, and L1PA removed)\n\n")

# ============================================================================
# MODEL 1: Simple regression (knee flexion only)
# ============================================================================
cat("=== Model 1: Simple Regression (Knee Flexion Only) ===\n")
model1 <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = df_clean)
summary1 <- summary(model1)
cat("Coefficient for Preop Knee Flexion:\n")
cat(sprintf("  Estimate: %.4f\n", summary1$coefficients[2, 1]))
cat(sprintf("  P-value: %.4e\n", summary1$coefficients[2, 4]))
cat(sprintf("  R²: %.4f\n", summary1$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n\n", summary1$adj.r.squared))

# ============================================================================
# MODEL 2: Multiple regression (DAG-based covariates)
# ============================================================================
cat("=== Model 2: Multiple Regression (DAG-Based Covariates) ===\n")
if (!is.null(age_var)) {
  cat("Adjusting for: S1PI (preop) + Thoracic Kyphosis (T2-T12) + S1PT + SVA + Age\n")
  cat("(PI-LL, L1-S1 Lordosis, and L1PA removed)\n\n")
  formula_str <- paste("change_lordosis ~ LATpre_LL_KneeAngle + LATpre_S1PI +", 
                       "LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 +", age_var)
  model2 <- lm(as.formula(formula_str), data = df_clean)
} else {
  cat("Adjusting for: S1PI (preop) + Thoracic Kyphosis (T2-T12) + S1PT + SVA\n")
  cat("(PI-LL, L1-S1 Lordosis, L1PA, and Age removed - age variable not found)\n\n")
  model2 <- lm(change_lordosis ~ LATpre_LL_KneeAngle + LATpre_S1PI + 
               LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1, data = df_clean)
}
summary2 <- summary(model2)

cat("Coefficient for Preop Knee Flexion (adjusted for all covariates):\n")
cat(sprintf("  Estimate: %.4f\n", summary2$coefficients[2, 1]))
cat(sprintf("  P-value: %.4e\n", summary2$coefficients[2, 4]))
cat(sprintf("  R²: %.4f\n", summary2$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n\n", summary2$adj.r.squared))

# Print all coefficients
cat("All Coefficients:\n")
coef_names <- rownames(summary2$coefficients)
for (i in 1:nrow(summary2$coefficients)) {
  var_name <- coef_names[i]
  if (var_name == "(Intercept)") {
    cat(sprintf("  %s: %.4f (p = %.4e)\n", var_name, 
                summary2$coefficients[i, 1], summary2$coefficients[i, 4]))
  } else {
    # Format variable name for readability
    var_label <- switch(var_name,
                       "LATpre_LL_KneeAngle" = "Knee Flexion",
                       "LATpre_S1PI" = "S1PI (preop)",
                       "LATpre_T2_T12" = "Thoracic Kyphosis (T2-T12)",
                       "LATpre_S1PT" = "S1PT",
                       "LATpre_SVA_C2_S1" = "SVA",
                       ifelse(!is.null(age_var) && var_name == age_var, "Age", var_name))
    cat(sprintf("  %s: %.4f (p = %.4e)\n", var_label, 
                summary2$coefficients[i, 1], summary2$coefficients[i, 4]))
  }
}
cat("\n")

# ============================================================================
# VIF ANALYSIS
# ============================================================================
cat("=== VIF Analysis ===\n\n")

# Calculate VIF for each variable
calculate_vif <- function(model) {
  X <- model.matrix(model)[, -1]  # Remove intercept
  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)
  
  for (i in 1:ncol(X)) {
    other_preds <- X[, -i, drop = FALSE]
    vif_model <- lm(X[, i] ~ other_preds)
    r_squared <- summary(vif_model)$r.squared
    
    # Handle perfect collinearity (R² = 1)
    if (r_squared >= 0.9999) {
      vif_values[i] <- Inf
    } else {
      vif_values[i] <- 1 / (1 - r_squared)
    }
  }
  
  return(vif_values)
}

vif_values <- calculate_vif(model2)

# Create VIF data frame with interpretation
get_vif_interpretation <- function(vif) {
  if (is.infinite(vif)) {
    return("Perfect collinearity (variable should be removed)")
  } else if (vif > 10) {
    return("Severe multicollinearity")
  } else if (vif > 5) {
    return("Moderate multicollinearity")
  } else if (vif > 2.5) {
    return("Mild multicollinearity")
  } else {
    return("Acceptable")
  }
}

vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = vif_values,
  Interpretation = sapply(vif_values, get_vif_interpretation)
)

# Format variable names for readability
vif_df$Variable_Label <- sapply(vif_df$Variable, function(x) {
  switch(x,
         "LATpre_LL_KneeAngle" = "Knee Flexion",
         "LATpre_S1PI" = "S1PI (preop)",
         "LATpre_T2_T12" = "Thoracic Kyphosis (T2-T12)",
         "LATpre_S1PT" = "S1PT",
         "LATpre_SVA_C2_S1" = "SVA",
         ifelse(!is.null(age_var) && x == age_var, "Age", x))
})

vif_df <- vif_df[order(vif_df$VIF, decreasing = TRUE), ]

cat("VIF values (higher = more collinear):\n")
print(vif_df[, c("Variable_Label", "VIF", "Interpretation")], row.names = FALSE)
cat("\n")

# Interpretation guidelines
cat("VIF Interpretation Guidelines:\n")
cat("  VIF < 2.5: Acceptable (minimal multicollinearity)\n")
cat("  2.5 ≤ VIF < 5: Mild multicollinearity (may be acceptable)\n")
cat("  5 ≤ VIF < 10: Moderate multicollinearity (concerning)\n")
cat("  VIF ≥ 10: Severe multicollinearity (problematic)\n")
cat("  VIF = Inf: Perfect collinearity (variable must be removed)\n\n")

# Count problematic variables
n_severe <- sum(vif_values > 10 & vif_values != Inf)
n_moderate <- sum(vif_values > 5 & vif_values <= 10)
n_mild <- sum(vif_values > 2.5 & vif_values <= 5)
n_acceptable <- sum(vif_values <= 2.5)

cat("Summary:\n")
cat(sprintf("  Acceptable (VIF < 2.5): %d variable(s)\n", n_acceptable))
cat(sprintf("  Mild (2.5 ≤ VIF < 5): %d variable(s)\n", n_mild))
cat(sprintf("  Moderate (5 ≤ VIF < 10): %d variable(s)\n", n_moderate))
cat(sprintf("  Severe (VIF ≥ 10): %d variable(s)\n", n_severe + sum(is.infinite(vif_values))))

if (n_severe > 0 || any(is.infinite(vif_values))) {
  cat("\n⚠️  WARNING: Found variable(s) with severe/perfect multicollinearity\n")
  cat("   These variables should be removed or the model should be reconsidered.\n\n")
} else if (n_moderate > 0) {
  cat("\n⚠️  CAUTION: Found variable(s) with moderate multicollinearity\n")
  cat("   Consider removing or combining these variables.\n\n")
} else if (n_mild > 0) {
  cat("\n✓ Mild multicollinearity detected. This is generally acceptable.\n\n")
} else {
  cat("\n✓ No multicollinearity issues detected. All VIF values are acceptable.\n\n")
}

# ============================================================================
# CORRELATION MATRIX
# ============================================================================
cat("=== Correlation Matrix ===\n\n")

# Select predictor variables
if (!is.null(age_var)) {
  predictors <- df_clean %>%
    select(
      LATpre_LL_KneeAngle,
      LATpre_S1PI,
      LATpre_T2_T12,
      LATpre_S1PT,
      LATpre_SVA_C2_S1,
      all_of(age_var)
    )
  cor_labels <- c("Knee Flexion", "S1PI (preop)", "T2-T12", "S1PT", "SVA", "Age")
} else {
  predictors <- df_clean %>%
    select(
      LATpre_LL_KneeAngle,
      LATpre_S1PI,
      LATpre_T2_T12,
      LATpre_S1PT,
      LATpre_SVA_C2_S1
    )
  cor_labels <- c("Knee Flexion", "S1PI (preop)", "T2-T12", "S1PT", "SVA")
}

# Create correlation matrix
cor_matrix <- cor(predictors, use = "complete.obs")
colnames(cor_matrix) <- cor_labels
rownames(cor_matrix) <- colnames(cor_matrix)

cat("Correlation matrix (rounded to 3 decimals):\n")
print(round(cor_matrix, 3))
cat("\n")

# Identify highly correlated pairs (|r| > 0.7)
cat("Highly Correlated Pairs (|r| > 0.7):\n")
high_cor_pairs <- data.frame()
for (i in 1:(nrow(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    r_val <- cor_matrix[i, j]
    if (abs(r_val) > 0.7) {
      high_cor_pairs <- rbind(high_cor_pairs, 
                             data.frame(
                               Var1 = rownames(cor_matrix)[i],
                               Var2 = colnames(cor_matrix)[j],
                               Correlation = r_val
                             ))
    }
  }
}

if (nrow(high_cor_pairs) > 0) {
  high_cor_pairs <- high_cor_pairs[order(abs(high_cor_pairs$Correlation), decreasing = TRUE), ]
  print(high_cor_pairs)
  cat("\n⚠️  WARNING: Found", nrow(high_cor_pairs), "highly correlated pairs (|r| > 0.7)\n")
  cat("   This suggests potential multicollinearity issues.\n\n")
} else {
  cat("No highly correlated pairs found (|r| > 0.7).\n\n")
}

# ============================================================================
# MODEL COMPARISON
# ============================================================================
cat("=== Model Comparison ===\n\n")

comparison_df <- data.frame(
  Model = c("Simple (KF only)", "Multiple (DAG-based covariates)"),
  Knee_Flexion_Coef = c(
    summary1$coefficients[2, 1],
    summary2$coefficients[2, 1]
  ),
  Knee_Flexion_P = c(
    summary1$coefficients[2, 4],
    summary2$coefficients[2, 4]
  ),
  R_squared = c(
    summary1$r.squared,
    summary2$r.squared
  ),
  Adj_R_squared = c(
    summary1$adj.r.squared,
    summary2$adj.r.squared
  ),
  Max_VIF = c(
    NA,
    max(vif_values[names(vif_values) != "LATpre_LL_KneeAngle"])
  )
)

print(comparison_df, row.names = FALSE)
cat("\n")

# ============================================================================
# VISUALIZATION
# ============================================================================
cat("=== Creating Visualizations ===\n\n")

if (!dir.exists("results")) {
  dir.create("results")
}

# Plot 1: VIF bar plot
png("results/analysis6_vif_plot.png", width = 10, height = 6, units = "in", res = 300)
vif_plot_df <- vif_df
vif_plot_df$VIF_plot <- ifelse(is.infinite(vif_plot_df$VIF), 50, vif_plot_df$VIF)  # Cap at 50 for plotting
vif_plot_df$Variable_Label <- factor(vif_plot_df$Variable_Label, 
                                     levels = vif_plot_df$Variable_Label[order(vif_plot_df$VIF_plot)])

p_vif <- ggplot(vif_plot_df, aes(x = Variable_Label, y = VIF_plot, fill = Interpretation)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("Acceptable" = "green", 
                               "Mild multicollinearity" = "yellow",
                               "Moderate multicollinearity" = "orange",
                               "Severe multicollinearity" = "red",
                               "Perfect collinearity (variable should be removed)" = "darkred")) +
  labs(
    x = "Variable",
    y = "Variance Inflation Factor (VIF)",
    title = "Multicollinearity Diagnostic: VIF Values\n(DAG-Based Model)",
    fill = "Interpretation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text", x = Inf, y = 5, label = "VIF = 5", hjust = 1.1, vjust = -0.5, color = "orange") +
  annotate("text", x = Inf, y = 10, label = "VIF = 10", hjust = 1.1, vjust = -0.5, color = "red")

ggsave("results/analysis6_vif_plot.png", plot = p_vif, width = 10, height = 6, dpi = 300)
cat("Saved VIF plot to results/analysis6_vif_plot.png\n")

# Plot 2: Correlation heatmap
png("results/analysis6_correlation_heatmap.png", width = 10, height = 8, units = "in", res = 300)
if (requireNamespace("corrplot", quietly = TRUE)) {
  library(corrplot)
  corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
           tl.cex = 0.8, tl.col = "black", addCoef.col = "black",
           number.cex = 0.7, col = colorRampPalette(c("blue", "white", "red"))(200),
           title = "Correlation Matrix of Predictor Variables\n(DAG-Based Model)", mar = c(0, 0, 2, 0))
  dev.off()
  cat("Saved correlation heatmap to results/analysis6_correlation_heatmap.png\n")
} else {
  dev.off()
  cat("corrplot package not available. Skipping correlation heatmap.\n")
}

# Plot 3: Model comparison
png("results/analysis6_coefficient_comparison.png", width = 10, height = 6, units = "in", res = 300)
coef_comparison <- data.frame(
  Model = comparison_df$Model,
  Coefficient = comparison_df$Knee_Flexion_Coef,
  Lower_CI = c(confint(model1)[2, 1], confint(model2)[2, 1]),
  Upper_CI = c(confint(model1)[2, 2], confint(model2)[2, 2]),
  P_value = comparison_df$Knee_Flexion_P
)

p3 <- ggplot(coef_comparison, aes(x = Model, y = Coefficient)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Model",
    y = "Knee Flexion Coefficient",
    title = "Model Comparison: Knee Flexion Coefficient\n(DAG-Based Model)"
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

ggsave("results/analysis6_coefficient_comparison.png", plot = p3, width = 10, height = 6, dpi = 300)
cat("Saved coefficient comparison plot to results/analysis6_coefficient_comparison.png\n")

# Plot 4: Unadjusted relationship (same style as analysis4)
cat("\nCreating regression line plots (analysis4 style)...\n")
r_kf_lordosis <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$change_lordosis, use = "complete.obs")

p4 <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = change_lordosis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    x = "Preoperative Knee Flexion",
    y = "Change in Lordosis (6-week)",
    title = "Model 1: Unadjusted\n(Preop Knee Flexion Only, 6-week follow-up)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("r = ", round(r_kf_lordosis, 3), 
                   "\nR² = ", round(summary1$r.squared, 3),
                   "\np = ", formatC(summary1$coefficients[2, 4], format = "e", digits = 2)),
    hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold"
  )

ggsave("results/analysis6_unadjusted.png", plot = p4, width = 8, height = 6, dpi = 300)
cat("Saved unadjusted plot to results/analysis6_unadjusted.png\n")

# Plot 5: Adjusted relationship (residuals from covariates model vs knee flexion)
# Build model without knee flexion to get residuals
if (!is.null(age_var)) {
  covariates_model <- lm(as.formula(paste("change_lordosis ~ LATpre_S1PI +", 
                                          "LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 +", age_var)), 
                        data = df_clean)
} else {
  covariates_model <- lm(change_lordosis ~ LATpre_S1PI + LATpre_T2_T12 + 
                        LATpre_S1PT + LATpre_SVA_C2_S1, data = df_clean)
}
df_clean$resid_from_covariates <- residuals(covariates_model)

kf_model_resid <- lm(resid_from_covariates ~ LATpre_LL_KneeAngle, data = df_clean)
summary_resid <- summary(kf_model_resid)

# Calculate partial correlation
t_stat_kf <- summary2$coefficients[2, 3]  # t-statistic for knee flexion
df_residual <- summary2$df[2]  # residual degrees of freedom
r_partial <- t_stat_kf / sqrt(t_stat_kf^2 + df_residual)

# Calculate partial R²
summary_no_kf <- summary(covariates_model)
r2_no_kf <- summary_no_kf$r.squared
r2_full <- summary2$r.squared
r2_partial <- r2_full - r2_no_kf  # Partial R² = R² change

# Create covariate list for label
if (!is.null(age_var)) {
  covar_list <- "S1PI, Thoracic Kyphosis,\nS1PT, SVA, Age"
} else {
  covar_list <- "S1PI, Thoracic Kyphosis,\nS1PT, SVA"
}

p5 <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = resid_from_covariates)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Preoperative Knee Flexion",
    y = paste0("Residuals from Covariates Model\n(Change in Lordosis after removing\neffects of ", covar_list, ")"),
    title = "Model 2: Adjusted Relationship (6-week follow-up)\n(Knee Flexion ~ Residuals after controlling for all covariates)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("Partial r = ", round(r_partial, 3),
                   "\nPartial R² = ", round(r2_partial, 3),
                   "\nCoefficient = ", round(summary2$coefficients[2, 1], 4),
                   "\np = ", formatC(summary2$coefficients[2, 4], format = "e", digits = 2)),
    hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold"
  )

ggsave("results/analysis6_adjusted.png", plot = p5, width = 8, height = 6, dpi = 300)
cat("Saved adjusted (residual) plot to results/analysis6_adjusted.png\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n=== SUMMARY ===\n\n")

cat("DAG-Based Multiple Regression Results:\n\n")

cat("1. Model Specification:\n")
cat("   - Main predictor: Knee Flexion\n")
if (!is.null(age_var)) {
  cat("   - Covariates: S1PI (preop), Thoracic Kyphosis (T2-T12), S1PT, SVA, Age\n")
} else {
  cat("   - Covariates: S1PI (preop), Thoracic Kyphosis (T2-T12), S1PT, SVA\n")
}
cat("   - Note: PI-LL, L1-S1 Lordosis, and L1PA removed\n")
cat("   - Outcome: Change in Lordosis (6-week)\n\n")

cat("2. Knee Flexion Effect:\n")
cat(sprintf("   - Simple model: %.4f (p = %.4e)\n", 
            summary1$coefficients[2, 1], summary1$coefficients[2, 4]))
cat(sprintf("   - Adjusted model: %.4f (p = %.4e)\n", 
            summary2$coefficients[2, 1], summary2$coefficients[2, 4]))
cat("\n")

cat("3. Model Fit:\n")
cat(sprintf("   - Simple model R²: %.4f\n", summary1$r.squared))
cat(sprintf("   - Adjusted model R²: %.4f\n", summary2$r.squared))
cat("\n")

cat("4. Multicollinearity Assessment:\n")
max_vif <- max(vif_values[names(vif_values) != "LATpre_LL_KneeAngle"])
cat(sprintf("   - Max VIF: %.2f\n", max_vif))
if (max_vif < 2.5) {
  cat("   - Status: ✓ No multicollinearity concerns\n")
} else if (max_vif < 5) {
  cat("   - Status: ⚠ Mild multicollinearity (generally acceptable)\n")
} else if (max_vif < 10) {
  cat("   - Status: ⚠ Moderate multicollinearity (concerning)\n")
} else {
  cat("   - Status: ⚠ Severe multicollinearity (problematic)\n")
}
cat("\n")

cat("=== Analysis Complete ===\n")
cat("Check the plots in results/ directory for visualizations.\n\n")

