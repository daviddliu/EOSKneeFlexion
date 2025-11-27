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
# If TRUE, excludes patients from surgeons with 5 or fewer cases
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

# Check what preoperative PI-LL columns are available
# Common names: LATpre_PI_LL, PI_LL_pre, etc.
# If preop PI-LL doesn't exist, calculate it: PI - LL_preop
# PI would be a constant for each patient, so we might need to calculate: PI - LATpre_L1_S1
# Or look for a column that might contain this

cat("\n=== Searching for preoperative PI-LL mismatch variable ===\n")
pi_ll_cols <- grep("PI.*LL|LL.*PI", names(df), ignore.case = TRUE, value = TRUE)
cat("Potential PI-LL columns found:\n")
for (col in pi_ll_cols) {
  cat(sprintf("  %s: %d non-NA values\n", col, sum(!is.na(df[[col]]))))
}

# Try to find or calculate preop PI-LL mismatch
# First check if LATpre_PI_LL exists
if ("LATpre_PI_LL" %in% names(df)) {
  df$preop_PI_LL <- df$LATpre_PI_LL
  cat("Using LATpre_PI_LL as preoperative PI-LL mismatch\n")
} else {
  # Check if PI exists separately
  if ("PI" %in% names(df)) {
    # Calculate preop PI-LL mismatch as PI - LATpre_L1_S1
    df$preop_PI_LL <- df$PI - df$LATpre_L1_S1
    cat("Calculating preoperative PI-LL mismatch as PI - LATpre_L1_S1\n")
  } else {
    cat("WARNING: Could not find PI or preoperative PI-LL mismatch. Using 6-week postoperative PI-LL as proxy.\n")
    cat("This analysis may not accurately test your hypothesis.\n")
    df$preop_PI_LL <- df$LAT6W_PI_LL  # Use 6-week postop as proxy (not ideal)
  }
}

# Calculate T4-L1 PA (LATpre_L1PA - LATpre_T4PA)
if ("LATpre_L1PA" %in% names(df) && "LATpre_T4PA" %in% names(df)) {
  df$LATpre_T4_L1_PA <- df$LATpre_L1PA - df$LATpre_T4PA
  cat("\nCalculated LATpre_T4_L1_PA = LATpre_L1PA - LATpre_T4PA\n")
} else {
  cat("\nWARNING: Could not find LATpre_L1PA or LATpre_T4PA columns.\n")
  cat("Skipping LATpre_T4_L1_PA from analysis.\n")
  df$LATpre_T4_L1_PA <- NA
}

# Drop patients with missing data for confounder analysis
# A priori confounder set based on causal reasoning:
# PI-LL mismatch, Preop Lordosis (L1-S1), L4-S1, Thoracic Kyphosis, S1PT, SVA, T4-L1 PA
df_clean <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(preop_PI_LL) &
         !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) & 
         !is.na(LATpre_S1PT) & !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA))

cat(paste("\nAnalysis includes", nrow(df_clean), "patients with complete data (using 6-week follow-up)\n"))
cat("Variables included in full model (a priori confounder set):\n")
cat("  - Preop Knee Flexion (LATpre_LL_KneeAngle)\n")
cat("  - Preop PI-LL Mismatch (preop_PI_LL)\n")
cat("  - Preop Lordosis L1-S1 (LATpre_L1_S1)\n")
cat("  - Preop L4-S1 (LATpre_L4_S1)\n")
cat("  - Preop Thoracic Kyphosis (LATpre_T2_T12)\n")
cat("  - Preop S1PT - Pelvic Tilt (LATpre_S1PT)\n")
cat("  - Preop SVA (LATpre_SVA_C2_S1)\n")
cat("  - Preop T4-L1 PA (LATpre_T4_L1_PA = LATpre_L1PA - LATpre_T4PA)\n")
cat("  - Change in Lordosis (LAT6W_L1_S1 - LATpre_L1_S1, 6-week follow-up)\n")

# ============================================================================
# ANALYSIS 1: Simple regression (preop knee flexion ~ change in lordosis)
# ============================================================================
cat("\n=== Model 1: Simple Regression (Preop Knee Flexion Only) ===\n")
model1 <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = df_clean)
summary1 <- summary(model1)
cat("Coefficient for Preop Knee Flexion:\n")
cat(sprintf("  Estimate: %.4f\n", summary1$coefficients[2, 1]))
cat(sprintf("  P-value: %.4e\n", summary1$coefficients[2, 4]))
cat(sprintf("  R²: %.4f\n", summary1$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n", summary1$adj.r.squared))

# ============================================================================
# ANALYSIS 2: Multiple regression (adjusting for all confounders)
# ============================================================================
cat("\n=== Model 2: Multiple Regression (A Priori Confounders) ===\n")
cat("Adjusting for: Preop PI-LL Mismatch + Preop Lordosis (L1-S1) + Preop L4-S1 + Preop Thoracic Kyphosis + Preop S1PT + Preop SVA + Preop T4-L1 PA\n")
cat("(Confounders selected a priori based on causal reasoning)\n")
model2 <- lm(change_lordosis ~ LATpre_LL_KneeAngle + preop_PI_LL + LATpre_L1_S1 + 
             LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA, data = df_clean)
summary2 <- summary(model2)

cat("\nCoefficient for Preop Knee Flexion (adjusted for all covariates):\n")
cat(sprintf("  Estimate: %.4f\n", summary2$coefficients[2, 1]))
cat(sprintf("  P-value: %.4e\n", summary2$coefficients[2, 4]))
# Print coefficients, finding them by name to handle dropped variables
coef_names <- rownames(summary2$coefficients)
get_coef <- function(name) {
  idx <- which(coef_names == name)
  if (length(idx) == 1) {
    return(c(summary2$coefficients[idx, 1], summary2$coefficients[idx, 4]))
  }
  return(NULL)
}

cat("\nCoefficient for Preop PI-LL Mismatch:\n")
coef_pill <- get_coef("preop_PI_LL")
if (!is.null(coef_pill)) {
  cat(sprintf("  Estimate: %.4f\n", coef_pill[1]))
  cat(sprintf("  P-value: %.4e\n", coef_pill[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity)\n")
}

cat("\nCoefficient for Preop Lordosis L1-S1 (LATpre_L1_S1):\n")
coef_l1s1 <- get_coef("LATpre_L1_S1")
if (!is.null(coef_l1s1)) {
  cat(sprintf("  Estimate: %.4f\n", coef_l1s1[1]))
  cat(sprintf("  P-value: %.4e\n", coef_l1s1[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity)\n")
}

cat("\nCoefficient for Preop L4-S1 (LATpre_L4_S1):\n")
coef_l4s1 <- get_coef("LATpre_L4_S1")
if (!is.null(coef_l4s1)) {
  cat(sprintf("  Estimate: %.4f\n", coef_l4s1[1]))
  cat(sprintf("  P-value: %.4e\n", coef_l4s1[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity)\n")
}

cat("\nCoefficient for Preop Thoracic Kyphosis (LATpre_T2_T12):\n")
coef_t2t12 <- get_coef("LATpre_T2_T12")
if (!is.null(coef_t2t12)) {
  cat(sprintf("  Estimate: %.4f\n", coef_t2t12[1]))
  cat(sprintf("  P-value: %.4e\n", coef_t2t12[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity)\n")
}

cat("\nCoefficient for Preop S1PT - Pelvic Tilt (LATpre_S1PT):\n")
coef_s1pt <- get_coef("LATpre_S1PT")
if (!is.null(coef_s1pt)) {
  cat(sprintf("  Estimate: %.4f\n", coef_s1pt[1]))
  cat(sprintf("  P-value: %.4e\n", coef_s1pt[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity)\n")
}

cat("\nCoefficient for Preop SVA (LATpre_SVA_C2_S1):\n")
coef_sva <- get_coef("LATpre_SVA_C2_S1")
if (!is.null(coef_sva)) {
  cat(sprintf("  Estimate: %.4f\n", coef_sva[1]))
  cat(sprintf("  P-value: %.4e\n", coef_sva[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity or insufficient data)\n")
}

cat("\nCoefficient for Preop T4-L1 PA (LATpre_T4_L1_PA):\n")
coef_t4l1pa <- get_coef("LATpre_T4_L1_PA")
if (!is.null(coef_t4l1pa)) {
  cat(sprintf("  Estimate: %.4f\n", coef_t4l1pa[1]))
  cat(sprintf("  P-value: %.4e\n", coef_t4l1pa[2]))
} else {
  cat("  Variable was dropped from model (likely due to collinearity or insufficient data)\n")
}

cat(sprintf("\nR²: %.4f\n", summary2$r.squared))
cat(sprintf("Adjusted R²: %.4f\n", summary2$adj.r.squared))
cat(sprintf("\nR² change from Model 1 to Model 2: %.4f (%.2f%% increase)\n", 
            summary2$r.squared - summary1$r.squared,
            (summary2$r.squared - summary1$r.squared) / summary1$r.squared * 100))

# ============================================================================
# COMPARISON: Check if confounder effect
# ============================================================================
cat("\n=== Confounder Analysis ===\n")
coef1_kf <- summary1$coefficients[2, 1]
coef2_kf <- summary2$coefficients[2, 1]
pval1_kf <- summary1$coefficients[2, 4]
pval2_kf <- summary2$coefficients[2, 4]

# Get PI-LL coefficient by name
coef2_pill <- get_coef("preop_PI_LL")
if (!is.null(coef2_pill)) {
  pval2_pill <- coef2_pill[2]
} else {
  coef2_pill <- c(NA, NA)
  pval2_pill <- NA
}

change_percent <- abs((coef2_kf - coef1_kf) / coef1_kf * 100)

cat(sprintf("Change in Knee Flexion coefficient:\n"))
cat(sprintf("  Unadjusted: %.4f (p = %.4e)\n", coef1_kf, pval1_kf))
cat(sprintf("  Adjusted for all covariates: %.4f (p = %.4e)\n", coef2_kf, pval2_kf))
cat(sprintf("  Percent change: %.2f%%\n", change_percent))

cat(sprintf("\nAll covariate coefficients (adjusted model):\n"))
if (!is.null(get_coef("preop_PI_LL"))) {
  coef_pill <- get_coef("preop_PI_LL")
  cat(sprintf("  PI-LL Mismatch: %.4f (p = %.4e)\n", coef_pill[1], coef_pill[2]))
}
if (!is.null(coef_l1s1)) {
  cat(sprintf("  Preop Lordosis L1-S1: %.4f (p = %.4e)\n", coef_l1s1[1], coef_l1s1[2]))
}
if (!is.null(coef_l4s1)) {
  cat(sprintf("  Preop L4-S1: %.4f (p = %.4e)\n", coef_l4s1[1], coef_l4s1[2]))
}
if (!is.null(coef_t2t12)) {
  cat(sprintf("  Preop Thoracic Kyphosis: %.4f (p = %.4e)\n", coef_t2t12[1], coef_t2t12[2]))
}
if (!is.null(coef_s1pt)) {
  cat(sprintf("  Preop S1PT - Pelvic Tilt: %.4f (p = %.4e)\n", coef_s1pt[1], coef_s1pt[2]))
}
if (!is.null(coef_sva)) {
  cat(sprintf("  Preop SVA: %.4f (p = %.4e)\n", coef_sva[1], coef_sva[2]))
}
if (!is.null(coef_t4l1pa)) {
  cat(sprintf("  Preop T4-L1 PA: %.4f (p = %.4e)\n", coef_t4l1pa[1], coef_t4l1pa[2]))
}

cat("\n=== Interpretation ===\n")
if (abs(change_percent) > 10 && pval2_kf > 0.05) {
  cat("CONCLUSION: Knee flexion effect is explained by confounders.\n")
  cat("  - The knee flexion coefficient changed substantially and is no longer significant\n")
  cat("  - After adjusting for covariates, knee flexion does not independently predict change in lordosis\n")
} else if (abs(change_percent) > 10 && pval2_kf < 0.05) {
  cat("CONCLUSION: Knee flexion effect is attenuated but persists.\n")
  cat("  - The knee flexion coefficient changed substantially but remains significant\n")
  cat("  - Knee flexion has an independent effect on change in lordosis\n")
} else if (pval2_kf < 0.05) {
  cat("CONCLUSION: Knee flexion effect persists after adjustment.\n")
  cat("  - Knee flexion remains significant after controlling for all covariates\n")
  cat("  - Knee flexion has an independent effect on change in lordosis\n")
} else {
  cat("CONCLUSION: Knee flexion effect is not significant after adjustment.\n")
  cat("  - After adjusting for covariates, knee flexion does not independently predict change in lordosis\n")
}

# ============================================================================
# PARTIAL CORRELATION ANALYSIS
# ============================================================================
cat("\n=== Partial Correlation Analysis ===\n")

# Zero-order correlation (unadjusted)
r_kf_lordosis <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$change_lordosis, use = "complete.obs")
r_pill_lordosis <- cor(df_clean$preop_PI_LL, df_clean$change_lordosis, use = "complete.obs")
r_kf_pill <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$preop_PI_LL, use = "complete.obs")

cat(sprintf("Zero-order correlations:\n"))
cat(sprintf("  Knee Flexion ~ Change in Lordosis: r = %.4f\n", r_kf_lordosis))
cat(sprintf("  PI-LL Mismatch ~ Change in Lordosis: r = %.4f\n", r_pill_lordosis))
cat(sprintf("  Knee Flexion ~ PI-LL Mismatch: r = %.4f\n", r_kf_pill))

# Calculate correlations with all included variables
r_l1s1_lordosis <- cor(df_clean$LATpre_L1_S1, df_clean$change_lordosis, use = "complete.obs")
r_l4s1_lordosis <- cor(df_clean$LATpre_L4_S1, df_clean$change_lordosis, use = "complete.obs")
r_t2t12_lordosis <- cor(df_clean$LATpre_T2_T12, df_clean$change_lordosis, use = "complete.obs")
r_s1pt_lordosis <- cor(df_clean$LATpre_S1PT, df_clean$change_lordosis, use = "complete.obs")
r_sva_lordosis <- cor(df_clean$LATpre_SVA_C2_S1, df_clean$change_lordosis, use = "complete.obs")
r_t4l1pa_lordosis <- cor(df_clean$LATpre_T4_L1_PA, df_clean$change_lordosis, use = "complete.obs")

cat(sprintf("  Preop Lordosis L1-S1 ~ Change in Lordosis: r = %.4f\n", r_l1s1_lordosis))
cat(sprintf("  Preop L4-S1 ~ Change in Lordosis: r = %.4f\n", r_l4s1_lordosis))
cat(sprintf("  Preop Thoracic Kyphosis ~ Change in Lordosis: r = %.4f\n", r_t2t12_lordosis))
cat(sprintf("  Preop S1PT ~ Change in Lordosis: r = %.4f\n", r_s1pt_lordosis))
cat(sprintf("  Preop SVA ~ Change in Lordosis: r = %.4f\n", r_sva_lordosis))
cat(sprintf("  Preop T4-L1 PA ~ Change in Lordosis: r = %.4f\n", r_t4l1pa_lordosis))

# Partial correlation: r(knee flexion, change in lordosis | all other variables)
# With multiple variables, we use the regression coefficient approach
# Partial correlation from regression: r_partial = t / sqrt(t² + df)
# where t is the t-statistic from the regression
t_stat_kf <- summary2$coefficients[2, 3]  # t-statistic for knee flexion
df_residual <- summary2$df[2]  # residual degrees of freedom
r_partial <- t_stat_kf / sqrt(t_stat_kf^2 + df_residual)

# Partial R² = R² change when adding knee flexion to model (standard definition)
# This is the increase in R² from model without KF to model with KF
model_no_kf <- lm(change_lordosis ~ preop_PI_LL + LATpre_L1_S1 + 
                  LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA, 
                  data = df_clean)
summary_no_kf <- summary(model_no_kf)
r2_no_kf <- summary_no_kf$r.squared
r2_full <- summary2$r.squared
r2_partial <- r2_full - r2_no_kf  # Partial R² = R² change

cat(sprintf("\nPartial correlation (Knee Flexion ~ Change in Lordosis, controlling for all covariates):\n"))
cat(sprintf("  Partial r = %.4f\n", r_partial))
cat(sprintf("  Partial R² = R² change = %.4f (%.2f%% of variance)\n", r2_partial, r2_partial * 100))
cat(sprintf("    Full model R² = %.4f\n", r2_full))
cat(sprintf("    Model without KF R² = %.4f\n", r2_no_kf))

if (abs(r_partial) < abs(r_kf_lordosis) * 0.7) {
  cat("\n  -> Partial correlation is substantially reduced compared to zero-order correlation\n")
  cat("     This suggests the covariates (PI-LL, Lordosis L1-S1, L4-S1, Thoracic Kyphosis, S1PT, SVA, T4-L1 PA) are confounders\n")
}

# ============================================================================
# VISUALIZATION: Create plots
# ============================================================================
if (!dir.exists("results")) {
  dir.create("results")
}

# Plot 1: Unadjusted relationship
p1 <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = change_lordosis)) +
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
                   "\np = ", formatC(pval1_kf, format = "e", digits = 2)),
    hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold"
  )

ggsave("results/analysis4_unadjusted.png", plot = p1, width = 8, height = 6, dpi = 300)
cat("\nSaved unadjusted plot to results/analysis4_unadjusted.png\n")

# Plot 2: Residuals from full model (excluding knee flexion) vs Knee Flexion
# This shows the relationship between knee flexion and lordosis AFTER removing effects of all other covariates
# Using a priori confounder set: PI-LL, Lordosis L1-S1, L4-S1, Thoracic Kyphosis, S1PT, SVA, T4-L1 PA
covariates_model <- lm(change_lordosis ~ preop_PI_LL + LATpre_L1_S1 + 
                       LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA, data = df_clean)
df_clean$resid_from_covariates <- residuals(covariates_model)

kf_model_resid <- lm(resid_from_covariates ~ LATpre_LL_KneeAngle, data = df_clean)
summary_resid <- summary(kf_model_resid)

p2 <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = resid_from_covariates)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Preoperative Knee Flexion",
    y = "Residuals from Covariates Model\n(Change in Lordosis after removing\neffects of PI-LL, Lordosis L1-S1, L4-S1,\nThoracic Kyphosis, S1PT, SVA, T4-L1 PA)",
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
                   "\np = ", formatC(pval2_kf, format = "e", digits = 2)),
    hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold"
  )

ggsave("results/analysis4_adjusted.png", plot = p2, width = 8, height = 6, dpi = 300)
cat("Saved adjusted (residual) plot to results/analysis4_adjusted.png\n")

# Plot 3: Comparison of coefficients
coef_df <- data.frame(
  Model = c("Unadjusted", "Adjusted for PI-LL"),
  Coefficient = c(coef1_kf, coef2_kf),
  P_value = c(pval1_kf, pval2_kf),
  Lower_CI = c(confint(model1)[2, 1], confint(model2)[2, 1]),
  Upper_CI = c(confint(model1)[2, 2], confint(model2)[2, 2])
)

p3 <- ggplot(coef_df, aes(x = Model, y = Coefficient)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Model",
    y = "Coefficient for Preop Knee Flexion",
    title = "Comparison: Unadjusted vs Adjusted Coefficient (6-week follow-up)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text", x = coef_df$Model, y = coef_df$Coefficient,
    label = paste0("p = ", formatC(coef_df$P_value, format = "e", digits = 2)),
    vjust = -1.5, size = 3
  )

ggsave("results/analysis4_coefficient_comparison.png", plot = p3, width = 8, height = 6, dpi = 300)
cat("Saved coefficient comparison plot to results/analysis4_coefficient_comparison.png\n")

cat("\n=== Analysis Complete ===\n")
cat("Check the plots and statistics above to evaluate if PI-LL mismatch is a confounder (using 6-week follow-up data).\n")

