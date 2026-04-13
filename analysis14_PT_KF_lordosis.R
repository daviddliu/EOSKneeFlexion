#!/usr/bin/env Rscript
#
# Analysis 14: PT vs Change in Lordosis, with and without KF
#
# This script compares:
#   Model 1: Change in Lordosis ~ Pelvic Tilt (PT only)
#   Model 2: Change in Lordosis ~ Pelvic Tilt + Knee Flexion (PT + KF)
#
# Hypothesis: Adding KF should not dramatically improve the model fit,
#             suggesting PT is the primary predictor.

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

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Calculate change in lordosis using 6-week data
df$change_lordosis <- df$LAT6W_L1_S1 - df$LATpre_L1_S1

# Drop patients with missing data
df_clean <- df %>%
  filter(!is.na(LATpre_S1PT) & !is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis))

cat(paste("\n=== Analysis 14: PT vs Change in Lordosis (with and without KF) ===\n"))
cat(paste("Sample size:", nrow(df_clean), "patients with complete data (using 6-week follow-up)\n\n"))

# ============================================================================
# Model 1: Simple regression (PT only)
# ============================================================================
cat("=== Model 1: Change in Lordosis ~ Pelvic Tilt (PT only) ===\n\n")

model1 <- lm(change_lordosis ~ LATpre_S1PT, data = df_clean)
summary1 <- summary(model1)

cat("Regression Results:\n")
cat(sprintf("  Pelvic Tilt coefficient: %.4f (SE = %.4f)\n", 
            summary1$coefficients[2, 1], summary1$coefficients[2, 2]))
cat(sprintf("  p-value: %.4e\n", summary1$coefficients[2, 4]))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
            confint(model1)[2, 1], confint(model1)[2, 2]))
cat(sprintf("  R²: %.4f\n", summary1$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n", summary1$adj.r.squared))
cat(sprintf("  Residual SE: %.4f\n", summary1$sigma))

# Calculate correlation
r_pt_lordosis <- cor(df_clean$LATpre_S1PT, df_clean$change_lordosis, use = "complete.obs")
cat(sprintf("  Pearson r: %.4f\n\n", r_pt_lordosis))

# ============================================================================
# Model 2: Multiple regression (PT + KF)
# ============================================================================
cat("=== Model 2: Change in Lordosis ~ Pelvic Tilt + Knee Flexion (PT + KF) ===\n\n")

model2 <- lm(change_lordosis ~ LATpre_S1PT + LATpre_LL_KneeAngle, data = df_clean)
summary2 <- summary(model2)

cat("Regression Results:\n")
cat(sprintf("  Pelvic Tilt coefficient: %.4f (SE = %.4f)\n", 
            summary2$coefficients[2, 1], summary2$coefficients[2, 2]))
cat(sprintf("  Pelvic Tilt p-value: %.4e\n", summary2$coefficients[2, 4]))
cat(sprintf("  Pelvic Tilt 95%% CI: [%.4f, %.4f]\n", 
            confint(model2)[2, 1], confint(model2)[2, 2]))
cat(sprintf("  Knee Flexion coefficient: %.4f (SE = %.4f)\n", 
            summary2$coefficients[3, 1], summary2$coefficients[3, 2]))
cat(sprintf("  Knee Flexion p-value: %.4e\n", summary2$coefficients[3, 4]))
cat(sprintf("  Knee Flexion 95%% CI: [%.4f, %.4f]\n", 
            confint(model2)[3, 1], confint(model2)[3, 2]))
cat(sprintf("  R²: %.4f\n", summary2$r.squared))
cat(sprintf("  Adjusted R²: %.4f\n", summary2$adj.r.squared))
cat(sprintf("  Residual SE: %.4f\n\n", summary2$sigma))

# ============================================================================
# Model Comparison
# ============================================================================
cat("=== Model Comparison ===\n\n")

# Calculate R² change
r2_change <- summary2$r.squared - summary1$r.squared
r2_change_pct <- (r2_change / summary1$r.squared) * 100

cat("Model 1 (PT only):\n")
cat(sprintf("  R² = %.4f\n", summary1$r.squared))
cat(sprintf("  Adjusted R² = %.4f\n", summary1$adj.r.squared))
cat("\n")

cat("Model 2 (PT + KF):\n")
cat(sprintf("  R² = %.4f\n", summary2$r.squared))
cat(sprintf("  Adjusted R² = %.4f\n", summary2$adj.r.squared))
cat("\n")

cat("Improvement from adding KF:\n")
cat(sprintf("  R² change: %.4f (absolute)\n", r2_change))
cat(sprintf("  R² change: %.2f%% (relative to Model 1)\n", r2_change_pct))
cat("\n")

# F-test for model comparison
anova_result <- anova(model1, model2)
cat("F-test for model comparison:\n")
cat(sprintf("  F-statistic: %.4f\n", anova_result$F[2]))
cat(sprintf("  p-value: %.4e\n", anova_result$`Pr(>F)`[2]))
if (anova_result$`Pr(>F)`[2] < 0.05) {
  cat("  → Adding KF significantly improves the model (p < 0.05)\n")
} else {
  cat("  → Adding KF does NOT significantly improve the model (p >= 0.05)\n")
}
cat("\n")

# Partial R² for KF (contribution of KF after controlling for PT)
# This is the R² change when adding KF to the PT-only model
cat("Partial R² for Knee Flexion (contribution after controlling for PT):\n")
cat(sprintf("  Partial R² = R² change = %.4f (%.2f%% of total variance)\n", 
            r2_change, r2_change * 100))
cat("\n")

# ============================================================================
# Visualization
# ============================================================================
cat("=== Creating Visualizations ===\n\n")

if (!dir.exists("results")) {
  dir.create("results")
}

# Plot 1: PT vs Change in Lordosis (Model 1)
png("results/analysis14_PT_vs_change_lordosis.png", width = 10, height = 8, units = "in", res = 300)
p1 <- ggplot(df_clean, aes(x = LATpre_S1PT, y = change_lordosis)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative Pelvic Tilt (degrees)",
    y = "Change in Lordosis (6-week - Preop, degrees)",
    title = "Model 1: Pelvic Tilt vs Change in Lordosis"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("r = ", round(r_pt_lordosis, 3),
                   "\nR² = ", round(summary1$r.squared, 3),
                   "\nCoefficient = ", round(summary1$coefficients[2, 1], 4),
                   "\np = ", formatC(summary1$coefficients[2, 4], format = "e", digits = 2),
                   "\nn = ", nrow(df_clean)),
    hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
  )

print(p1)
dev.off()
cat("Saved plot to results/analysis14_PT_vs_change_lordosis.png\n")

# Plot 2: Residuals from PT model vs KF (showing what KF adds)
png("results/analysis14_residuals_PT_vs_KF.png", width = 10, height = 8, units = "in", res = 300)

# Calculate residuals from PT-only model
df_clean$resid_from_PT <- residuals(model1)

# Fit regression of residuals on KF
resid_model <- lm(resid_from_PT ~ LATpre_LL_KneeAngle, data = df_clean)
summary_resid <- summary(resid_model)
r_resid_kf <- cor(df_clean$resid_from_PT, df_clean$LATpre_LL_KneeAngle, use = "complete.obs")

p2 <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = resid_from_PT)) +
  geom_point(alpha = 0.6, color = "darkgrey") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Preoperative Knee Flexion (degrees)",
    y = "Residuals from PT Model\n(Change in Lordosis after removing effect of PT)",
    title = "What Does KF Add?\n(Residuals from PT model vs KF)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("r = ", round(r_resid_kf, 3),
                   "\nR² = ", round(summary_resid$r.squared, 3),
                   "\nCoefficient = ", round(summary_resid$coefficients[2, 1], 4),
                   "\np = ", formatC(summary_resid$coefficients[2, 4], format = "e", digits = 2),
                   "\n\nPartial R² (KF) = ", round(r2_change, 4)),
    hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
  )

print(p2)
dev.off()
cat("Saved plot to results/analysis14_residuals_PT_vs_KF.png\n")

# Plot 3: Model comparison (coefficients with confidence intervals)
png("results/analysis14_model_comparison.png", width = 12, height = 6, units = "in", res = 300)

comparison_df <- data.frame(
  Model = rep(c("Model 1 (PT only)", "Model 2 (PT + KF)"), each = 2),
  Variable = c("Pelvic Tilt", "Intercept", "Pelvic Tilt", "Knee Flexion"),
  Coefficient = c(
    summary1$coefficients[2, 1],
    summary1$coefficients[1, 1],
    summary2$coefficients[2, 1],
    summary2$coefficients[3, 1]
  ),
  Lower_CI = c(
    confint(model1)[2, 1],
    confint(model1)[1, 1],
    confint(model2)[2, 1],
    confint(model2)[3, 1]
  ),
  Upper_CI = c(
    confint(model1)[2, 2],
    confint(model1)[1, 2],
    confint(model2)[2, 2],
    confint(model2)[3, 2]
  ),
  P_value = c(
    summary1$coefficients[2, 4],
    summary1$coefficients[1, 4],
    summary2$coefficients[2, 4],
    summary2$coefficients[3, 4]
  )
)

# Only plot PT coefficients for comparison
pt_comparison <- comparison_df %>%
  filter(Variable == "Pelvic Tilt")

p3 <- ggplot(pt_comparison, aes(x = Model, y = Coefficient)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Model",
    y = "Pelvic Tilt Coefficient",
    title = "Model Comparison: Pelvic Tilt Coefficient\n(with and without KF in the model)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate(
    "text", x = pt_comparison$Model, y = pt_comparison$Coefficient,
    label = paste0("p = ", formatC(pt_comparison$P_value, format = "e", digits = 2)),
    vjust = -1.5, size = 3
  )

print(p3)
dev.off()
cat("Saved plot to results/analysis14_model_comparison.png\n")

# Plot 4: R² comparison
png("results/analysis14_R2_comparison.png", width = 10, height = 6, units = "in", res = 300)

r2_df <- data.frame(
  Model = c("Model 1\n(PT only)", "Model 2\n(PT + KF)"),
  R_squared = c(summary1$r.squared, summary2$r.squared),
  Adj_R_squared = c(summary1$adj.r.squared, summary2$adj.r.squared)
)

p4 <- ggplot(r2_df, aes(x = Model)) +
  geom_bar(aes(y = R_squared), stat = "identity", fill = "steelblue", alpha = 0.7, width = 0.5) +
  geom_bar(aes(y = Adj_R_squared), stat = "identity", fill = "darkblue", alpha = 0.7, width = 0.5, position = position_nudge(x = 0.5)) +
  geom_text(aes(y = R_squared, label = sprintf("%.3f", R_squared)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  geom_text(aes(y = Adj_R_squared, label = sprintf("%.3f", Adj_R_squared)), 
            vjust = -0.5, size = 4, fontface = "bold", position = position_nudge(x = 0.5)) +
  labs(
    x = "Model",
    y = "R²",
    title = "Model Comparison: R² Values"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = 1.5, y = max(r2_df$R_squared) * 0.95,
           label = paste0("R² change = ", round(r2_change, 4), "\n(", round(r2_change_pct, 2), "% increase)"),
           size = 4, fontface = "bold", color = "darkred")

print(p4)
dev.off()
cat("Saved plot to results/analysis14_R2_comparison.png\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n=== SUMMARY ===\n\n")

cat("Model 1 (PT only):\n")
cat(sprintf("  - Pelvic Tilt explains %.2f%% of variance in change in lordosis\n", summary1$r.squared * 100))
cat(sprintf("  - PT coefficient: %.4f (p = %.4e)\n", 
            summary1$coefficients[2, 1], summary1$coefficients[2, 4]))
cat("\n")

cat("Model 2 (PT + KF):\n")
cat(sprintf("  - Combined model explains %.2f%% of variance\n", summary2$r.squared * 100))
cat(sprintf("  - PT coefficient: %.4f (p = %.4e)\n", 
            summary2$coefficients[2, 1], summary2$coefficients[2, 4]))
cat(sprintf("  - KF coefficient: %.4f (p = %.4e)\n", 
            summary2$coefficients[3, 1], summary2$coefficients[3, 4]))
cat("\n")

cat("Contribution of KF:\n")
cat(sprintf("  - KF adds %.2f%% to R² (absolute: %.4f)\n", r2_change_pct, r2_change))
cat(sprintf("  - Partial R² for KF: %.4f (%.2f%% of total variance)\n", r2_change, r2_change * 100))
if (r2_change < 0.05) {
  cat("  → KF adds relatively little (< 5% absolute increase in R²)\n")
  cat("  → This suggests PT is the primary predictor\n")
} else {
  cat("  → KF adds meaningful variance explanation (>= 5% absolute increase)\n")
}
cat("\n")

if (anova_result$`Pr(>F)`[2] < 0.05) {
  cat("Statistical test:\n")
  cat("  → Adding KF significantly improves the model (F-test p < 0.05)\n")
  cat("  → However, the practical improvement may still be modest\n")
} else {
  cat("Statistical test:\n")
  cat("  → Adding KF does NOT significantly improve the model (F-test p >= 0.05)\n")
  cat("  → This supports that PT is the primary predictor\n")
}
cat("\n")

cat("=== Analysis Complete ===\n")
cat("Check the plots in results/ directory for visualizations.\n\n")
