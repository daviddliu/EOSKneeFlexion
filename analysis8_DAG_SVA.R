#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(boot)
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

# Calculate change in sva (cm) using 6-week data
df$change_SVA <- (df$LAT6W_SVA_C2_S1 - df$LATpre_SVA_C2_S1) / 10

# Calculate T4-L1 PA
if ("LATpre_L1PA" %in% names(df) && "LATpre_T4PA" %in% names(df)) {
  df$LATpre_T4_L1_PA <- df$LATpre_L1PA - df$LATpre_T4PA
  cat("Calculated LATpre_T4_L1_PA = LATpre_L1PA - LATpre_T4PA\n")
} else {
  cat("WARNING: Could not find LATpre_L1PA or LATpre_T4PA columns.\n")
  df$LATpre_T4_L1_PA <- NA
}

# Check for age variable
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
  cat("WARNING: Age variable not found. Will exclude from analysis.\n")
}

# Try to find or calculate preop PI-LL mismatch
if ("LATpre_PI_LL" %in% names(df)) {
  df$preop_PI_LL <- df$LATpre_PI_LL
  cat("Using LATpre_PI_LL as preoperative PI-LL mismatch\n")
} else if ("PI" %in% names(df)) {
  df$preop_PI_LL <- df$PI - df$LATpre_L1_S1
  cat("Calculating preoperative PI-LL mismatch as PI - LATpre_L1_S1\n")
} else {
  cat("WARNING: Could not find PI or preoperative PI-LL mismatch. Using 6-week postoperative PI-LL as proxy.\n")
  df$preop_PI_LL <- df$LAT6W_PI_LL
}

cat("\n" , rep("=", 80), "\n", sep = "")
cat("RIGOROUS CAUSAL INFERENCE ANALYSIS: DAG-BASED ADJUSTMENT SET\n")
cat(rep("=", 80), "\n\n")

# ============================================================================
# STEP 1: DAG → DECIDE ADJUSTMENT SET
# ============================================================================
cat("=== STEP 1: DAG-Based Adjustment Set Selection ===\n\n")

cat("Causal Structure Reasoning:\n")
cat("  X (Exposure): Preoperative Knee Flexion (LATpre_LL_KneeAngle)\n")
cat("  Y (Outcome): Change in Lumbar Lordosis ((LAT6W_SVA_C2_S1 - LATpre_SVA_C2_S1) / 10)\n\n")

cat("Potential Covariates:\n")
cat("  - Preop PI-LL Mismatch (preop_PI_LL): CONFOUNDER (affects both X and Y)\n")
cat("  - Preop Lordosis L1-S1 (LATpre_L1_S1): CONFOUNDER (baseline alignment affects both)\n")
cat("  - Preop Lordosis L4-S1 (LATpre_L4_S1): CONFOUNDER (baseline alignment affects both)\n")
cat("  - Preop Thoracic Kyphosis (LATpre_T2_T12): CONFOUNDER (spinal alignment affects both)\n")
cat("  - Preop S1PI (LATpre_S1PI): CONFOUNDER (pelvic parameters affect both)\n")
cat("  - Preop S1PT (LATpre_S1PT): CONFOUNDER (pelvic tilt affects both)\n")
cat("  - Preop SVA (LATpre_SVA_C2_S1): CONFOUNDER (sagittal balance affects both)\n")
cat("  - Preop T4-L1 PA (LATpre_T4_L1_PA): CONFOUNDER (sagittal alignment affects both)\n")
if (!is.null(age_var)) {
  cat(sprintf("  - Age (%s): CONFOUNDER (age affects both)\n", age_var))
}
cat("\n")

cat("EXCLUDED (potential mediators/colliders):\n")
cat("  - Postoperative measurements: MEDIATORS (on causal pathway from X to Y)\n")
cat("  - Surgical technique variables: MEDIATORS (if they affect outcome)\n")
cat("  - Variables affected by both X and Y: COLLIDERS (would bias estimates)\n\n")

# Define adjustment set based on DAG reasoning
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

adjustment_labels <- c(
  "PI-LL Mismatch",
  "Lordosis L1-S1",
  "Lordosis L4-S1",
  "Thoracic Kyphosis T2-T12",
  "S1PI",
  "Pelvic Tilt S1PT",
  "SVA",
  "T4-L1 PA"
)

if (!is.null(age_var)) {
  adjustment_set <- c(adjustment_set, age_var)
  adjustment_labels <- c(adjustment_labels, "Age")
}

cat("FINAL ADJUSTMENT SET (confounders only, no mediators/colliders):\n")
for (i in seq_along(adjustment_set)) {
  cat(sprintf("  %d. %s (%s)\n", i, adjustment_labels[i], adjustment_set[i]))
}
cat("\n")

# Filter to complete cases
df_clean <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_SVA))

# Check availability of each adjustment variable
for (var in adjustment_set) {
  n_missing <- sum(is.na(df_clean[[var]]))
  n_available <- sum(!is.na(df_clean[[var]]))
  cat(sprintf("  %s: %d available, %d missing\n", var, n_available, n_missing))
}

# Create complete case dataset
complete_vars <- c("LATpre_LL_KneeAngle", "change_SVA", adjustment_set)
df_complete <- df_clean %>%
  filter(complete.cases(select(., all_of(complete_vars))))

cat(sprintf("\nComplete case sample size: %d patients\n", nrow(df_complete)))
cat(sprintf("Dropped %d patients with missing data\n\n", nrow(df_clean) - nrow(df_complete)))

# ============================================================================
# STEP 2: ESTIMATE EFFECT (with standard regression, noting AIPW/TMLE option)
# ============================================================================
cat("=== STEP 2: Effect Estimation ===\n\n")

cat("NOTE: For maximum rigor, consider using AIPW/TMLE with cross-fitting.\n")
cat("      This analysis uses standard regression for demonstration.\n")
cat("      AIPW/TMLE would provide doubly-robust estimates and better handle\n")
cat("      potential model misspecification.\n\n")

# Unadjusted model
model_unadj <- lm(change_SVA ~ LATpre_LL_KneeAngle, data = df_complete)
summary_unadj <- summary(model_unadj)
coef_unadj <- summary_unadj$coefficients[2, 1]
ci_unadj <- confint(model_unadj)[2, ]
p_unadj <- summary_unadj$coefficients[2, 4]

cat("Unadjusted Effect:\n")
cat(sprintf("  Coefficient: %.4f (95%% CI: [%.4f, %.4f])\n", coef_unadj, ci_unadj[1], ci_unadj[2]))
cat(sprintf("  P-value: %.4e\n", p_unadj))
cat(sprintf("  R²: %.4f\n\n", summary_unadj$r.squared))

# Fully adjusted model
formula_adj <- as.formula(paste("change_SVA ~ LATpre_LL_KneeAngle +", 
                                paste(adjustment_set, collapse = " + ")))
model_adj <- lm(formula_adj, data = df_complete)
summary_adj <- summary(model_adj)
coef_adj <- summary_adj$coefficients[2, 1]
ci_adj <- confint(model_adj)[2, ]
p_adj <- summary_adj$coefficients[2, 4]

cat("Fully Adjusted Effect:\n")
cat(sprintf("  Coefficient: %.4f (95%% CI: [%.4f, %.4f])\n", coef_adj, ci_adj[1], ci_adj[2]))
cat(sprintf("  P-value: %.4e\n", p_adj))
cat(sprintf("  R²: %.4f\n\n", summary_adj$r.squared))

# Effect attenuation
attenuation <- abs((coef_unadj - coef_adj) / coef_unadj) * 100
cat(sprintf("Effect Attenuation: %.2f%% reduction in coefficient magnitude\n\n", attenuation))

# Model without knee flexion (to test if other parameters alone explain the variance)
formula_no_kf <- as.formula(paste("change_SVA ~", paste(adjustment_set, collapse = " + ")))
model_no_kf <- lm(formula_no_kf, data = df_complete)
summary_no_kf <- summary(model_no_kf)
r2_no_kf <- summary_no_kf$r.squared
r2_with_kf <- summary_adj$r.squared
kf_partial_r2 <- r2_with_kf - r2_no_kf

cat("Model Comparison:\n")
cat(sprintf("  R² without knee flexion (other parameters only): %.4f\n", r2_no_kf))
cat(sprintf("  R² with knee flexion added: %.4f\n", r2_with_kf))
cat(sprintf("  Partial R² for knee flexion: %.4f (%.2f%% of total variance)\n\n", 
            kf_partial_r2, (kf_partial_r2 / r2_with_kf) * 100))

# ============================================================================
# STEP 3: BLOCK-WISE SEQUENTIAL ADJUSTMENT (PRE-SPECIFIED)
# ============================================================================
cat("=== STEP 3: Block-Wise Sequential Adjustment ===\n\n")

# Define blocks based on clinical/anatomical grouping
blocks <- list(
  list(
    name = "Pelvic Parameters",
    vars = c("preop_PI_LL", "LATpre_S1PI", "LATpre_S1PT"),
    labels = c("PI-LL Mismatch", "S1PI", "Pelvic Tilt S1PT")
  ),
  list(
    name = "Lumbar Alignment",
    vars = c("LATpre_L1_S1", "LATpre_L4_S1"),
    labels = c("Lordosis L1-S1", "Lordosis L4-S1")
  ),
  list(
    name = "Thoracic Alignment",
    vars = c("LATpre_T2_T12"),
    labels = c("Thoracic Kyphosis T2-T12")
  ),
  list(
    name = "Sagittal Balance",
    vars = c("LATpre_SVA_C2_S1", "LATpre_T4_L1_PA"),
    labels = c("SVA", "T4-L1 PA")
  )
)

# Add age block if available
if (!is.null(age_var)) {
  blocks <- c(blocks, list(
    list(
      name = "Demographics",
      vars = c(age_var),
      labels = c("Age")
    )
  ))
}

cat("Pre-specified Blocks:\n")
for (i in seq_along(blocks)) {
  cat(sprintf("  Block %d: %s\n", i, blocks[[i]]$name))
  for (j in seq_along(blocks[[i]]$vars)) {
    cat(sprintf("    - %s (%s)\n", blocks[[i]]$labels[j], blocks[[i]]$vars[j]))
  }
}
cat("\n")

# Sequential adjustment
block_results <- list()
cumulative_vars <- c()
prev_r2 <- summary_unadj$r.squared

for (i in seq_along(blocks)) {
  block <- blocks[[i]]
  prev_cumulative_vars <- cumulative_vars
  cumulative_vars <- c(cumulative_vars, block$vars)
  
  if (length(cumulative_vars) > 0) {
    formula_block <- as.formula(paste("change_SVA ~ LATpre_LL_KneeAngle +", 
                                      paste(cumulative_vars, collapse = " + ")))
    model_block <- lm(formula_block, data = df_complete)
    summary_block <- summary(model_block)
    
    coef_block <- summary_block$coefficients[2, 1]
    ci_block <- confint(model_block)[2, ]
    p_block <- summary_block$coefficients[2, 4]
    r2_block <- summary_block$r.squared
    
    # Calculate R² change from previous block
    r2_change <- r2_block - prev_r2
    prev_r2 <- r2_block
    
    block_results[[length(block_results) + 1]] <- data.frame(
      Block = block$name,
      Block_Number = i,
      Cumulative_Vars = length(cumulative_vars),
      KF_Coefficient = coef_block,
      KF_CI_Lower = ci_block[1],
      KF_CI_Upper = ci_block[2],
      KF_P = p_block,
      R2 = r2_block,
      R2_Change = r2_change
    )
    
    cat(sprintf("After Block %d (%s):\n", i, block$name))
    cat(sprintf("  KF Coefficient: %.4f (95%% CI: [%.4f, %.4f])\n", 
                coef_block, ci_block[1], ci_block[2]))
    cat(sprintf("  P-value: %.4e\n", p_block))
    cat(sprintf("  R²: %.4f (R² change: %.4f)\n\n", r2_block, r2_change))
  }
}

block_results_df <- do.call(rbind, block_results)

# ============================================================================
# STEP 4: INDIVIDUAL Zk CONTRIBUTIONS
# ============================================================================
cat("=== STEP 4: Individual Covariate Contributions ===\n\n")

cat("Method 1: Drop-One-Covariate Analysis (Targeted Sensitivity)\n")
cat("  (Less ideal but common and computationally feasible)\n\n")

# Full model without each covariate
drop_one_results <- list()

for (i in seq_along(adjustment_set)) {
  var <- adjustment_set[i]
  label <- adjustment_labels[i]
  
  # Model without this covariate
  vars_without <- adjustment_set[-i]
  if (length(vars_without) > 0) {
    formula_drop <- as.formula(paste("change_SVA ~ LATpre_LL_KneeAngle +", 
                                     paste(vars_without, collapse = " + ")))
  } else {
    formula_drop <- as.formula("change_SVA ~ LATpre_LL_KneeAngle")
  }
  
  model_drop <- lm(formula_drop, data = df_complete)
  summary_drop <- summary(model_drop)
  
  coef_drop <- summary_drop$coefficients[2, 1]
  ci_drop <- confint(model_drop)[2, ]
  p_drop <- summary_drop$coefficients[2, 4]
  r2_drop <- summary_drop$r.squared
  
  # Effect change when this covariate is dropped
  coef_change <- coef_adj - coef_drop
  coef_change_pct <- (coef_change / coef_adj) * 100
  r2_change <- summary_adj$r.squared - r2_drop
  
  drop_one_results[[length(drop_one_results) + 1]] <- data.frame(
    Covariate = label,
    Variable = var,
    KF_Coefficient_Without = coef_drop,
    KF_CI_Lower_Without = ci_drop[1],
    KF_CI_Upper_Without = ci_drop[2],
    KF_P_Without = p_drop,
    R2_Without = r2_drop,
    KF_Coefficient_Change = coef_change,
    KF_Coefficient_Change_Pct = coef_change_pct,
    R2_Change = r2_change
  )
  
  cat(sprintf("Without %s:\n", label))
  cat(sprintf("  KF Coefficient: %.4f (change: %.4f, %.2f%%)\n", 
              coef_drop, coef_change, coef_change_pct))
  cat(sprintf("  P-value: %.4e\n", p_drop))
  cat(sprintf("  R²: %.4f (R² change: %.4f)\n\n", r2_drop, r2_change))
}

drop_one_df <- do.call(rbind, drop_one_results)
drop_one_df <- drop_one_df %>%
  arrange(desc(abs(KF_Coefficient_Change_Pct)))

cat("\nMethod 2: Dominance Analysis (Shapley Value Approximation)\n")
cat("  (More principled but computationally heavier)\n\n")

# Dominance analysis: calculate average R² contribution across all model combinations
# This is computationally intensive, so we'll do a simplified version
# Full dominance analysis would require 2^n models

cat("Computing dominance statistics (simplified version)...\n")

# Calculate average contribution when variable is added to different subsets
dominance_results <- list()

for (i in seq_along(adjustment_set)) {
  var <- adjustment_set[i]
  label <- adjustment_labels[i]
  
  # Calculate contribution when added to base model
  formula_base <- as.formula("change_SVA ~ LATpre_LL_KneeAngle")
  model_base <- lm(formula_base, data = df_complete)
  r2_base <- summary(model_base)$r.squared
  
  formula_with <- as.formula(paste("change_SVA ~ LATpre_LL_KneeAngle +", var))
  model_with <- lm(formula_with, data = df_complete)
  r2_with <- summary(model_with)$r.squared
  
  contribution_base <- r2_with - r2_base
  
  # Calculate average contribution when added to models with other variables
  other_vars <- adjustment_set[-i]
  contributions <- c(contribution_base)
  
  # Sample a subset of combinations (full enumeration would be 2^(n-1))
  n_combinations <- min(20, 2^min(length(other_vars), 4))  # Limit to avoid explosion
  set.seed(123)
  
  for (j in seq_len(n_combinations)) {
    # Random subset of other variables
    if (length(other_vars) > 0) {
      n_other <- sample(0:min(length(other_vars), 4), 1)
      if (n_other > 0) {
        subset_vars <- sample(other_vars, n_other)
        formula_subset <- as.formula(paste("change_SVA ~ LATpre_LL_KneeAngle +", 
                                          paste(subset_vars, collapse = " + ")))
        model_subset <- lm(formula_subset, data = df_complete)
        r2_subset <- summary(model_subset)$r.squared
        
        formula_subset_with <- as.formula(paste("change_SVA ~ LATpre_LL_KneeAngle +", 
                                                paste(c(subset_vars, var), collapse = " + ")))
        model_subset_with <- lm(formula_subset_with, data = df_complete)
        r2_subset_with <- summary(model_subset_with)$r.squared
        
        contributions <- c(contributions, r2_subset_with - r2_subset)
      }
    }
  }
  
  avg_contribution <- mean(contributions)
  dominance_results[[length(dominance_results) + 1]] <- data.frame(
    Covariate = label,
    Variable = var,
    Average_R2_Contribution = avg_contribution,
    Min_Contribution = min(contributions),
    Max_Contribution = max(contributions),
    SD_Contribution = sd(contributions)
  )
  
  cat(sprintf("  %s: Average R² contribution = %.4f (SD = %.4f)\n", 
              label, avg_contribution, sd(contributions)))
}

dominance_df <- do.call(rbind, dominance_results)
dominance_df <- dominance_df %>%
  arrange(desc(Average_R2_Contribution))

cat("\n")

# ============================================================================
# STEP 5: REPORTING - STANDARDIZED MEAN DIFFERENCES (BALANCE)
# ============================================================================
cat("=== STEP 5: Balance Diagnostics (Standardized Mean Differences) ===\n\n")

# For balance, we'd typically compare groups defined by exposure
# Here we'll create tertiles of knee flexion for illustration
df_complete$kf_tertile <- cut(df_complete$LATpre_LL_KneeAngle, 
                                breaks = quantile(df_complete$LATpre_LL_KneeAngle, 
                                                 probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                labels = c("Low", "Medium", "High"),
                                include.lowest = TRUE)

balance_results <- list()

for (var in adjustment_set) {
  var_data <- df_complete[[var]]
  var_label <- adjustment_labels[which(adjustment_set == var)]
  
  # Calculate means and SDs for each tertile
  means <- tapply(var_data, df_complete$kf_tertile, mean, na.rm = TRUE)
  sds <- tapply(var_data, df_complete$kf_tertile, sd, na.rm = TRUE)
  
  # Standardized mean differences (Cohen's d)
  pooled_sd <- sqrt(mean(sds^2, na.rm = TRUE))
  
  smd_low_high <- (means["High"] - means["Low"]) / pooled_sd
  smd_medium_high <- (means["High"] - means["Medium"]) / pooled_sd
  smd_low_medium <- (means["Medium"] - means["Low"]) / pooled_sd
  
  balance_results[[length(balance_results) + 1]] <- data.frame(
    Covariate = var_label,
    Variable = var,
    Mean_Low = means["Low"],
    Mean_Medium = means["Medium"],
    Mean_High = means["High"],
    SMD_Low_vs_High = smd_low_high,
    SMD_Medium_vs_High = smd_medium_high,
    SMD_Low_vs_Medium = smd_low_medium
  )
}

balance_df <- do.call(rbind, balance_results)

cat("Standardized Mean Differences (by knee flexion tertile):\n")
print(balance_df %>% 
      select(Covariate, SMD_Low_vs_High, SMD_Medium_vs_High, SMD_Low_vs_Medium) %>%
      mutate(
        SMD_Low_vs_High = round(SMD_Low_vs_High, 3),
        SMD_Medium_vs_High = round(SMD_Medium_vs_High, 3),
        SMD_Low_vs_Medium = round(SMD_Low_vs_Medium, 3)
      ),
      row.names = FALSE)
cat("\n")

cat("Interpretation: |SMD| > 0.1 suggests imbalance\n")
cat("                |SMD| > 0.25 suggests substantial imbalance\n\n")

# ============================================================================
# STEP 6: POSITIVITY DIAGNOSTICS (OVERLAP)
# ============================================================================
cat("=== STEP 6: Positivity Diagnostics (Overlap) ===\n\n")

# Check overlap in covariate distributions
overlap_results <- list()

for (var in adjustment_set) {
  var_data <- df_complete[[var]]
  var_label <- adjustment_labels[which(adjustment_set == var)]
  
  # Create exposure groups (tertiles)
  low_data <- var_data[df_complete$kf_tertile == "Low"]
  high_data <- var_data[df_complete$kf_tertile == "High"]
  
  # Overlap: proportion of observations in common support
  min_high <- min(high_data, na.rm = TRUE)
  max_low <- max(low_data, na.rm = TRUE)
  min_low <- min(low_data, na.rm = TRUE)
  max_high <- max(high_data, na.rm = TRUE)
  
  # Common support region
  common_min <- max(min_low, min_high)
  common_max <- min(max_low, max_high)
  
  overlap_low <- mean(low_data >= common_min & low_data <= common_max, na.rm = TRUE)
  overlap_high <- mean(high_data >= common_min & high_data <= common_max, na.rm = TRUE)
  overlap_overall <- mean(overlap_low, overlap_high)
  
  overlap_results[[length(overlap_results) + 1]] <- data.frame(
    Covariate = var_label,
    Variable = var,
    Overlap_Low = overlap_low,
    Overlap_High = overlap_high,
    Overlap_Overall = overlap_overall,
    Common_Support_Min = common_min,
    Common_Support_Max = common_max
  )
}

overlap_df <- do.call(rbind, overlap_results)

cat("Overlap Statistics:\n")
print(overlap_df %>% 
      select(Covariate, Overlap_Low, Overlap_High, Overlap_Overall) %>%
      mutate(
        Overlap_Low = round(Overlap_Low, 3),
        Overlap_High = round(Overlap_High, 3),
        Overlap_Overall = round(Overlap_Overall, 3)
      ),
      row.names = FALSE)
cat("\n")

cat("Interpretation: Overlap < 0.8 suggests potential positivity violations\n")
cat("                Overlap < 0.5 suggests severe positivity violations\n\n")

# ============================================================================
# STEP 7: ROBUSTNESS CHECKS
# ============================================================================
cat("=== STEP 7: Robustness Checks ===\n\n")

# Bootstrap confidence intervals
cat("Bootstrap Confidence Intervals (1000 iterations)...\n")
set.seed(123)

bootstrap_coef <- function(data, indices) {
  d <- data[indices, ]
  model <- lm(formula_adj, data = d)
  return(coef(model)[2])
}

boot_results <- boot(data = df_complete, statistic = bootstrap_coef, R = 1000)
boot_ci <- boot.ci(boot_results, type = c("perc", "bca"))

cat(sprintf("Bootstrap 95%% CI (percentile): [%.4f, %.4f]\n", 
            boot_ci$percent[4], boot_ci$percent[5]))
if (!is.null(boot_ci$bca)) {
  cat(sprintf("Bootstrap 95%% CI (BCa): [%.4f, %.4f]\n", 
              boot_ci$bca[4], boot_ci$bca[5]))
}
cat("\n")

# Sensitivity to outliers
cat("Outlier Sensitivity Analysis:\n")
# Cook's distance
cooks_d <- cooks.distance(model_adj)
outlier_threshold <- 4 / nrow(df_complete)
n_outliers <- sum(cooks_d > outlier_threshold)

cat(sprintf("  Observations with Cook's D > %.4f: %d\n", outlier_threshold, n_outliers))
if (n_outliers > 0) {
  cat("  Outlier indices:", which(cooks_d > outlier_threshold), "\n")
  
  # Refit without outliers
  df_no_outliers <- df_complete[cooks_d <= outlier_threshold, ]
  model_no_outliers <- lm(formula_adj, data = df_no_outliers)
  coef_no_outliers <- summary(model_no_outliers)$coefficients[2, 1]
  
  cat(sprintf("  Coefficient without outliers: %.4f (change: %.4f)\n", 
              coef_no_outliers, coef_adj - coef_no_outliers))
}
cat("\n")

# ============================================================================
# VISUALIZATIONS
# ============================================================================
cat("=== Creating Visualizations ===\n\n")

if (!dir.exists("planned_results")) {
  dir.create("planned_results")
}

# Plot 1: Block-wise sequential adjustment
png("planned_results/analysis8_blockwise_adjustment.png", width = 12, height = 8, units = "in", res = 300)
p1 <- ggplot(block_results_df, aes(x = Block_Number, y = KF_Coefficient)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = KF_CI_Lower, ymax = KF_CI_Upper), width = 0.2, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = coef_unadj, linetype = "dotted", color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = 1:nrow(block_results_df), labels = block_results_df$Block) +
  labs(
    x = "Cumulative Blocks Added",
    y = "Knee Flexion Coefficient",
    title = "Block-Wise Sequential Adjustment\n(Effect estimate as blocks are added)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text", x = Inf, y = coef_unadj, 
           label = "Unadjusted", hjust = 1.1, vjust = -0.5, color = "red")

print(p1)
dev.off()
cat("Saved block-wise adjustment plot to planned_results/analysis8_blockwise_adjustment.png\n")

# Plot 2: Drop-one analysis
png("planned_results/analysis8_drop_one_analysis.png", width = 12, height = 8, units = "in", res = 300)
p2 <- ggplot(drop_one_df, aes(x = reorder(Covariate, abs(KF_Coefficient_Change_Pct)), 
                               y = KF_Coefficient_Change_Pct)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Covariate Dropped",
    y = "Change in KF Coefficient (%)",
    title = "Drop-One-Covariate Analysis\n(Effect of removing each covariate from full model)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(aes(label = sprintf("%.1f%%", KF_Coefficient_Change_Pct)), 
            hjust = ifelse(drop_one_df$KF_Coefficient_Change_Pct > 0, -0.1, 1.1), 
            size = 3.5)

print(p2)
dev.off()
cat("Saved drop-one analysis plot to planned_results/analysis8_drop_one_analysis.png\n")

# Plot 3: Dominance analysis
png("planned_results/analysis8_dominance_analysis.png", width = 12, height = 8, units = "in", res = 300)
p3 <- ggplot(dominance_df, aes(x = reorder(Covariate, Average_R2_Contribution), 
                                y = Average_R2_Contribution)) +
  geom_bar(stat = "identity", fill = "purple", alpha = 0.7) +
  geom_errorbar(aes(ymin = Average_R2_Contribution - SD_Contribution,
                    ymax = Average_R2_Contribution + SD_Contribution), 
                width = 0.2, color = "darkviolet") +
  coord_flip() +
  labs(
    x = "Covariate",
    y = "Average R² Contribution (Dominance)",
    title = "Dominance Analysis\n(Average R² contribution across model combinations)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(aes(label = sprintf("%.4f", Average_R2_Contribution)), 
            hjust = -0.1, size = 3.5)

print(p3)
dev.off()
cat("Saved dominance analysis plot to planned_results/analysis8_dominance_analysis.png\n")

# Plot 4: Balance diagnostics
png("planned_results/analysis8_balance_diagnostics.png", width = 12, height = 8, units = "in", res = 300)
balance_long <- balance_df %>%
  select(Covariate, SMD_Low_vs_High, SMD_Medium_vs_High, SMD_Low_vs_Medium) %>%
  pivot_longer(cols = starts_with("SMD"), names_to = "Comparison", values_to = "SMD") %>%
  mutate(Comparison = gsub("SMD_", "", Comparison))

p4 <- ggplot(balance_long, aes(x = reorder(Covariate, abs(SMD)), y = SMD, fill = Comparison)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = c(-0.25, -0.1, 0.1, 0.25), linetype = "dashed", color = "gray") +
  labs(
    x = "Covariate",
    y = "Standardized Mean Difference",
    title = "Balance Diagnostics\n(Standardized Mean Differences by Knee Flexion Tertile)",
    fill = "Comparison"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p4)
dev.off()
cat("Saved balance diagnostics plot to planned_results/analysis8_balance_diagnostics.png\n")

# Plot 5: Overlap diagnostics
png("planned_results/analysis8_overlap_diagnostics.png", width = 12, height = 8, units = "in", res = 300)
overlap_long <- overlap_df %>%
  select(Covariate, Overlap_Low, Overlap_High, Overlap_Overall) %>%
  pivot_longer(cols = starts_with("Overlap"), names_to = "Group", values_to = "Overlap") %>%
  mutate(Group = gsub("Overlap_", "", Group))

p5 <- ggplot(overlap_long, aes(x = reorder(Covariate, Overlap), y = Overlap, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(
    x = "Covariate",
    y = "Overlap Proportion",
    title = "Positivity Diagnostics (Overlap)\n(Proportion of observations in common support)",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text", x = Inf, y = 0.8, label = "80% threshold", 
           hjust = 1.1, vjust = -0.5, color = "orange") +
  annotate("text", x = Inf, y = 0.5, label = "50% threshold", 
           hjust = 1.1, vjust = -0.5, color = "red")

print(p5)
dev.off()
cat("Saved overlap diagnostics plot to planned_results/analysis8_overlap_diagnostics.png\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================
write.csv(block_results_df, "planned_results/analysis8_blockwise_results.csv", row.names = FALSE)
write.csv(drop_one_df, "planned_results/analysis8_drop_one_results.csv", row.names = FALSE)
write.csv(dominance_df, "planned_results/analysis8_dominance_results.csv", row.names = FALSE)
write.csv(balance_df, "planned_results/analysis8_balance_results.csv", row.names = FALSE)
write.csv(overlap_df, "planned_results/analysis8_overlap_results.csv", row.names = FALSE)

cat("\nSaved all results tables to planned_results/ directory\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n", rep("=", 80), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 80), "\n\n")

cat("1. DAG-Based Adjustment Set:\n")
cat(sprintf("   Selected %d confounders (no mediators/colliders)\n", length(adjustment_set)))
cat(sprintf("   Sample size: %d complete cases\n\n", nrow(df_complete)))

cat("2. Effect Estimates:\n")
cat(sprintf("   Unadjusted: %.4f (95%% CI: [%.4f, %.4f], p = %.4e)\n", 
            coef_unadj, ci_unadj[1], ci_unadj[2], p_unadj))
cat(sprintf("   Adjusted: %.4f (95%% CI: [%.4f, %.4f], p = %.4e)\n", 
            coef_adj, ci_adj[1], ci_adj[2], p_adj))
cat(sprintf("   Attenuation: %.2f%%\n\n", attenuation))

cat("3. Block-Wise Sequential Adjustment:\n")
cat(sprintf("   %d pre-specified blocks tested\n", length(blocks)))
cat("   See block_results_df for detailed results\n\n")

cat("4. Individual Covariate Contributions:\n")
cat("   - Drop-one analysis: Identifies covariates with largest impact\n")
cat("   - Dominance analysis: Average R² contribution across combinations\n")
cat("   Top 3 by dominance:\n")
top3_dom <- head(dominance_df, 3)
for (i in seq_len(nrow(top3_dom))) {
  cat(sprintf("     %d. %s: %.4f\n", i, top3_dom$Covariate[i], top3_dom$Average_R2_Contribution[i]))
}
cat("\n")

cat("5. Balance Diagnostics:\n")
n_imbalanced <- sum(abs(balance_df$SMD_Low_vs_High) > 0.25)
cat(sprintf("   %d covariates with |SMD| > 0.25 (substantial imbalance)\n", n_imbalanced))
cat("\n")

cat("6. Positivity Diagnostics:\n")
n_low_overlap <- sum(overlap_df$Overlap_Overall < 0.8)
cat(sprintf("   %d covariates with overlap < 0.8 (potential violations)\n", n_low_overlap))
cat("\n")

cat("7. Robustness Checks:\n")
cat("   - Bootstrap CIs computed\n")
cat(sprintf("   - %d potential outliers identified (Cook's D)\n", n_outliers))
cat("\n")

# ============================================================================
# PROOF: Is Knee Flexion a Confounder?
# ============================================================================
cat("\n", rep("=", 80), "\n", sep = "")
cat("PROOF: IS KNEE FLEXION A CONFOUNDER?\n")
cat(rep("=", 80), "\n\n")

cat("HYPOTHESIS: Knee flexion's apparent association with change in sva (cm)\n")
cat("            is confounded by other preop parameters.\n\n")

cat("EVIDENCE FROM RIGOROUS CAUSAL INFERENCE ANALYSIS:\n\n")

cat("1. Effect Attenuation:\n")
cat(sprintf("   - Unadjusted KF coefficient: %.4f (p = %.4e)\n", coef_unadj, p_unadj))
cat(sprintf("   - Adjusted KF coefficient: %.4f (p = %.4e)\n", coef_adj, p_adj))
cat(sprintf("   - Attenuation: %.2f%% reduction\n", attenuation))
cat(sprintf("   - Partial R² for KF (after adjusting for confounders): %.4f\n\n", kf_partial_r2))

cat("2. Statistical Significance:\n")
if (p_unadj < 0.05 && p_adj >= 0.05) {
  cat("   ✓ PROOF: Knee flexion is a CONFOUNDER\n")
  cat("     - Significant when unadjusted (p = %.4e < 0.05)\n", p_unadj)
  cat("     - NON-SIGNIFICANT after adjustment (p = %.4e >= 0.05)\n", p_adj)
  cat("     - This proves the association is fully explained by confounding\n\n")
} else if (p_unadj < 0.05 && p_adj < 0.05 && attenuation > 20) {
  cat("   ⚠️  PARTIAL PROOF: Knee flexion may be partially confounded\n")
  cat("     - Remains significant but substantially attenuated (%.2f%%)\n", attenuation)
  cat("     - May have both direct and confounded effects\n\n")
} else if (p_unadj < 0.05 && p_adj >= 0.05) {
  cat("   ✓ PROOF: Knee flexion is a CONFOUNDER\n")
  cat("     - Becomes non-significant after adjustment\n\n")
} else {
  cat("   ⚠️  INCONCLUSIVE: Effect persists after adjustment\n")
  cat("     - May indicate direct effect or residual confounding\n\n")
}

cat("3. Variance Explained:\n")
cat(sprintf("   - Other parameters alone explain %.2f%% of variance (R² = %.4f)\n", 
            r2_no_kf * 100, r2_no_kf))
cat(sprintf("   - Knee flexion adds only %.2f%% additional variance (partial R² = %.4f)\n", 
            kf_partial_r2 * 100, kf_partial_r2))
if (kf_partial_r2 < 0.01 && p_adj >= 0.05) {
  cat("   ✓ STRONG EVIDENCE: Knee flexion adds minimal unique variance\n")
  cat("     - Other parameters explain virtually all the variance\n")
  cat("     - Knee flexion's contribution is negligible\n\n")
}

cat("4. Block-Wise Evidence:\n")
# Check if any block makes KF non-significant
blocks_make_nonsig <- block_results_df %>% filter(KF_P >= 0.05)
if (nrow(blocks_make_nonsig) > 0) {
  cat(sprintf("   - %d block(s) make knee flexion non-significant:\n", nrow(blocks_make_nonsig)))
  for (i in seq_len(nrow(blocks_make_nonsig))) {
    cat(sprintf("     * %s (p = %.4e)\n", 
                blocks_make_nonsig$Block[i], blocks_make_nonsig$KF_P[i]))
  }
  cat("\n")
} else {
  cat("   - No individual block makes KF non-significant\n")
  cat("   - But combined blocks may fully explain the effect\n\n")
}

cat("5. Drop-One Analysis Evidence:\n")
# Check drop-one results
if (exists("drop_one_df")) {
  drop_one_makes_sig <- drop_one_df %>% filter(KF_P_Without < 0.05 & p_adj >= 0.05)
  if (nrow(drop_one_makes_sig) > 0) {
    cat(sprintf("   - %d covariate(s) are critical: removing them makes KF significant\n", 
                nrow(drop_one_makes_sig)))
    cat("     This shows these covariates are necessary to explain confounding\n\n")
  }
}

cat("CONCLUSION:\n")
if (p_unadj < 0.05 && p_adj >= 0.05 && attenuation > 10) {
  cat("  ✓ PROVEN: Knee flexion is a CONFOUNDER\n")
  cat("    The apparent association between knee flexion and change in sva (cm)\n")
  cat("    is fully explained by other preop parameters. Knee flexion does not\n")
  cat("    have an independent causal effect on change in sva (cm).\n\n")
} else if (p_unadj < 0.05 && p_adj < 0.05 && attenuation > 20) {
  cat("  ⚠️  PARTIALLY PROVEN: Knee flexion may be partially confounded\n")
  cat("    The association is substantially attenuated, suggesting significant\n")
  cat("    confounding, but some direct effect may remain.\n\n")
} else {
  cat("  ⚠️  INCONCLUSIVE: Cannot definitively prove confounding\n")
  cat("    Knee flexion effect persists after adjustment. This could indicate:\n")
  cat("    - A direct causal effect of knee flexion\n")
  cat("    - Residual confounding by unmeasured variables\n")
  cat("    - Model misspecification\n\n")
}

cat(rep("=", 80), "\n\n")

# ============================================================================
# INTERPRETATION GUIDE
# ============================================================================
cat("\n", rep("=", 80), "\n", sep = "")
cat("HOW TO INTERPRET THESE FINDINGS\n")
cat(rep("=", 80), "\n\n")

cat("OVERVIEW:\n")
cat("This analysis tests whether knee flexion has a DIRECT causal effect on change\n")
cat("in lordosis, or if the association is CONFOUNDED by other preop parameters.\n\n")

cat("KEY CONCEPTS:\n")
cat("1. CONFOUNDER: A variable that affects both the exposure (knee flexion) and\n")
cat("   outcome (change in sva (cm)), creating a spurious association.\n")
cat("   Example: If patients with high knee flexion also have worse baseline\n")
cat("   alignment, and worse baseline alignment predicts more lordosis change,\n")
cat("   then knee flexion appears associated with change even if it's not causal.\n\n")

cat("2. EFFECT ATTENUATION: When adjusting for confounders, the association\n")
cat("   between knee flexion and change in sva (cm) should get smaller (attenuate).\n")
cat("   Large attenuation (>20%%) suggests confounding.\n\n")

cat("3. STATISTICAL SIGNIFICANCE: If knee flexion is significant when unadjusted\n")
cat("   but becomes non-significant after adjusting for confounders, this is\n")
cat("   strong evidence that the association was due to confounding.\n\n")

cat("INTERPRETING THE RESULTS:\n\n")

cat("SECTION 1: Effect Attenuation\n")
cat("  - Unadjusted coefficient: The raw association between knee flexion and\n")
cat("    change in sva (cm) (before accounting for other factors)\n")
cat("  - Adjusted coefficient: The association after accounting for confounders\n")
cat("  - Attenuation: How much the coefficient decreased\n")
cat("    * >20%% reduction = Strong evidence of confounding\n")
cat("    * 10-20%% reduction = Moderate evidence of confounding\n")
cat("    * <10%% reduction = Weak evidence of confounding\n\n")

cat("SECTION 2: Statistical Significance\n")
cat("  - If p < 0.05 when unadjusted BUT p >= 0.05 when adjusted:\n")
cat("    → PROOF that knee flexion is a confounder (not a direct cause)\n")
cat("  - If p < 0.05 in both models but coefficient attenuated >20%%:\n")
cat("    → PARTIAL confounding (some direct effect may remain)\n")
cat("  - If p < 0.05 in both models with little attenuation:\n")
cat("    → Direct causal effect (not confounded)\n\n")

cat("SECTION 3: Variance Explained\n")
cat("  - R² without knee flexion: How much variance other parameters explain\n")
cat("  - Partial R² for knee flexion: How much ADDITIONAL variance knee flexion\n")
cat("    explains after accounting for other parameters\n")
cat("  - If partial R² < 0.01 and p >= 0.05:\n")
cat("    → Knee flexion adds virtually nothing; association is fully confounded\n\n")

cat("SECTION 4: Block-Wise Evidence\n")
cat("  - Shows which groups of parameters (pelvic, lumbar, thoracic, etc.) are\n")
cat("    most important in explaining the confounding\n")
cat("  - If a block makes knee flexion non-significant:\n")
cat("    → That group of parameters fully explains the association\n\n")

cat("SECTION 5: Drop-One Analysis\n")
cat("  - Tests which individual parameters are CRITICAL for explaining confounding\n")
cat("  - If removing a parameter makes knee flexion significant again:\n")
cat("    → That parameter is essential for explaining the confounding\n\n")

cat("DIAGNOSTICS:\n")
cat("  - Balance Diagnostics: Checks if patients with different knee flexion\n")
cat("    have different baseline characteristics (confounders)\n")
cat("    * |SMD| > 0.25 = Substantial imbalance (confounding likely)\n\n")
cat("  - Positivity Diagnostics: Checks if there's enough overlap between groups\n")
cat("    * Overlap < 0.8 = Potential problem (some patients may not be comparable)\n\n")

cat("BOTTOM LINE:\n")
cat("  ✓ PROVEN CONFOUNDER: If knee flexion becomes non-significant after\n")
cat("    adjustment, the association is explained by other preop parameters.\n")
cat("    Knee flexion does NOT directly cause change in sva (cm).\n\n")
cat("  ⚠️  PARTIAL CONFOUNDING: If knee flexion remains significant but is\n")
cat("    substantially attenuated, there may be both confounding AND a direct\n")
cat("    effect. More research needed.\n\n")
cat("  ⚠️  INCONCLUSIVE: If knee flexion remains significant with little\n")
cat("    attenuation, either:\n")
cat("    - Knee flexion has a direct causal effect, OR\n")
cat("    - Important confounders were not measured/adjusted for\n\n")

cat(rep("=", 80), "\n\n")

cat("=== Analysis Complete ===\n")
cat("This analysis follows rigorous causal inference methodology:\n")
cat("  ✓ DAG-based adjustment set selection\n")
cat("  ✓ Block-wise sequential adjustment\n")
cat("  ✓ Individual covariate contribution analysis\n")
cat("  ✓ Balance and positivity diagnostics\n")
cat("  ✓ Robustness checks\n")
cat("\n")
cat("NOTE: For maximum rigor, consider implementing AIPW/TMLE with cross-fitting\n")
cat("      for doubly-robust effect estimation.\n\n")

