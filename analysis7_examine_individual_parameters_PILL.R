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
df <- attach_planned_dll(df, db_path)

# Filter out low-volume surgeons if enabled
if (EXCLUDE_LOW_VOLUME_SURGEONS) {
  df <- filter_low_volume_surgeons(df, min_surgeon_cases = 10)
}

df <- restrict_planned_dll_analysis_cohort(df)
cat(sprintf(
  "\nCohort restricted to |preop PI–LL| > %d° and preop SVA C2–C7 < %d (non-missing): n = %d rows\n",
  PREOP_ABS_PI_LL_GT,
  PREOP_SVA_C2_C7_LT,
  nrow(df)
))

# Outcome: Planned \u0394LL (utils/planned_pi_ll.R)

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

# Drop patients with missing data for analysis
# Base variables: knee flexion and Planned \u0394LL
# We'll check each parameter individually for missingness
df_base <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle) & !is.na(planned_DLL) & planned_DLL >= PLANNED_DLL_MIN_KEEP
  )

cat(paste("\n=== Individual Parameter Variance Analysis ===\n"))
cat(paste("Base sample size:", nrow(df_base), "patients with knee flexion and Planned \u0394LL\n\n"))

# ============================================================================
# Define parameters to test
# ============================================================================
parameters <- list(
  list(var = "LATpre_S1PI", label = "S1PI (preop)", type = "pelvic"),
  list(var = "LATpre_L1_S1", label = "Lordosis L1-S1", type = "lumbar"),
  list(var = "LATpre_L4_S1", label = "Lordosis L4-S1", type = "lumbar"),
  list(var = "LATpre_T2_T12", label = "Thoracic Kyphosis T2-T12", type = "thoracic"),
  list(var = "LATpre_S1PT", label = "Pelvic Tilt S1PT", type = "pelvic"),
  list(var = "LATpre_SVA_C2_S1", label = "SVA", type = "sagittal"),
  list(var = "LATpre_T4_L1_PA", label = "T4-L1 PA", type = "sagittal"),
  list(var = "preop_PI_LL", label = "PI-LL Mismatch", type = "pelvic")
)

# Add age if available
if (!is.null(age_var)) {
  parameters <- c(parameters, list(list(var = age_var, label = "Age", type = "demographic")))
}

cat("Parameters to be tested:\n")
for (i in seq_along(parameters)) {
  cat(sprintf("  %d. %s (%s)\n", i, parameters[[i]]$label, parameters[[i]]$var))
}
cat("\n")

# ============================================================================
# Base model: Knee flexion only
# ============================================================================
cat("=== Base Model: Knee Flexion Only ===\n")
model_base <- lm(planned_DLL ~ LATpre_LL_KneeAngle, data = df_base)
summary_base <- summary(model_base)

base_r2 <- summary_base$r.squared
base_kf_coef <- summary_base$coefficients[2, 1]
base_kf_p <- summary_base$coefficients[2, 4]

cat(sprintf("R²: %.4f\n", base_r2))
cat(sprintf("Knee Flexion coefficient: %.4f (p = %.4e)\n\n", base_kf_coef, base_kf_p))

# ============================================================================
# Test each parameter individually
# ============================================================================
cat("=== Testing Each Parameter Individually ===\n\n")

results_list <- list()

for (i in seq_along(parameters)) {
  param <- parameters[[i]]
  param_var <- param$var
  param_label <- param$label
  
  cat(sprintf("Testing %d/%d: %s (%s)\n", i, length(parameters), param_label, param_var))
  
  # Filter to patients with non-missing data for this parameter
  df_param <- df_base %>%
    filter(!is.na(!!sym(param_var)))
  
  n_available <- nrow(df_param)
  n_dropped <- nrow(df_base) - n_available
  
  if (n_available < 10) {
    cat(sprintf("  ⚠️  Skipping: Only %d patients have data for this parameter (need at least 10)\n\n", n_available))
    next
  }
  
  if (n_dropped > 0) {
    cat(sprintf("  Using %d patients (dropped %d with missing %s)\n", n_available, n_dropped, param_var))
  }
  
  # Base model on this subset
  model_base_subset <- lm(planned_DLL ~ LATpre_LL_KneeAngle, data = df_param)
  summary_base_subset <- summary(model_base_subset)
  base_r2_subset <- summary_base_subset$r.squared
  base_kf_coef_subset <- summary_base_subset$coefficients[2, 1]
  base_kf_p_subset <- summary_base_subset$coefficients[2, 4]
  
  # Model with parameter added
  formula_with_param <- as.formula(paste("planned_DLL ~ LATpre_LL_KneeAngle +", param_var))
  model_with_param <- lm(formula_with_param, data = df_param)
  summary_with_param <- summary(model_with_param)
  
  # Extract results
  r2_with_param <- summary_with_param$r.squared
  r2_change <- r2_with_param - base_r2_subset  # R² change (partial R² for this parameter)
  
  # Get knee flexion coefficient and p-value from model with parameter
  kf_coef_with_param <- summary_with_param$coefficients[2, 1]
  kf_p_with_param <- summary_with_param$coefficients[2, 4]
  
  # Get parameter coefficient and p-value
  param_coef <- summary_with_param$coefficients[param_var, 1]
  param_p <- summary_with_param$coefficients[param_var, 4]
  
  # Calculate change in knee flexion coefficient
  kf_coef_change <- kf_coef_with_param - base_kf_coef_subset
  kf_coef_change_pct <- (kf_coef_change / base_kf_coef_subset) * 100
  
  # Calculate partial correlation for the parameter
  # Partial r = t / sqrt(t² + df)
  t_stat_param <- summary_with_param$coefficients[param_var, 3]
  df_residual <- summary_with_param$df[2]
  partial_r_param <- t_stat_param / sqrt(t_stat_param^2 + df_residual)
  
  # Calculate percentage of TOTAL variance explained (not relative to base R²)
  # This is the R² change as a percentage of total variance (can't exceed 100%)
  # If base R² is very small, the old calculation could exceed 100%
  r2_change_pct_of_total <- r2_change * 100  # As percentage of total variance
  r2_change_relative_to_base <- ifelse(base_r2_subset > 0, 
                                       (r2_change / base_r2_subset) * 100, 
                                       NA)  # Relative to base (can exceed 100%)
  
  # Store results
  results_list[[length(results_list) + 1]] <- data.frame(
    Parameter = param_label,
    Variable = param_var,
    Type = param$type,
    N = n_available,
    Base_R2 = base_r2_subset,
    R2_With_Param = r2_with_param,
    R2_Change = r2_change,
    R2_Change_Pct_Of_Total = r2_change_pct_of_total,  # % of total variance (0-100%)
    R2_Change_Relative_To_Base = r2_change_relative_to_base,  # Relative to base (can >100%)
    Partial_R2 = r2_change,
    Partial_R = partial_r_param,
    Param_Coefficient = param_coef,
    Param_P = param_p,
    KF_Coefficient_Base = base_kf_coef_subset,
    KF_Coefficient_With_Param = kf_coef_with_param,
    KF_Coefficient_Change = kf_coef_change,
    KF_Coefficient_Change_Pct = kf_coef_change_pct,
    KF_P_Base = base_kf_p_subset,
    KF_P_With_Param = kf_p_with_param,
    KF_Significant_Base = base_kf_p_subset < 0.05,
    KF_Significant_With_Param = kf_p_with_param < 0.05
  )
  
  cat(sprintf("  Base R²: %.4f\n", base_r2_subset))
  cat(sprintf("  R² with %s: %.4f\n", param_label, r2_with_param))
  cat(sprintf("  R² change (partial R²): %.4f (%.2f%% of total variance)\n", 
              r2_change, r2_change_pct_of_total))
  if (!is.na(r2_change_relative_to_base)) {
    cat(sprintf("  R² change relative to base: %.2f%% (%.1fx the base R²)\n", 
                r2_change_relative_to_base, r2_change / base_r2_subset))
  }
  cat(sprintf("  %s coefficient: %.4f (p = %.4e)\n", param_label, param_coef, param_p))
  cat(sprintf("  Knee Flexion coefficient change: %.4f (%.2f%% change)\n", 
              kf_coef_change, kf_coef_change_pct))
  cat(sprintf("  Knee Flexion p-value: %.4e -> %.4e\n", base_kf_p_subset, kf_p_with_param))
  if (base_kf_p_subset < 0.05 && kf_p_with_param >= 0.05) {
    cat(sprintf("  ✓ %s makes knee flexion non-significant!\n", param_label))
  }
  cat("\n")
}

# Combine results
if (length(results_list) == 0) {
  cat("ERROR: No parameters could be tested. Check data availability.\n")
  quit(status = 1)
}

results_df <- do.call(rbind, results_list)

# Sort by R² change (descending)
results_df <- results_df %>%
  arrange(desc(R2_Change))

cat("=== Summary: Parameters Ranked by R² Contribution ===\n\n")
cat("NOTE: R² Change % shows percentage of TOTAL variance explained (0-100%%)\n")
cat("      Values >100%% in old calculation meant the parameter explained more\n")
cat("      variance than knee flexion alone (which is possible if base R² is small)\n\n")
print(results_df %>% 
      select(Parameter, N, R2_Change, R2_Change_Pct_Of_Total, R2_Change_Relative_To_Base, 
             Partial_R, Param_P, KF_Coefficient_Change_Pct, KF_P_With_Param) %>%
      mutate(
        R2_Change = round(R2_Change, 4),
        R2_Change_Pct_Of_Total = round(R2_Change_Pct_Of_Total, 2),
        R2_Change_Relative_To_Base = round(R2_Change_Relative_To_Base, 1),
        Partial_R = round(Partial_R, 3),
        Param_P = formatC(Param_P, format = "e", digits = 2),
        KF_Coefficient_Change_Pct = round(KF_Coefficient_Change_Pct, 2),
        KF_P_With_Param = formatC(KF_P_With_Param, format = "e", digits = 2)
      ) %>%
      rename(
        "R² Change" = R2_Change,
        "R² Change (% of Total)" = R2_Change_Pct_Of_Total,
        "R² Change (Relative to Base)" = R2_Change_Relative_To_Base
      ),
      row.names = FALSE)

cat("\n")

# ============================================================================
# Create visualizations
# ============================================================================
cat("=== Creating Visualizations ===\n\n")

pr_dir <- ensure_planned_results_dir()

# Plot 1: R² change (partial R²) for each parameter
png(file.path(pr_dir, "analysis7_r2_contribution.png"), width = 12, height = 8, units = "in", res = 300)
p1 <- ggplot(results_df, aes(x = reorder(Parameter, R2_Change), y = R2_Change)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(
    x = "Parameter",
    y = "Partial R² (R² Change)",
    title = "Variance Explained by Each Parameter\n(Partial R² when added to knee flexion model)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(aes(label = sprintf("%.4f", R2_Change)), 
            hjust = -0.1, size = 3.5)

print(p1)
dev.off()
cat("Saved R² contribution plot to planned_results/analysis7_r2_contribution.png\n")

# Plot 2: R² change as percentage of total variance (0-100%)
png(file.path(pr_dir, "analysis7_r2_contribution_pct.png"), width = 12, height = 8, units = "in", res = 300)
p2 <- ggplot(results_df, aes(x = reorder(Parameter, R2_Change_Pct_Of_Total), y = R2_Change_Pct_Of_Total)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  coord_flip() +
  labs(
    x = "Parameter",
    y = "R² Change (% of Total Variance)",
    title = "Variance Explained by Each Parameter\n(% of total variance explained, 0-100%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(aes(label = sprintf("%.1f%%", R2_Change_Pct_Of_Total)), 
            hjust = -0.1, size = 3.5)

print(p2)
dev.off()
cat("Saved R² contribution percentage plot to planned_results/analysis7_r2_contribution_pct.png\n")

# Plot 3: Effect on knee flexion coefficient
png(file.path(pr_dir, "analysis7_kf_coefficient_change.png"), width = 12, height = 8, units = "in", res = 300)
p3 <- ggplot(results_df, aes(x = reorder(Parameter, abs(KF_Coefficient_Change_Pct)), 
                             y = KF_Coefficient_Change_Pct)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Parameter",
    y = "Change in Knee Flexion Coefficient (%)",
    title = "Effect of Each Parameter on Knee Flexion Coefficient\n(% change when parameter is added to model)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(aes(label = sprintf("%.1f%%", KF_Coefficient_Change_Pct)), 
            hjust = ifelse(results_df$KF_Coefficient_Change_Pct > 0, -0.1, 1.1), 
            size = 3.5)

print(p3)
dev.off()
cat("Saved knee flexion coefficient change plot to planned_results/analysis7_kf_coefficient_change.png\n")

# Plot 4: Partial correlation coefficients
png(file.path(pr_dir, "analysis7_partial_correlations.png"), width = 12, height = 8, units = "in", res = 300)
p4 <- ggplot(results_df, aes(x = reorder(Parameter, abs(Partial_R)), y = Partial_R)) +
  geom_bar(stat = "identity", fill = "purple", alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Parameter",
    y = "Partial Correlation (r)",
    title = "Partial Correlations\n(Relationship with change in pi-ll mismatch, controlling for knee flexion)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(aes(label = sprintf("%.3f", Partial_R)), 
            hjust = ifelse(results_df$Partial_R > 0, -0.1, 1.1), 
            size = 3.5)

print(p4)
dev.off()
cat("Saved partial correlations plot to planned_results/analysis7_partial_correlations.png\n")

# Plot 5: Combined visualization - R² change vs effect on KF coefficient
png(file.path(pr_dir, "analysis7_combined_importance.png"), width = 12, height = 8, units = "in", res = 300)
p5 <- ggplot(results_df, aes(x = R2_Change, y = abs(KF_Coefficient_Change_Pct), 
                             label = Parameter, color = Type)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text(hjust = 0, vjust = 0, nudge_x = 0.001, nudge_y = 0.5, size = 3) +
  labs(
    x = "Partial R² (Variance Explained)",
    y = "Absolute % Change in Knee Flexion Coefficient",
    title = "Parameter Importance: Variance Explained vs Effect on Knee Flexion",
    color = "Parameter Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p5)
dev.off()
cat("Saved combined importance plot to planned_results/analysis7_combined_importance.png\n")

# Plot 6: Knee flexion p-value change
png(file.path(pr_dir, "analysis7_kf_pvalue_change.png"), width = 12, height = 8, units = "in", res = 300)
results_df$KF_P_Change <- results_df$KF_P_With_Param - results_df$KF_P_Base
results_df$KF_Significance_Change <- factor(
  ifelse(results_df$KF_Significant_Base & !results_df$KF_Significant_With_Param, "Made Non-Sig",
         ifelse(!results_df$KF_Significant_Base & results_df$KF_Significant_With_Param, "Made Sig",
                ifelse(results_df$KF_Significant_Base & results_df$KF_Significant_With_Param, "Remained Sig",
                       "Remained Non-Sig"))),
  levels = c("Made Non-Sig", "Remained Sig", "Remained Non-Sig", "Made Sig")
)

p6 <- ggplot(results_df, aes(x = reorder(Parameter, KF_P_Change), y = -log10(KF_P_With_Param))) +
  geom_bar(stat = "identity", aes(fill = KF_Significance_Change), alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "Parameter",
    y = "-log10(p-value) for Knee Flexion",
    title = "Effect of Each Parameter on Knee Flexion Significance",
    fill = "Significance\nChange"
  ) +
  scale_fill_manual(values = c("Made Non-Sig" = "green", 
                                "Remained Sig" = "orange",
                                "Remained Non-Sig" = "gray",
                                "Made Sig" = "blue")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text", x = Inf, y = -log10(0.05), 
           label = "p = 0.05 threshold", hjust = 1.1, vjust = -0.5, color = "red")

print(p6)
dev.off()
cat("Saved knee flexion p-value change plot to planned_results/analysis7_kf_pvalue_change.png\n")

# ============================================================================
# Summary statistics
# ============================================================================
cat("\n=== Summary Statistics ===\n\n")

cat("Top 3 parameters by R² contribution:\n")
top3_r2 <- head(results_df %>% arrange(desc(R2_Change)), 3)
for (i in seq_len(nrow(top3_r2))) {
  cat(sprintf("  %d. %s: R² change = %.4f (%.2f%% of total variance)\n", 
              i, top3_r2$Parameter[i], top3_r2$R2_Change[i], top3_r2$R2_Change_Pct_Of_Total[i]))
  if (!is.na(top3_r2$R2_Change_Relative_To_Base[i])) {
    cat(sprintf("      This is %.1fx the base R² (%.2f%% relative to base)\n", 
                top3_r2$R2_Change[i] / top3_r2$Base_R2[i], 
                top3_r2$R2_Change_Relative_To_Base[i]))
  }
}
cat("\n")

cat("Parameters that make knee flexion non-significant (p > 0.05):\n")
non_sig_params <- results_df %>% filter(KF_P_With_Param >= 0.05 & KF_Significant_Base)
if (nrow(non_sig_params) > 0) {
  for (i in seq_len(nrow(non_sig_params))) {
    cat(sprintf("  - %s: p = %.4e, R² change = %.4f\n", 
                non_sig_params$Parameter[i], 
                non_sig_params$KF_P_With_Param[i],
                non_sig_params$R2_Change[i]))
  }
} else {
  cat("  None\n")
}
cat("\n")

cat("Parameters with largest effect on knee flexion coefficient (>10% change):\n")
large_effect <- results_df %>% filter(abs(KF_Coefficient_Change_Pct) > 10)
if (nrow(large_effect) > 0) {
  for (i in seq_len(nrow(large_effect))) {
    cat(sprintf("  - %s: %.2f%% change, R² change = %.4f\n", 
                large_effect$Parameter[i], 
                large_effect$KF_Coefficient_Change_Pct[i],
                large_effect$R2_Change[i]))
  }
} else {
  cat("  None\n")
}
cat("\n")

# Save results table
write.csv(results_df, file.path(pr_dir, "analysis7_individual_parameter_results.csv"), row.names = FALSE)
cat("Saved detailed results table to planned_results/analysis7_individual_parameter_results.csv\n")

# ============================================================================
# PROOF: Is Knee Flexion a Confounder?
# ============================================================================
cat("\n", rep("=", 80), "\n", sep = "")
cat("PROOF: IS KNEE FLEXION A CONFOUNDER?\n")
cat(rep("=", 80), "\n\n")

cat("HYPOTHESIS: Knee flexion's apparent association with change in pi-ll mismatch\n")
cat("            is confounded by other preop parameters.\n\n")

cat("EVIDENCE:\n\n")

# Test with model that includes ALL parameters (if we have complete data)
# First, identify which parameters we successfully tested
tested_params <- results_df$Variable
tested_labels <- results_df$Parameter

# Find subset with all tested parameters available
all_params_vars <- c("LATpre_LL_KneeAngle", "planned_DLL", tested_params)
df_all_params <- df_base %>%
  filter(complete.cases(select(., all_of(all_params_vars))))

if (nrow(df_all_params) >= 10) {
  cat(sprintf("1. Testing with ALL parameters simultaneously (n = %d):\n", nrow(df_all_params)))
  
  # Unadjusted model
  model_unadj_all <- lm(planned_DLL ~ LATpre_LL_KneeAngle, data = df_all_params)
  summary_unadj_all <- summary(model_unadj_all)
  coef_unadj_all <- summary_unadj_all$coefficients[2, 1]
  p_unadj_all <- summary_unadj_all$coefficients[2, 4]
  r2_unadj_all <- summary_unadj_all$r.squared
  
  # Fully adjusted model
  if (length(tested_params) > 0) {
    formula_adj_all <- as.formula(paste("planned_DLL ~ LATpre_LL_KneeAngle +", 
                                        paste(tested_params, collapse = " + ")))
    model_adj_all <- lm(formula_adj_all, data = df_all_params)
    summary_adj_all <- summary(model_adj_all)
    coef_adj_all <- summary_adj_all$coefficients[2, 1]
    p_adj_all <- summary_adj_all$coefficients[2, 4]
    r2_adj_all <- summary_adj_all$r.squared
    
    # Calculate attenuation
    attenuation_all <- abs((coef_unadj_all - coef_adj_all) / coef_unadj_all) * 100
    r2_explained_by_params <- r2_adj_all - r2_unadj_all
    
    cat(sprintf("   Unadjusted KF effect: %.4f (p = %.4e, R² = %.4f)\n", 
                coef_unadj_all, p_unadj_all, r2_unadj_all))
    cat(sprintf("   Adjusted KF effect: %.4f (p = %.4e, R² = %.4f)\n", 
                coef_adj_all, p_adj_all, r2_adj_all))
    cat(sprintf("   Coefficient attenuation: %.2f%%\n", attenuation_all))
    cat(sprintf("   Variance explained by other parameters: %.4f (%.2f%% of unadjusted R²)\n\n", 
                r2_explained_by_params, (r2_explained_by_params / r2_unadj_all) * 100))
    
    # Conclusion
    if (p_unadj_all < 0.05 && p_adj_all >= 0.05 && attenuation_all > 10) {
      cat("   ✓ PROOF: Knee flexion is a CONFOUNDER\n")
      cat("     - Knee flexion is significant when unadjusted (p < 0.05)\n")
      cat("     - Knee flexion becomes NON-SIGNIFICANT after adjusting for other parameters (p >= 0.05)\n")
      cat("     - Coefficient attenuated by >10% after adjustment\n")
      cat("     - This proves the association is explained by confounding\n\n")
    } else if (p_unadj_all < 0.05 && p_adj_all < 0.05 && attenuation_all > 20) {
      cat("   ⚠️  PARTIAL PROOF: Knee flexion may be partially confounded\n")
      cat("     - Knee flexion remains significant but is substantially attenuated (>20%%)\n")
      cat("     - Other parameters explain substantial variance\n")
      cat("     - Knee flexion may have both direct and confounded effects\n\n")
    } else if (p_unadj_all < 0.05 && p_adj_all >= 0.05) {
      cat("   ✓ PROOF: Knee flexion is a CONFOUNDER\n")
      cat("     - Knee flexion becomes non-significant after adjustment\n")
      cat("     - Association is explained by confounding\n\n")
    } else {
      cat("   ⚠️  INCONCLUSIVE: Knee flexion effect persists after adjustment\n")
      cat("     - May indicate direct effect or residual confounding\n\n")
    }
  }
}

# Evidence from individual parameter tests
cat("2. Evidence from individual parameter tests:\n")
n_makes_nonsig <- sum(results_df$KF_Significant_Base & !results_df$KF_Significant_With_Param, na.rm = TRUE)
n_large_effect <- sum(abs(results_df$KF_Coefficient_Change_Pct) > 20, na.rm = TRUE)

cat(sprintf("   - %d parameter(s) make knee flexion non-significant when added individually\n", n_makes_nonsig))
cat(sprintf("   - %d parameter(s) cause >20%% change in knee flexion coefficient\n", n_large_effect))
cat(sprintf("   - Total variance explained by individual parameters: %.4f\n\n", sum(results_df$R2_Change, na.rm = TRUE)))

if (n_makes_nonsig > 0) {
  cat("   ✓ STRONG EVIDENCE: Individual parameters can fully explain knee flexion effect\n")
  cat("     Parameters that make KF non-significant:\n")
  nonsig_params <- results_df %>% filter(KF_Significant_Base & !KF_Significant_With_Param)
  for (i in seq_len(nrow(nonsig_params))) {
    cat(sprintf("       - %s (p changes from %.4e to %.4e)\n", 
                nonsig_params$Parameter[i], 
                nonsig_params$KF_P_Base[i],
                nonsig_params$KF_P_With_Param[i]))
  }
  cat("\n")
}

cat("3. Summary:\n")
cat("   If knee flexion becomes non-significant after adjusting for other preop parameters,\n")
cat("   this PROVES that knee flexion is a confounder - its apparent association with\n")
cat("   change in pi-ll mismatch is explained by its correlation with other preop parameters\n")
cat("   that are the true predictors of change in pi-ll mismatch.\n\n")

cat(rep("=", 80), "\n\n")

cat("\n=== Analysis Complete ===\n")
cat("This analysis shows how much variance each individual preop parameter explains\n")
cat("in the preop knee flexion vs change in lumbar lordosis relationship.\n")
cat("Check the plots in planned_results/ directory for visualizations.\n\n")

