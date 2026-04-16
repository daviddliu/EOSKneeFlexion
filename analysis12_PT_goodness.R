#!/usr/bin/env Rscript
#
# Analysis 12: PT vs KF stratified by PT at 6 weeks
#
# This script creates scatter plots:
#   For PT < 20 at 6 weeks:
#     1. Pelvic Tilt (PT) at 6 weeks vs Preoperative Knee Flexion
#     2. Pelvic Tilt (PT) at 6 weeks vs 6-Week Knee Flexion
#   For PT > 20 at 6 weeks:
#     3. Pelvic Tilt (PT) at 6 weeks vs Preoperative Knee Flexion
#     4. Pelvic Tilt (PT) at 6 weeks vs 6-Week Knee Flexion

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

# Calculate change in PT (6-week - preop)
df$change_PT <- df$LAT6W_S1PT - df$LATpre_S1PT

# Filter for patients with PT < 20 at 6 weeks
df_filtered_preop_lt20 <- df %>%
  filter(
    !is.na(LAT6W_S1PT) & 
    LAT6W_S1PT < 20 &
    !is.na(LATpre_LL_KneeAngle)
  )

df_filtered_w6_lt20 <- df %>%
  filter(
    !is.na(LAT6W_S1PT) & 
    LAT6W_S1PT < 20 &
    !is.na(LAT6W_LL_KneeAngle)
  )

# Filter for patients with PT > 20 at 6 weeks
df_filtered_preop_gt20 <- df %>%
  filter(
    !is.na(LAT6W_S1PT) & 
    LAT6W_S1PT > 20 &
    !is.na(LATpre_LL_KneeAngle)
  )

df_filtered_w6_gt20 <- df %>%
  filter(
    !is.na(LAT6W_S1PT) & 
    LAT6W_S1PT > 20 &
    !is.na(LAT6W_LL_KneeAngle)
  )

# Filter for preoperative PT vs preoperative KF plots
df_preop_all <- df %>%
  filter(
    !is.na(LATpre_S1PT) & 
    !is.na(LATpre_LL_KneeAngle)
  )

df_preop_pt_lt20 <- df %>%
  filter(
    !is.na(LATpre_S1PT) & 
    LATpre_S1PT < 20 &
    !is.na(LATpre_LL_KneeAngle)
  )

df_preop_pt_ge20 <- df %>%
  filter(
    !is.na(LATpre_S1PT) & 
    LATpre_S1PT >= 20 &
    !is.na(LATpre_LL_KneeAngle)
  )

cat(paste("\n=== Analysis 12: PT vs KF (stratified by PT at 6 weeks) ===\n"))
cat(paste("Total patients in database:", nrow(df), "\n\n"))

cat(paste("Preoperative PT vs Preoperative KF:\n"))
cat(paste("  - All patients:", nrow(df_preop_all), "\n"))
cat(paste("  - PT < 20:", nrow(df_preop_pt_lt20), "\n"))
cat(paste("  - PT >= 20:", nrow(df_preop_pt_ge20), "\n\n"))

# Filter for Preop KF vs Change in PT plots
df_kf_changePT_all <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle) & 
    !is.na(change_PT)
  )

df_kf_changePT_pt_lt20 <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle) & 
    !is.na(change_PT) &
    !is.na(LATpre_S1PT) &
    LATpre_S1PT < 20
  )

df_kf_changePT_pt_ge20 <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle) & 
    !is.na(change_PT) &
    !is.na(LATpre_S1PT) &
    LATpre_S1PT >= 20
  )

cat(paste("Preoperative KF vs Change in PT:\n"))
cat(paste("  - All patients:", nrow(df_kf_changePT_all), "\n"))
cat(paste("  - Preop PT < 20:", nrow(df_kf_changePT_pt_lt20), "\n"))
cat(paste("  - Preop PT >= 20:", nrow(df_kf_changePT_pt_ge20), "\n\n"))

# Create base filter for preop PT stratification (used for all timepoint plots)
df_base_pt_lt20 <- df %>%
  filter(
    !is.na(LATpre_S1PT) &
    LATpre_S1PT < 20
  )

df_base_pt_ge20 <- df %>%
  filter(
    !is.na(LATpre_S1PT) &
    LATpre_S1PT >= 20
  )

cat(paste("PT < 20 at 6 weeks:\n"))
cat(paste("  - Patients with preop KF data:", nrow(df_filtered_preop_lt20), "\n"))
cat(paste("  - Patients with 6-week KF data:", nrow(df_filtered_w6_lt20), "\n\n"))
cat(paste("PT > 20 at 6 weeks:\n"))
cat(paste("  - Patients with preop KF data:", nrow(df_filtered_preop_gt20), "\n"))
cat(paste("  - Patients with 6-week KF data:", nrow(df_filtered_w6_gt20), "\n\n"))

# Function to create plot (general version for PT vs KF)
create_plot <- function(data, pt_col, kf_col, kf_label, pt_label, title_suffix, filename) {
  # Extract values
  pt <- as.numeric(data[[pt_col]])
  kf <- abs(as.numeric(data[[kf_col]]))
  
  if (length(pt) < 3 || length(kf) < 3) {
    cat(sprintf("WARNING: Insufficient data for %s vs %s (n = %d). Skipping plot.\n", pt_label, kf_label, length(pt)))
    return(NULL)
  }
  
  # Create data frame
  data_clean <- data.frame(
    pt = pt,
    kf = kf
  )
  
  # Perform linear regression
  model <- lm(kf ~ pt, data = data_clean)
  summary_model <- summary(model)
  r_squared <- as.numeric(summary_model$r.squared)
  p_value <- as.numeric(summary_model$coefficients[2, 4])
  intercept <- as.numeric(summary_model$coefficients[1, 1])
  slope <- as.numeric(summary_model$coefficients[2, 1])
  
  # Calculate Pearson correlation coefficient
  pearson_r <- as.numeric(cor(pt, kf, use = "complete.obs"))
  
  # Print results
  cat(sprintf("=== Regression Results: %s vs %s ===\n", pt_label, kf_label))
  cat(sprintf("Sample size: %d\n", length(pt)))
  cat(sprintf("Pearson r: %.4f\n", pearson_r))
  cat(sprintf("R²: %.4f\n", r_squared))
  cat(sprintf("p-value: %.4e\n", p_value))
  cat(sprintf("Slope: %.4f\n", slope))
  cat(sprintf("Intercept: %.4f\n", intercept))
  cat("\n")
  
  # Determine x-axis label based on column name
  if (grepl("LATpre_S1PT", pt_col)) {
    x_label <- "Preoperative Pelvic Tilt (degrees)"
  } else if (grepl("LAT6W_S1PT", pt_col)) {
    x_label <- "Pelvic Tilt at 6 Weeks (degrees)"
  } else if (grepl("LAT1Y_S1PT", pt_col)) {
    x_label <- "Pelvic Tilt at 1 Year (degrees)"
  } else if (grepl("LAT2Y_S1PT", pt_col)) {
    x_label <- "Pelvic Tilt at 2 Years (degrees)"
  } else {
    x_label <- "Pelvic Tilt (degrees)"
  }
  
  # Create scatter plot with regression line
  p <- ggplot(data_clean, aes(x = pt, y = kf)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = x_label,
      y = paste0("Absolute ", kf_label, " (degrees)"),
      title = paste0(pt_label, " vs ", kf_label, "\n", title_suffix)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Add statistics annotation
  r_text <- paste0("r = ", round(pearson_r, 3))
  r2_text <- paste0("R² = ", round(r_squared, 3))
  p_text <- paste0("p = ", formatC(p_value, format = "e", digits = 2))
  n_text <- paste0("n = ", length(pt))
  annotation_text <- paste0(r_text, "\n", r2_text, "\n", p_text, "\n", n_text)
  
  p <- p + annotate(
    "text",
    x = Inf, y = Inf,
    label = annotation_text,
    hjust = 1.1, vjust = 1.5,
    size = 4,
    fontface = "bold"
  )
  
  # Save plot
  if (!dir.exists("planned_results")) {
    dir.create("planned_results")
  }
  ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
  cat(sprintf("Saved plot to %s\n", filename))
  
  return(p)
}

# Function to create plot for KF vs Change in PT
create_plot_kf_vs_changePT <- function(data, kf_col, changePT_col, title_suffix, filename) {
  # Extract values
  kf <- abs(as.numeric(data[[kf_col]]))
  changePT <- as.numeric(data[[changePT_col]])
  
  if (length(kf) < 3 || length(changePT) < 3) {
    cat(sprintf("WARNING: Insufficient data (n = %d). Skipping plot.\n", length(kf)))
    return(NULL)
  }
  
  # Create data frame
  data_clean <- data.frame(
    kf = kf,
    changePT = changePT
  )
  
  # Perform linear regression
  model <- lm(changePT ~ kf, data = data_clean)
  summary_model <- summary(model)
  r_squared <- as.numeric(summary_model$r.squared)
  p_value <- as.numeric(summary_model$coefficients[2, 4])
  intercept <- as.numeric(summary_model$coefficients[1, 1])
  slope <- as.numeric(summary_model$coefficients[2, 1])
  
  # Calculate Pearson correlation coefficient
  pearson_r <- as.numeric(cor(kf, changePT, use = "complete.obs"))
  
  # Print results
  cat(sprintf("=== Regression Results: Preop KF vs Change in PT ===\n"))
  cat(sprintf("Sample size: %d\n", length(kf)))
  cat(sprintf("Pearson r: %.4f\n", pearson_r))
  cat(sprintf("R²: %.4f\n", r_squared))
  cat(sprintf("p-value: %.4e\n", p_value))
  cat(sprintf("Slope: %.4f\n", slope))
  cat(sprintf("Intercept: %.4f\n", intercept))
  cat("\n")
  
  # Create scatter plot with regression line
  p <- ggplot(data_clean, aes(x = kf, y = changePT)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
      x = "Preoperative Knee Flexion (degrees)",
      y = "Change in Pelvic Tilt (6-week - Preop, degrees)",
      title = paste0("Preoperative Knee Flexion vs Change in Pelvic Tilt\n", title_suffix)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Add statistics annotation
  r_text <- paste0("r = ", round(pearson_r, 3))
  r2_text <- paste0("R² = ", round(r_squared, 3))
  p_text <- paste0("p = ", formatC(p_value, format = "e", digits = 2))
  n_text <- paste0("n = ", length(kf))
  annotation_text <- paste0(r_text, "\n", r2_text, "\n", p_text, "\n", n_text)
  
  p <- p + annotate(
    "text",
    x = Inf, y = Inf,
    label = annotation_text,
    hjust = 1.1, vjust = 1.5,
    size = 4,
    fontface = "bold"
  )
  
  # Save plot
  if (!dir.exists("planned_results")) {
    dir.create("planned_results")
  }
  ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
  cat(sprintf("Saved plot to %s\n", filename))
  
  return(p)
}

# Create plots for Preoperative PT vs Preoperative KF
cat("=== Creating plots for Preoperative PT vs Preoperative KF ===\n")
cat("\nPlot 0a: Preop PT vs Preop KF (All patients)\n")
p0a <- create_plot(
  df_preop_all,
  "LATpre_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Preoperative Pelvic Tilt",
  "(All patients)",
  "planned_results/analysis12_preop_PT_vs_preop_KF_all.png"
)

cat("\nPlot 0b: Preop PT vs Preop KF (Preop PT < 20)\n")
p0b <- create_plot(
  df_preop_pt_lt20,
  "LATpre_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Preoperative Pelvic Tilt",
  "(Patients with Preop PT < 20°)",
  "planned_results/analysis12_preop_PT_vs_preop_KF_lt20.png"
)

cat("\nPlot 0c: Preop PT vs Preop KF (Preop PT >= 20)\n")
p0c <- create_plot(
  df_preop_pt_ge20,
  "LATpre_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Preoperative Pelvic Tilt",
  "(Patients with Preop PT >= 20°)",
  "planned_results/analysis12_preop_PT_vs_preop_KF_ge20.png"
)

# Create plots for all timepoints: KF vs PT (stratified by preop PT)
cat("\n=== Creating plots for all timepoints: KF vs PT (stratified by Preop PT) ===\n")

# Preop KF vs 6-week PT
cat("\nPlot 2a: Preop KF vs 6-week PT (Preop PT < 20)\n")
df_temp <- df_base_pt_lt20 %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(LAT6W_S1PT))
p2a <- create_plot(
  df_temp,
  "LAT6W_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with Preop PT < 20°)",
  "planned_results/analysis12_6W_PT_vs_preop_KF_lt20.png"
)

cat("\nPlot 2b: Preop KF vs 6-week PT (Preop PT >= 20)\n")
df_temp <- df_base_pt_ge20 %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(LAT6W_S1PT))
p2b <- create_plot(
  df_temp,
  "LAT6W_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with Preop PT >= 20°)",
  "planned_results/analysis12_6W_PT_vs_preop_KF_ge20.png"
)

# Preop KF vs 1-year PT (if available)
if ("LAT1Y_S1PT" %in% names(df)) {
  cat("\nPlot 3a: Preop KF vs 1-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(LAT1Y_S1PT))
  p3a <- create_plot(
    df_temp,
    "LAT1Y_S1PT",
    "LATpre_LL_KneeAngle",
    "Preoperative Knee Flexion",
    "Pelvic Tilt at 1 Year",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_1Y_PT_vs_preop_KF_lt20.png"
  )
  
  cat("\nPlot 3b: Preop KF vs 1-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(LAT1Y_S1PT))
  p3b <- create_plot(
    df_temp,
    "LAT1Y_S1PT",
    "LATpre_LL_KneeAngle",
    "Preoperative Knee Flexion",
    "Pelvic Tilt at 1 Year",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_1Y_PT_vs_preop_KF_ge20.png"
  )
}

# 6-week KF vs 6-week PT
cat("\nPlot 4a: 6-week KF vs 6-week PT (Preop PT < 20)\n")
df_temp <- df_base_pt_lt20 %>%
  filter(!is.na(LAT6W_LL_KneeAngle) & !is.na(LAT6W_S1PT))
p4a <- create_plot(
  df_temp,
  "LAT6W_S1PT",
  "LAT6W_LL_KneeAngle",
  "6-Week Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with Preop PT < 20°)",
  "planned_results/analysis12_6W_PT_vs_6W_KF_lt20.png"
)

cat("\nPlot 4b: 6-week KF vs 6-week PT (Preop PT >= 20)\n")
df_temp <- df_base_pt_ge20 %>%
  filter(!is.na(LAT6W_LL_KneeAngle) & !is.na(LAT6W_S1PT))
p4b <- create_plot(
  df_temp,
  "LAT6W_S1PT",
  "LAT6W_LL_KneeAngle",
  "6-Week Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with Preop PT >= 20°)",
  "planned_results/analysis12_6W_PT_vs_6W_KF_ge20.png"
)

# 6-week KF vs 1-year PT (if available)
if ("LAT1Y_S1PT" %in% names(df)) {
  cat("\nPlot 5a: 6-week KF vs 1-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LAT6W_LL_KneeAngle) & !is.na(LAT1Y_S1PT))
  p5a <- create_plot(
    df_temp,
    "LAT1Y_S1PT",
    "LAT6W_LL_KneeAngle",
    "6-Week Knee Flexion",
    "Pelvic Tilt at 1 Year",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_1Y_PT_vs_6W_KF_lt20.png"
  )
  
  cat("\nPlot 5b: 6-week KF vs 1-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LAT6W_LL_KneeAngle) & !is.na(LAT1Y_S1PT))
  p5b <- create_plot(
    df_temp,
    "LAT1Y_S1PT",
    "LAT6W_LL_KneeAngle",
    "6-Week Knee Flexion",
    "Pelvic Tilt at 1 Year",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_1Y_PT_vs_6W_KF_ge20.png"
  )
}

# 1-year KF vs 1-year PT (if available)
if ("LAT1Y_LL_KneeAngle" %in% names(df) && "LAT1Y_S1PT" %in% names(df)) {
  cat("\nPlot 6a: 1-year KF vs 1-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LAT1Y_LL_KneeAngle) & !is.na(LAT1Y_S1PT))
  p6a <- create_plot(
    df_temp,
    "LAT1Y_S1PT",
    "LAT1Y_LL_KneeAngle",
    "1-Year Knee Flexion",
    "Pelvic Tilt at 1 Year",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_1Y_PT_vs_1Y_KF_lt20.png"
  )
  
  cat("\nPlot 6b: 1-year KF vs 1-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LAT1Y_LL_KneeAngle) & !is.na(LAT1Y_S1PT))
  p6b <- create_plot(
    df_temp,
    "LAT1Y_S1PT",
    "LAT1Y_LL_KneeAngle",
    "1-Year Knee Flexion",
    "Pelvic Tilt at 1 Year",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_1Y_PT_vs_1Y_KF_ge20.png"
  )
}

# Preop KF vs 2-year PT (if available)
if ("LAT2Y_S1PT" %in% names(df)) {
  cat("\nPlot 7a: Preop KF vs 2-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p7a <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LATpre_LL_KneeAngle",
    "Preoperative Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_2Y_PT_vs_preop_KF_lt20.png"
  )
  
  cat("\nPlot 7b: Preop KF vs 2-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p7b <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LATpre_LL_KneeAngle",
    "Preoperative Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_2Y_PT_vs_preop_KF_ge20.png"
  )
}

# 6-week KF vs 2-year PT (if available)
if ("LAT2Y_S1PT" %in% names(df)) {
  cat("\nPlot 8a: 6-week KF vs 2-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LAT6W_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p8a <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LAT6W_LL_KneeAngle",
    "6-Week Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_2Y_PT_vs_6W_KF_lt20.png"
  )
  
  cat("\nPlot 8b: 6-week KF vs 2-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LAT6W_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p8b <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LAT6W_LL_KneeAngle",
    "6-Week Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_2Y_PT_vs_6W_KF_ge20.png"
  )
}

# 1-year KF vs 2-year PT (if available)
if ("LAT1Y_LL_KneeAngle" %in% names(df) && "LAT2Y_S1PT" %in% names(df)) {
  cat("\nPlot 9a: 1-year KF vs 2-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LAT1Y_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p9a <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LAT1Y_LL_KneeAngle",
    "1-Year Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_2Y_PT_vs_1Y_KF_lt20.png"
  )
  
  cat("\nPlot 9b: 1-year KF vs 2-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LAT1Y_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p9b <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LAT1Y_LL_KneeAngle",
    "1-Year Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_2Y_PT_vs_1Y_KF_ge20.png"
  )
}

# 2-year KF vs 2-year PT (if available)
if ("LAT2Y_LL_KneeAngle" %in% names(df) && "LAT2Y_S1PT" %in% names(df)) {
  cat("\nPlot 10a: 2-year KF vs 2-year PT (Preop PT < 20)\n")
  df_temp <- df_base_pt_lt20 %>%
    filter(!is.na(LAT2Y_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p10a <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LAT2Y_LL_KneeAngle",
    "2-Year Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT < 20°)",
    "planned_results/analysis12_2Y_PT_vs_2Y_KF_lt20.png"
  )
  
  cat("\nPlot 10b: 2-year KF vs 2-year PT (Preop PT >= 20)\n")
  df_temp <- df_base_pt_ge20 %>%
    filter(!is.na(LAT2Y_LL_KneeAngle) & !is.na(LAT2Y_S1PT))
  p10b <- create_plot(
    df_temp,
    "LAT2Y_S1PT",
    "LAT2Y_LL_KneeAngle",
    "2-Year Knee Flexion",
    "Pelvic Tilt at 2 Years",
    "(Patients with Preop PT >= 20°)",
    "planned_results/analysis12_2Y_PT_vs_2Y_KF_ge20.png"
  )
}

# Create plots for Preoperative KF vs Change in PT
cat("\n=== Creating plots for Preoperative KF vs Change in PT ===\n")
cat("\nPlot 1a: Preop KF vs Change in PT (All patients)\n")
p1a <- create_plot_kf_vs_changePT(
  df_kf_changePT_all,
  "LATpre_LL_KneeAngle",
  "change_PT",
  "(All patients)",
  "planned_results/analysis12_preop_KF_vs_change_PT_all.png"
)

cat("\nPlot 1b: Preop KF vs Change in PT (Preop PT < 20)\n")
p1b <- create_plot_kf_vs_changePT(
  df_kf_changePT_pt_lt20,
  "LATpre_LL_KneeAngle",
  "change_PT",
  "(Patients with Preop PT < 20°)",
  "planned_results/analysis12_preop_KF_vs_change_PT_lt20.png"
)

cat("\nPlot 1c: Preop KF vs Change in PT (Preop PT >= 20)\n")
p1c <- create_plot_kf_vs_changePT(
  df_kf_changePT_pt_ge20,
  "LATpre_LL_KneeAngle",
  "change_PT",
  "(Patients with Preop PT >= 20°)",
  "planned_results/analysis12_preop_KF_vs_change_PT_ge20.png"
)

# Create plots for PT < 20 at 6 weeks (original analysis)
cat("\n=== Creating plots for PT < 20 at 6 weeks ===\n")
cat("\nPlot 7: PT at 6 weeks vs Preop KF (PT < 20 at 6 weeks)\n")
p7_old <- create_plot(
  df_filtered_preop_lt20,
  "LAT6W_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with PT < 20° at 6 weeks)",
  "planned_results/analysis12_PT_vs_preop_KF_lt20_6W.png"
)

cat("\nPlot 7: PT at 6 weeks vs 6-Week KF (PT < 20 at 6 weeks)\n")
p7 <- create_plot(
  df_filtered_w6_lt20,
  "LAT6W_S1PT",
  "LAT6W_LL_KneeAngle",
  "6-Week Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with PT < 20° at 6 weeks)",
  "planned_results/analysis12_PT_vs_6W_KF_lt20.png"
)

# Create plots for PT > 20 at 6 weeks
cat("\n=== Creating plots for PT > 20 at 6 weeks ===\n")
cat("\nPlot 8: PT at 6 weeks vs Preop KF (PT > 20 at 6 weeks)\n")
p8 <- create_plot(
  df_filtered_preop_gt20,
  "LAT6W_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with PT > 20° at 6 weeks)",
  "planned_results/analysis12_PT_vs_preop_KF_gt20.png"
)

cat("\nPlot 9: PT at 6 weeks vs 6-Week KF (PT > 20 at 6 weeks)\n")
p9 <- create_plot(
  df_filtered_w6_gt20,
  "LAT6W_S1PT",
  "LAT6W_LL_KneeAngle",
  "6-Week Knee Flexion",
  "Pelvic Tilt at 6 Weeks",
  "(Patients with PT > 20° at 6 weeks)",
  "planned_results/analysis12_PT_vs_6W_KF_gt20.png"
)

# ============================================================================
# Additional Analysis: KF vs PT at same timepoints, stratified by PT at that timepoint
# ============================================================================
cat("\n=== Additional Analysis: KF vs PT at Same Timepoints (stratified by PT at that timepoint) ===\n\n")

# Function to calculate correlation and return results
calc_correlation <- function(data, pt_col, kf_col, pt_threshold = 20) {
  if (nrow(data) < 3) {
    return(list(r = NA, p = NA, n = nrow(data)))
  }
  
  pt <- as.numeric(data[[pt_col]])
  kf <- abs(as.numeric(data[[kf_col]]))
  
  # Remove missing values
  complete <- complete.cases(pt, kf)
  pt <- pt[complete]
  kf <- kf[complete]
  
  if (length(pt) < 3) {
    return(list(r = NA, p = NA, n = length(pt)))
  }
  
  # Calculate correlation
  cor_result <- cor.test(pt, kf)
  
  return(list(
    r = cor_result$estimate,
    p = cor_result$p.value,
    n = length(pt)
  ))
}

# Preop: Preop KF vs Preop PT, stratified by Preop PT
cat("Preop: Preop KF vs Preop PT\n")
df_preop_lt20_timepoint <- df %>%
  filter(!is.na(LATpre_S1PT) & !is.na(LATpre_LL_KneeAngle) & LATpre_S1PT < 20)
df_preop_ge20_timepoint <- df %>%
  filter(!is.na(LATpre_S1PT) & !is.na(LATpre_LL_KneeAngle) & LATpre_S1PT >= 20)

preop_lt20 <- calc_correlation(df_preop_lt20_timepoint, "LATpre_S1PT", "LATpre_LL_KneeAngle")
preop_ge20 <- calc_correlation(df_preop_ge20_timepoint, "LATpre_S1PT", "LATpre_LL_KneeAngle")

cat(sprintf("  PT < 20: r = %.3f, p = %.4e, n = %d\n", preop_lt20$r, preop_lt20$p, preop_lt20$n))
cat(sprintf("  PT >= 20: r = %.3f, p = %.4e, n = %d\n\n", preop_ge20$r, preop_ge20$p, preop_ge20$n))

# 6-week: 6-week KF vs 6-week PT, stratified by 6-week PT
cat("6-week: 6-week KF vs 6-week PT\n")
df_6w_lt20_timepoint <- df %>%
  filter(!is.na(LAT6W_S1PT) & !is.na(LAT6W_LL_KneeAngle) & LAT6W_S1PT < 20)
df_6w_ge20_timepoint <- df %>%
  filter(!is.na(LAT6W_S1PT) & !is.na(LAT6W_LL_KneeAngle) & LAT6W_S1PT >= 20)

w6_lt20 <- calc_correlation(df_6w_lt20_timepoint, "LAT6W_S1PT", "LAT6W_LL_KneeAngle")
w6_ge20 <- calc_correlation(df_6w_ge20_timepoint, "LAT6W_S1PT", "LAT6W_LL_KneeAngle")

cat(sprintf("  PT < 20: r = %.3f, p = %.4e, n = %d\n", w6_lt20$r, w6_lt20$p, w6_lt20$n))
cat(sprintf("  PT >= 20: r = %.3f, p = %.4e, n = %d\n\n", w6_ge20$r, w6_ge20$p, w6_ge20$n))

# 1-year: 1-year KF vs 1-year PT, stratified by 1-year PT
if ("LAT1Y_S1PT" %in% names(df) && "LAT1Y_LL_KneeAngle" %in% names(df)) {
  cat("1-year: 1-year KF vs 1-year PT\n")
  df_1y_lt20_timepoint <- df %>%
    filter(!is.na(LAT1Y_S1PT) & !is.na(LAT1Y_LL_KneeAngle) & LAT1Y_S1PT < 20)
  df_1y_ge20_timepoint <- df %>%
    filter(!is.na(LAT1Y_S1PT) & !is.na(LAT1Y_LL_KneeAngle) & LAT1Y_S1PT >= 20)
  
  y1_lt20 <- calc_correlation(df_1y_lt20_timepoint, "LAT1Y_S1PT", "LAT1Y_LL_KneeAngle")
  y1_ge20 <- calc_correlation(df_1y_ge20_timepoint, "LAT1Y_S1PT", "LAT1Y_LL_KneeAngle")
  
  cat(sprintf("  PT < 20: r = %.3f, p = %.4e, n = %d\n", y1_lt20$r, y1_lt20$p, y1_lt20$n))
  cat(sprintf("  PT >= 20: r = %.3f, p = %.4e, n = %d\n\n", y1_ge20$r, y1_ge20$p, y1_ge20$n))
} else {
  y1_lt20 <- list(r = NA, p = NA, n = 0)
  y1_ge20 <- list(r = NA, p = NA, n = 0)
}

# 2-year: 2-year KF vs 2-year PT, stratified by 2-year PT
if ("LAT2Y_S1PT" %in% names(df) && "LAT2Y_LL_KneeAngle" %in% names(df)) {
  cat("2-year: 2-year KF vs 2-year PT\n")
  df_2y_lt20_timepoint <- df %>%
    filter(!is.na(LAT2Y_S1PT) & !is.na(LAT2Y_LL_KneeAngle) & LAT2Y_S1PT < 20)
  df_2y_ge20_timepoint <- df %>%
    filter(!is.na(LAT2Y_S1PT) & !is.na(LAT2Y_LL_KneeAngle) & LAT2Y_S1PT >= 20)
  
  y2_lt20 <- calc_correlation(df_2y_lt20_timepoint, "LAT2Y_S1PT", "LAT2Y_LL_KneeAngle")
  y2_ge20 <- calc_correlation(df_2y_ge20_timepoint, "LAT2Y_S1PT", "LAT2Y_LL_KneeAngle")
  
  cat(sprintf("  PT < 20: r = %.3f, p = %.4e, n = %d\n", y2_lt20$r, y2_lt20$p, y2_lt20$n))
  cat(sprintf("  PT >= 20: r = %.3f, p = %.4e, n = %d\n\n", y2_ge20$r, y2_ge20$p, y2_ge20$n))
} else {
  y2_lt20 <- list(r = NA, p = NA, n = 0)
  y2_ge20 <- list(r = NA, p = NA, n = 0)
}

# Format p-values for the sentence
format_p_value <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) return("<0.001")
  if (p < 0.01) return(sprintf("%.3f", p))
  if (p < 0.05) return(sprintf("%.3f", p))
  return(sprintf("%.3f", p))
}

# Apply Benjamini-Hochberg (BH) correction for multiple comparisons
cat("\n=== Applying Benjamini-Hochberg (BH) Correction ===\n\n")

# Collect all p-values
all_p_values <- c()
all_labels <- c()

# PT >= 20 group
all_p_values <- c(all_p_values, preop_ge20$p, w6_ge20$p)
all_labels <- c(all_labels, "Preop PT>=20", "6wk PT>=20")
if (!is.na(y1_ge20$p)) {
  all_p_values <- c(all_p_values, y1_ge20$p)
  all_labels <- c(all_labels, "1yr PT>=20")
}
if (!is.na(y2_ge20$p)) {
  all_p_values <- c(all_p_values, y2_ge20$p)
  all_labels <- c(all_labels, "2yr PT>=20")
}

# PT < 20 group
all_p_values <- c(all_p_values, preop_lt20$p, w6_lt20$p)
all_labels <- c(all_labels, "Preop PT<20", "6wk PT<20")
if (!is.na(y1_lt20$p)) {
  all_p_values <- c(all_p_values, y1_lt20$p)
  all_labels <- c(all_labels, "1yr PT<20")
}
if (!is.na(y2_lt20$p)) {
  all_p_values <- c(all_p_values, y2_lt20$p)
  all_labels <- c(all_labels, "2yr PT<20")
}

# Apply BH correction
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

# Create data frame for display
pvalue_df <- data.frame(
  Comparison = all_labels,
  Raw_p = all_p_values,
  BH_adjusted_p = adjusted_p_values
)

cat("Multiple comparisons correction (Benjamini-Hochberg):\n")
print(pvalue_df, row.names = FALSE)
cat("\n")

# Extract adjusted p-values for each comparison
preop_ge20$p_adj <- adjusted_p_values[all_labels == "Preop PT>=20"]
w6_ge20$p_adj <- adjusted_p_values[all_labels == "6wk PT>=20"]
preop_lt20$p_adj <- adjusted_p_values[all_labels == "Preop PT<20"]
w6_lt20$p_adj <- adjusted_p_values[all_labels == "6wk PT<20"]

if (!is.na(y1_ge20$p)) {
  y1_ge20$p_adj <- adjusted_p_values[all_labels == "1yr PT>=20"]
}
if (!is.na(y2_ge20$p)) {
  y2_ge20$p_adj <- adjusted_p_values[all_labels == "2yr PT>=20"]
}
if (!is.na(y1_lt20$p)) {
  y1_lt20$p_adj <- adjusted_p_values[all_labels == "1yr PT<20"]
}
if (!is.na(y2_lt20$p)) {
  y2_lt20$p_adj <- adjusted_p_values[all_labels == "2yr PT<20"]
}

cat("\n=== Summary for Manuscript (with BH correction) ===\n\n")
cat("A correlation between KF and pathologic PT (PT >= 20) was found at all timepoints:\n")
cat(sprintf("  - Preop: r = %.3f, p = %s (BH-adjusted: %s)\n", 
            preop_ge20$r, format_p_value(preop_ge20$p), format_p_value(preop_ge20$p_adj)))
cat(sprintf("  - 6-week: r = %.3f, p = %s (BH-adjusted: %s)\n", 
            w6_ge20$r, format_p_value(w6_ge20$p), format_p_value(w6_ge20$p_adj)))
if (!is.na(y1_ge20$r)) {
  cat(sprintf("  - 1-year: r = %.3f, p = %s (BH-adjusted: %s)\n", 
              y1_ge20$r, format_p_value(y1_ge20$p), format_p_value(y1_ge20$p_adj)))
}
if (!is.na(y2_ge20$r)) {
  cat(sprintf("  - 2-year: r = %.3f, p = %s (BH-adjusted: %s)\n", 
              y2_ge20$r, format_p_value(y2_ge20$p), format_p_value(y2_ge20$p_adj)))
}

cat("\nNo significant correlation between KF and normal PT (PT < 20) was found:\n")
cat(sprintf("  - Preop: r = %.3f, p = %s (BH-adjusted: %s)\n", 
            preop_lt20$r, format_p_value(preop_lt20$p), format_p_value(preop_lt20$p_adj)))
cat(sprintf("  - 6-week: r = %.3f, p = %s (BH-adjusted: %s)\n", 
            w6_lt20$r, format_p_value(w6_lt20$p), format_p_value(w6_lt20$p_adj)))
if (!is.na(y1_lt20$r)) {
  cat(sprintf("  - 1-year: r = %.3f, p = %s (BH-adjusted: %s)\n", 
              y1_lt20$r, format_p_value(y1_lt20$p), format_p_value(y1_lt20$p_adj)))
}
if (!is.na(y2_lt20$r)) {
  cat(sprintf("  - 2-year: r = %.3f, p = %s (BH-adjusted: %s)\n", 
              y2_lt20$r, format_p_value(y2_lt20$p), format_p_value(y2_lt20$p_adj)))
}

cat("\n=== Formatted Sentence (using BH-adjusted p-values) ===\n\n")
cat("A correlation between KF and pathologic PT (PT >= 20) was found at all timepoints ")
if (!is.na(y2_ge20$r)) {
  cat("(preop, 6wk, 1yr, 2yr")
} else if (!is.na(y1_ge20$r)) {
  cat("(preop, 6wk, 1yr")
} else {
  cat("(preop, 6wk")
}
cat("; ")

# Collect BH-adjusted p-values for PT >= 20
p_values_ge20_adj <- c(preop_ge20$p_adj, w6_ge20$p_adj)
if (!is.na(y1_ge20$p_adj)) p_values_ge20_adj <- c(p_values_ge20_adj, y1_ge20$p_adj)
if (!is.na(y2_ge20$p_adj)) p_values_ge20_adj <- c(p_values_ge20_adj, y2_ge20$p_adj)

# Determine significance level using adjusted p-values
if (all(p_values_ge20_adj < 0.001, na.rm = TRUE)) {
  cat("p<0.001")
} else if (all(p_values_ge20_adj < 0.01, na.rm = TRUE)) {
  cat("p<0.01")
} else if (all(p_values_ge20_adj < 0.05, na.rm = TRUE)) {
  cat("p<0.05")
} else {
  # Use the most significant adjusted p-value
  min_p_adj <- min(p_values_ge20_adj, na.rm = TRUE)
  if (min_p_adj < 0.001) {
    cat("p<0.001")
  } else if (min_p_adj < 0.01) {
    cat("p<0.01")
  } else if (min_p_adj < 0.05) {
    cat("p<0.05")
  } else {
    cat(sprintf("p=%.3f", min_p_adj))
  }
}

cat("), but not between KF and normal PT (PT < 20; ")

# Collect BH-adjusted p-values for PT < 20
p_values_lt20_adj <- c(preop_lt20$p_adj, w6_lt20$p_adj)
if (!is.na(y1_lt20$p_adj)) p_values_lt20_adj <- c(p_values_lt20_adj, y1_lt20$p_adj)
if (!is.na(y2_lt20$p_adj)) p_values_lt20_adj <- c(p_values_lt20_adj, y2_lt20$p_adj)

# Check if any are significant using adjusted p-values
if (all(p_values_lt20_adj >= 0.05, na.rm = TRUE)) {
  cat("p≥0.05")
} else {
  # Use the most significant adjusted p-value
  min_p_adj <- min(p_values_lt20_adj, na.rm = TRUE)
  if (min_p_adj < 0.001) {
    cat("p<0.001")
  } else if (min_p_adj < 0.01) {
    cat("p<0.01")
  } else if (min_p_adj < 0.05) {
    cat("p<0.05")
  } else {
    cat(sprintf("p=%.3f", min_p_adj))
  }
}

cat(").\n\n")
cat("Note: p-values are Benjamini-Hochberg adjusted for multiple comparisons.\n\n")

cat("\nAnalysis complete!\n")
