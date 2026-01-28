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

cat(paste("\n=== Analysis 12: PT vs KF (stratified by PT at 6 weeks) ===\n"))
cat(paste("Total patients in database:", nrow(df), "\n\n"))
cat(paste("PT < 20 at 6 weeks:\n"))
cat(paste("  - Patients with preop KF data:", nrow(df_filtered_preop_lt20), "\n"))
cat(paste("  - Patients with 6-week KF data:", nrow(df_filtered_w6_lt20), "\n\n"))
cat(paste("PT > 20 at 6 weeks:\n"))
cat(paste("  - Patients with preop KF data:", nrow(df_filtered_preop_gt20), "\n"))
cat(paste("  - Patients with 6-week KF data:", nrow(df_filtered_w6_gt20), "\n\n"))

# Function to create plot
create_plot <- function(data, pt_col, kf_col, kf_label, title_suffix, filename) {
  # Extract values
  pt <- as.numeric(data[[pt_col]])
  kf <- abs(as.numeric(data[[kf_col]]))
  
  if (length(pt) < 3 || length(kf) < 3) {
    cat(sprintf("WARNING: Insufficient data for %s (n = %d). Skipping plot.\n", kf_label, length(pt)))
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
  cat(sprintf("=== Regression Results: %s ===\n", kf_label))
  cat(sprintf("Sample size: %d\n", length(pt)))
  cat(sprintf("Pearson r: %.4f\n", pearson_r))
  cat(sprintf("R²: %.4f\n", r_squared))
  cat(sprintf("p-value: %.4e\n", p_value))
  cat(sprintf("Slope: %.4f\n", slope))
  cat(sprintf("Intercept: %.4f\n", intercept))
  cat("\n")
  
  # Create scatter plot with regression line
  p <- ggplot(data_clean, aes(x = pt, y = kf)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = "Pelvic Tilt at 6 Weeks (degrees)",
      y = paste0("Absolute ", kf_label, " (degrees)"),
      title = paste0("Pelvic Tilt vs ", kf_label, "\n", title_suffix)
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
  if (!dir.exists("results")) {
    dir.create("results")
  }
  ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
  cat(sprintf("Saved plot to %s\n", filename))
  
  return(p)
}

# Create plots for PT < 20
cat("=== Creating plots for PT < 20 at 6 weeks ===\n")
cat("\nPlot 1: PT at 6 weeks vs Preop KF (PT < 20)\n")
p1 <- create_plot(
  df_filtered_preop_lt20,
  "LAT6W_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "(Patients with PT < 20° at 6 weeks)",
  "results/analysis12_PT_vs_preop_KF_lt20.png"
)

cat("\nPlot 2: PT at 6 weeks vs 6-Week KF (PT < 20)\n")
p2 <- create_plot(
  df_filtered_w6_lt20,
  "LAT6W_S1PT",
  "LAT6W_LL_KneeAngle",
  "6-Week Knee Flexion",
  "(Patients with PT < 20° at 6 weeks)",
  "results/analysis12_PT_vs_6W_KF_lt20.png"
)

# Create plots for PT > 20
cat("\n=== Creating plots for PT > 20 at 6 weeks ===\n")
cat("\nPlot 3: PT at 6 weeks vs Preop KF (PT > 20)\n")
p3 <- create_plot(
  df_filtered_preop_gt20,
  "LAT6W_S1PT",
  "LATpre_LL_KneeAngle",
  "Preoperative Knee Flexion",
  "(Patients with PT > 20° at 6 weeks)",
  "results/analysis12_PT_vs_preop_KF_gt20.png"
)

cat("\nPlot 4: PT at 6 weeks vs 6-Week KF (PT > 20)\n")
p4 <- create_plot(
  df_filtered_w6_gt20,
  "LAT6W_S1PT",
  "LAT6W_LL_KneeAngle",
  "6-Week Knee Flexion",
  "(Patients with PT > 20° at 6 weeks)",
  "results/analysis12_PT_vs_6W_KF_gt20.png"
)

cat("\nAnalysis complete!\n")
