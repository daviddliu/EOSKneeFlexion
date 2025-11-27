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

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Filter patients with satisfactory PI-LL mismatch (between -10 and 10)
# Drop patients with missing PI-LL data (using 6-week data)
df_filtered <- df %>%
  filter(LAT6W_PI_LL >= -10 & LAT6W_PI_LL <= 10, !is.na(LAT6W_PI_LL))

cat(paste("Filtered to", nrow(df_filtered), "patients with PI-LL between -10 and 10 (6-week)\n"))

# Calculate difference in knee angle (using 6-week data)
# Note: If either value is missing, knee_angle_diff will be NA (these will be dropped in each plot)
df_filtered$knee_angle_diff <- df_filtered$LAT6W_LL_KneeAngle - df_filtered$LATpre_LL_KneeAngle

# Function to create plot with regression line, R², and p-value
create_plot <- function(data, x_var, y_var, x_label, y_label, title, filename) {
  # Drop patients with missing data for this specific variable pair (listwise deletion)
  # No imputation is performed - patients with missing values are excluded
  data_clean <- data %>%
    filter(!is.na(!!sym(x_var)) & !is.na(!!sym(y_var)))
  
  n_dropped <- nrow(data) - nrow(data_clean)
  if (n_dropped > 0) {
    cat(paste("  Dropped", n_dropped, "patients with missing", x_var, "or", y_var, "data\n"))
  }
  
  if (nrow(data_clean) == 0) {
    cat(paste("Warning: No valid data for", title, "\n"))
    return(NULL)
  }
  
  # Perform linear regression
  formula <- as.formula(paste(y_var, "~", x_var))
  model <- lm(formula, data = data_clean)
  summary_model <- summary(model)
  r_squared <- summary_model$r.squared
  p_value <- summary_model$coefficients[2, 4]
  
  # Calculate Pearson correlation coefficient
  pearson_r <- cor(data_clean[[x_var]], data_clean[[y_var]], use = "complete.obs")
  
  # Create plot
  p <- ggplot(data_clean, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = x_label,
      y = y_label,
      title = title
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Add Pearson r, R² and p-value annotation
  r_text <- paste0("r = ", round(pearson_r, 3))
  r2_text <- paste0("R² = ", round(r_squared, 3))
  p_text <- paste0("p = ", formatC(p_value, format = "e", digits = 2))
  annotation_text <- paste0(r_text, "\n", r2_text, "\n", p_text)
  
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
  filepath <- file.path("results", filename)
  ggsave(filepath, plot = p, width = 10, height = 8, dpi = 300)
  cat(paste("Saved plot to", filepath, "\n"))
  
  return(p)
}

# Ensure results directory exists
if (!dir.exists("results")) {
  dir.create("results")
}

# Plot 1: PI-LL vs Preop Knee Angle
create_plot(
  df_filtered,
  "LAT6W_PI_LL",
  "LATpre_LL_KneeAngle",
  "PI-LL Mismatch (6-week)",
  "Preoperative Knee Flexion",
  "PI-LL Mismatch (6-week) vs Preoperative Knee Flexion\n(Patients with PI-LL between -10 and 10)",
  "analysis2_correctedLL_kf_1.png"
)

# Plot 2: PI-LL vs Postop Knee Angle
create_plot(
  df_filtered,
  "LAT6W_PI_LL",
  "LAT6W_LL_KneeAngle",
  "PI-LL Mismatch (6-week)",
  "Postoperative Knee Flexion (6-week)",
  "PI-LL Mismatch (6-week) vs Postoperative Knee Flexion (6-week)\n(Patients with PI-LL between -10 and 10)",
  "analysis2_correctedLL_kf_2.png"
)

# Plot 3: PI-LL vs Change in Knee Angle
create_plot(
  df_filtered,
  "LAT6W_PI_LL",
  "knee_angle_diff",
  "PI-LL Mismatch (6-week)",
  "Change in Knee Flexion (6-week)",
  "PI-LL Mismatch (6-week) vs Change in Knee Flexion (6-week)\n(Patients with PI-LL between -10 and 10)",
  "analysis2_correctedLL_kf_3.png"
)

cat("All plots completed\n")

