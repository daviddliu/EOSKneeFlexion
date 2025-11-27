#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(ordinal)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration: Toggle PJK exclusion on/off
EXCLUDE_PJK <- TRUE

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
cat("\n")

# Print statistics on data availability for variables of interest
cat("\n=== Data Availability Statistics ===\n")
cat(paste("Total number of patients:", nrow(df), "\n\n"))

# Variables of interest
variables_of_interest <- list(
  "LATpre_LL_KneeAngle" = "Preoperative Knee Flexion",
  "LAT1Y_LL_KneeAngle" = "Postoperative Knee Flexion (1 year)",
  "LATpre_L1_S1" = "Preoperative Lordosis (L1-S1)",
  "LAT1Y_L1_S1" = "Postoperative Lordosis (L1-S1, 1 year)",
  "LAT1Y_PI_LL" = "PI-LL Mismatch (1 year)",
  "SxI_Gen_Surgeon" = "Surgeon"
)

for (var in names(variables_of_interest)) {
  var_name <- variables_of_interest[[var]]
  n_available <- sum(!is.na(df[[var]]))
  n_total <- nrow(df)
  pct_available <- round(100 * n_available / n_total, 1)
  n_missing <- n_total - n_available
  
  cat(sprintf("%s (%s):\n", var_name, var))
  cat(sprintf("  Available: %d (%.1f%%)\n", n_available, pct_available))
  cat(sprintf("  Missing: %d (%.1f%%)\n", n_missing, 100 - pct_available))
  cat("\n")
}

# Calculate derived variables availability
cat("=== Derived Variables ===\n")

# Change in lordosis
df$change_lordosis_temp <- df$LAT1Y_L1_S1 - df$LATpre_L1_S1
n_change_lordosis <- sum(!is.na(df$change_lordosis_temp))
cat(sprintf("Change in Lordosis (LAT1Y_L1_S1 - LATpre_L1_S1):\n"))
cat(sprintf("  Available: %d (%.1f%%)\n", n_change_lordosis, round(100 * n_change_lordosis / nrow(df), 1)))
cat("\n")

# Change in knee flexion
df$change_knee_temp <- df$LAT1Y_LL_KneeAngle - df$LATpre_LL_KneeAngle
n_change_knee <- sum(!is.na(df$change_knee_temp))
cat(sprintf("Change in Knee Flexion (LAT1Y_LL_KneeAngle - LATpre_LL_KneeAngle):\n"))
cat(sprintf("  Available: %d (%.1f%%)\n", n_change_knee, round(100 * n_change_knee / nrow(df), 1)))
cat("\n")

# Analysis-specific combinations
cat("=== Analysis-Specific Data Availability ===\n")

# Analysis 1: Preop knee flexion vs change in lordosis
# (using change_lordosis_temp already calculated above)
n_analysis1 <- sum(!is.na(df$LATpre_LL_KneeAngle) & !is.na(df$change_lordosis_temp))
cat(sprintf("Analysis 1 (Preop Knee Flexion vs Change in Lordosis):\n"))
cat(sprintf("  Available: %d (%.1f%%)\n", n_analysis1, round(100 * n_analysis1 / nrow(df), 1)))
cat("\n")

# Analysis 2: PI-LL < 10 patients
n_pi_ll_available <- sum(!is.na(df$LAT1Y_PI_LL))
n_pi_ll_less10 <- sum(df$LAT1Y_PI_LL < 10, na.rm = TRUE)
cat(sprintf("Analysis 2 (PI-LL < 10):\n"))
cat(sprintf("  Patients with PI-LL data: %d (%.1f%%)\n", n_pi_ll_available, round(100 * n_pi_ll_available / nrow(df), 1)))
cat(sprintf("  Patients with PI-LL < 10: %d (%.1f%% of total, %.1f%% of those with PI-LL data)\n", 
            n_pi_ll_less10, 
            round(100 * n_pi_ll_less10 / nrow(df), 1),
            round(100 * n_pi_ll_less10 / n_pi_ll_available, 1)))
cat("\n")

# Clean up temporary variables
df$change_lordosis_temp <- NULL
df$change_knee_temp <- NULL

cat("=== End of Data Availability Statistics ===\n\n")
