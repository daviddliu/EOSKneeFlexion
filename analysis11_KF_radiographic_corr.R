#!/usr/bin/env Rscript
#
# Analysis 11: KF correlation with radiographic parameters at all timepoints
#
# Hypothesis: PT-KF decoupling is the biggest change from pre-op to post-op
#
# This script correlates KF with all radiographic parameters used in the PCA model
# at all timepoints (Pre, 6W, 1Y, 2Y) and creates a table where:
#   - Rows = timepoints
#   - Columns = radiographic parameters
#   - Values = correlation coefficients (r)
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(knitr)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration
EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

# Load data
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Check for age variable
age_var <- NULL
if ("demo_Age" %in% names(df)) {
  age_var <- "demo_Age"
} else if ("demo_age" %in% names(df)) {
  age_var <- "demo_age"
} else if ("age" %in% names(df)) {
  age_var <- "age"
} else if ("Age" %in% names(df)) {
  age_var <- "Age"
}

# Calculate T4-L1 PA for all timepoints
timepoint_prefixes <- c("LATpre", "LAT6W", "LAT1Y", "LAT2Y")
for (prefix in timepoint_prefixes) {
  l1pa_col <- paste0(prefix, "_L1PA")
  t4pa_col <- paste0(prefix, "_T4PA")
  t4l1pa_col <- paste0(prefix, "_T4_L1_PA")
  
  if (l1pa_col %in% names(df) && t4pa_col %in% names(df)) {
    df[[t4l1pa_col]] <- df[[l1pa_col]] - df[[t4pa_col]]
  } else {
    df[[t4l1pa_col]] <- NA
  }
}

cat("\n")
cat("================================================================================\n")
cat("Analysis 11: KF correlation with radiographic parameters at all timepoints\n")
cat("================================================================================\n\n")

# Define radiographic parameters (same as PCA model)
# Note: Age is constant across timepoints, so we'll handle it separately
radiographic_params <- list(
  list(name = "S1PI", 
       pre = "LATpre_S1PI", w6 = "LAT6W_S1PI", y1 = "LAT1Y_S1PI", y2 = "LAT2Y_S1PI",
       label = "S1PI"),
  list(name = "L1_S1", 
       pre = "LATpre_L1_S1", w6 = "LAT6W_L1_S1", y1 = "LAT1Y_L1_S1", y2 = "LAT2Y_L1_S1",
       label = "Lordosis L1-S1"),
  list(name = "L4_S1", 
       pre = "LATpre_L4_S1", w6 = "LAT6W_L4_S1", y1 = "LAT1Y_L4_S1", y2 = "LAT2Y_L4_S1",
       label = "Lordosis L4-S1"),
  list(name = "T2_T12", 
       pre = "LATpre_T2_T12", w6 = "LAT6W_T2_T12", y1 = "LAT1Y_T2_T12", y2 = "LAT2Y_T2_T12",
       label = "Thoracic Kyphosis T2-T12"),
  list(name = "S1PT", 
       pre = "LATpre_S1PT", w6 = "LAT6W_S1PT", y1 = "LAT1Y_S1PT", y2 = "LAT2Y_S1PT",
       label = "Pelvic Tilt S1PT"),
  list(name = "SVA", 
       pre = "LATpre_SVA_C2_S1", w6 = "LAT6W_SVA_C2_S1", y1 = "LAT1Y_SVA_C2_S1", y2 = "LAT2Y_SVA_C2_S1",
       label = "SVA"),
  list(name = "T4_L1_PA", 
       pre = "LATpre_T4_L1_PA", w6 = "LAT6W_T4_L1_PA", y1 = "LAT1Y_T4_L1_PA", y2 = "LAT2Y_T4_L1_PA",
       label = "T4-L1 PA")
)

# KF columns for each timepoint
kf_cols <- list(
  pre = "LATpre_LL_KneeAngle",
  w6 = "LAT6W_LL_KneeAngle",
  y1 = "LAT1Y_LL_KneeAngle",
  y2 = "LAT2Y_LL_KneeAngle"
)

# Initialize correlation matrix
timepoints <- c("Pre", "6W", "1Y", "2Y")
param_names <- sapply(radiographic_params, function(x) x$name)
if (!is.null(age_var)) {
  param_names <- c(param_names, "Age")
}

cor_matrix <- matrix(NA, nrow = length(timepoints), ncol = length(param_names))
rownames(cor_matrix) <- timepoints
colnames(cor_matrix) <- param_names

# Also store sample sizes
n_matrix <- matrix(NA, nrow = length(timepoints), ncol = length(param_names))
rownames(n_matrix) <- timepoints
colnames(n_matrix) <- param_names

# Calculate correlations for each timepoint and parameter
cat("Calculating correlations...\n\n")

for (tp_idx in 1:length(timepoints)) {
  tp_name <- timepoints[tp_idx]
  tp_short <- c("pre", "w6", "y1", "y2")[tp_idx]
  kf_col <- kf_cols[[tp_short]]
  
  if (!kf_col %in% names(df)) {
    cat(sprintf("  %s: KF column not found, skipping\n", tp_name))
    next
  }
  
  cat(sprintf("  %s:\n", tp_name))
  
  # Correlate KF with each radiographic parameter
  for (param_idx in 1:length(radiographic_params)) {
    param <- radiographic_params[[param_idx]]
    param_col <- switch(tp_short,
                       "pre" = param$pre,
                       "w6" = param$w6,
                       "y1" = param$y1,
                       "y2" = param$y2)
    
    if (param_col %in% names(df)) {
      # Get valid pairs
      kf_vec <- as.numeric(df[[kf_col]])
      param_vec <- as.numeric(df[[param_col]])
      valid <- !is.na(kf_vec) & !is.na(param_vec)
      
      if (sum(valid) >= 3) {
        r <- cor(kf_vec[valid], param_vec[valid], use = "complete.obs")
        cor_matrix[tp_idx, param_idx] <- r
        n_matrix[tp_idx, param_idx] <- sum(valid)
        cat(sprintf("    %s: r = %.4f (n = %d)\n", param$label, r, sum(valid)))
      } else {
        cat(sprintf("    %s: insufficient data (n = %d)\n", param$label, sum(valid)))
      }
    } else {
      cat(sprintf("    %s: column not found\n", param$label))
    }
  }
  
  # Handle Age separately (constant across timepoints)
  if (!is.null(age_var) && age_var %in% names(df)) {
    kf_vec <- as.numeric(df[[kf_col]])
    age_vec <- as.numeric(df[[age_var]])
    valid <- !is.na(kf_vec) & !is.na(age_vec)
    
    if (sum(valid) >= 3) {
      r <- cor(kf_vec[valid], age_vec[valid], use = "complete.obs")
      cor_matrix[tp_idx, length(param_names)] <- r
      n_matrix[tp_idx, length(param_names)] <- sum(valid)
      cat(sprintf("    Age: r = %.4f (n = %d)\n", r, sum(valid)))
    }
  }
  
  cat("\n")
}

# Create formatted table
cat("================================================================================\n")
cat("Correlation Table: KF vs Radiographic Parameters\n")
cat("================================================================================\n")
cat("Rows = Timepoints, Columns = Radiographic Parameters\n")
cat("Values = Pearson correlation coefficient (r)\n\n")

# Print correlation matrix
cat("Correlation Coefficients:\n")
print(round(cor_matrix, 4))
cat("\n")

# Print sample sizes
cat("Sample Sizes:\n")
print(n_matrix)
cat("\n")

# Calculate changes from pre-op to post-op
cat("================================================================================\n")
cat("Change in Correlation from Pre-op to Post-op (Δr = Pre - Post)\n")
cat("================================================================================\n\n")

if (sum(!is.na(cor_matrix[1,])) > 0) {
  for (param_idx in 1:ncol(cor_matrix)) {
    param_name <- colnames(cor_matrix)[param_idx]
    r_pre <- cor_matrix[1, param_idx]
    
    if (!is.na(r_pre)) {
      cat(sprintf("%s:\n", param_name))
      cat(sprintf("  Pre-op: r = %.4f\n", r_pre))
      
      for (tp_idx in 2:length(timepoints)) {
        tp_name <- timepoints[tp_idx]
        r_post <- cor_matrix[tp_idx, param_idx]
        
        if (!is.na(r_post)) {
          delta_r <- r_pre - r_post
          pct_change <- (delta_r / abs(r_pre)) * 100
          cat(sprintf("  %s: r = %.4f, Δr = %.4f (%.1f%% %s)\n", 
                      tp_name, r_post, delta_r, abs(pct_change),
                      ifelse(delta_r > 0, "decrease", "increase")))
        }
      }
      cat("\n")
    }
  }
}

# Identify which parameter shows the biggest decoupling (largest decrease in correlation)
cat("================================================================================\n")
cat("Biggest Decoupling (Largest Decrease in Correlation from Pre-op)\n")
cat("================================================================================\n\n")

max_decreases <- data.frame(
  Parameter = character(),
  Timepoint = character(),
  Pre_r = numeric(),
  Post_r = numeric(),
  Delta_r = numeric(),
  Pct_decrease = numeric(),
  stringsAsFactors = FALSE
)

for (param_idx in 1:ncol(cor_matrix)) {
  param_name <- colnames(cor_matrix)[param_idx]
  r_pre <- cor_matrix[1, param_idx]
  
  if (!is.na(r_pre)) {
    for (tp_idx in 2:length(timepoints)) {
      tp_name <- timepoints[tp_idx]
      r_post <- cor_matrix[tp_idx, param_idx]
      
      if (!is.na(r_post)) {
        delta_r <- r_pre - r_post
        if (delta_r > 0) {  # Only consider decreases
          pct_decrease <- (delta_r / abs(r_pre)) * 100
          max_decreases <- rbind(max_decreases, data.frame(
            Parameter = param_name,
            Timepoint = tp_name,
            Pre_r = r_pre,
            Post_r = r_post,
            Delta_r = delta_r,
            Pct_decrease = pct_decrease,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

if (nrow(max_decreases) > 0) {
  # Sort by absolute decrease
  max_decreases <- max_decreases[order(-max_decreases$Delta_r), ]
  
  cat("Top 5 largest decreases in correlation:\n")
  print(head(max_decreases, 5))
  cat("\n")
  
  # Overall biggest
  biggest <- max_decreases[1, ]
  cat(sprintf("BIGGEST DECOUPLING: %s at %s\n", biggest$Parameter, biggest$Timepoint))
  cat(sprintf("  Pre-op: r = %.4f\n", biggest$Pre_r))
  cat(sprintf("  %s: r = %.4f\n", biggest$Timepoint, biggest$Post_r))
  cat(sprintf("  Decrease: Δr = %.4f (%.1f%% decrease)\n", biggest$Delta_r, biggest$Pct_decrease))
  cat("\n")
  
  # Test hypothesis: Is PT-KF the biggest?
  pt_kf_decreases <- max_decreases[max_decreases$Parameter == "S1PT", ]
  if (nrow(pt_kf_decreases) > 0) {
    max_pt_kf <- pt_kf_decreases[which.max(pt_kf_decreases$Delta_r), ]
    cat("PT-KF decoupling:\n")
    cat(sprintf("  Largest decrease: %s, Δr = %.4f (%.1f%%)\n", 
                max_pt_kf$Timepoint, max_pt_kf$Delta_r, max_pt_kf$Pct_decrease))
    
    if (biggest$Parameter == "S1PT") {
      cat("  → HYPOTHESIS CONFIRMED: PT-KF is the biggest decoupling!\n")
    } else {
      cat(sprintf("  → HYPOTHESIS NOT CONFIRMED: %s shows bigger decoupling (Δr = %.4f)\n",
                  biggest$Parameter, biggest$Delta_r))
    }
  }
} else {
  cat("No decreases in correlation found.\n")
}

# Save correlation matrix to CSV
if (!dir.exists("planned_results")) {
  dir.create("planned_results")
}

# Save as CSV
write.csv(cor_matrix, "planned_results/analysis11_KF_radiographic_correlations.csv", row.names = TRUE)
cat("\nSaved correlation matrix to: planned_results/analysis11_KF_radiographic_correlations.csv\n")

# Save sample sizes
write.csv(n_matrix, "planned_results/analysis11_KF_radiographic_correlations_n.csv", row.names = TRUE)
cat("Saved sample sizes to: planned_results/analysis11_KF_radiographic_correlations_n.csv\n")

# Create a formatted table for display (using knitr::kable if available)
cat("\n================================================================================\n")
cat("Formatted Correlation Table\n")
cat("================================================================================\n\n")

# Convert to data frame for better display
cor_df <- as.data.frame(cor_matrix)
cor_df$Timepoint <- rownames(cor_df)
cor_df <- cor_df[, c("Timepoint", setdiff(colnames(cor_df), "Timepoint"))]

print(cor_df)
cat("\n")
