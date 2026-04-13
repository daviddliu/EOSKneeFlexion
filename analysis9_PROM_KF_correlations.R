#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(readxl)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration: Toggle PJK exclusion on/off
EXCLUDE_PJK <- TRUE

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Load HRQL sheet
cat("Loading HRQL sheet...\n")
hrql_df <- readxl::read_excel(
  db_path,
  sheet = "HRQL",
  .name_repair = ~make.unique(.x, sep = "__")
)
cat("Loaded HRQL sheet\n")

# Find demo_id column in main dataframe
demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || length(demo_id_col) == 0) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
  if (is.na(demo_id_col) || length(demo_id_col) == 0) {
    demo_id_col <- "demo_id"
  }
}

# Find BL_demo_id column in HRQL dataframe (could be renamed by .name_repair)
hrql_demo_id_col <- grep("^BL_demo_id$", names(hrql_df), value = TRUE)[1]
if (is.na(hrql_demo_id_col) || length(hrql_demo_id_col) == 0) {
  hrql_demo_id_col <- grep("^BL_demo_id", names(hrql_df), value = TRUE)[1]
  if (is.na(hrql_demo_id_col) || length(hrql_demo_id_col) == 0) {
    hrql_demo_id_col <- "BL_demo_id"
  }
}

cat(sprintf("Using %s from main dataframe and %s from HRQL sheet for join\n", demo_id_col, hrql_demo_id_col))

# Merge HRQL data with main dataframe
df <- df %>%
  left_join(hrql_df, by = setNames(hrql_demo_id_col, demo_id_col))
cat("Merged HRQL data with main dataframe\n\n")

# Debug: Show all HRQL columns (in case they were renamed by .name_repair)
hrql_cols <- grep("HRQL", names(df), value = TRUE)
cat("All columns containing 'HRQL':\n")
for (col in hrql_cols) {
  n_available <- sum(!is.na(df[[col]]))
  cat(sprintf("  %s: %d non-NA values\n", col, n_available))
}
cat("\n")

# Find PROM columns
# Pattern: *_HRQL_SRS_TOTAL and *_HRQL_ODI where * is timepoint
# Note: columns may have been renamed by .name_repair (e.g., BL_HRQL_ODI__1)
prom_cols <- grep("_HRQL_(SRS_TOTAL|ODI)", names(df), value = TRUE)
cat("Found PROM columns:\n")
for (col in prom_cols) {
  n_available <- sum(!is.na(df[[col]]))
  cat(sprintf("  %s: %d non-NA values\n", col, n_available))
}
cat("\n")

# Identify specific PROM columns for each timepoint
# Pattern: *_HRQL_SRS_TOTAL and *_HRQL_ODI where * is timepoint
# Use case-insensitive matching
# Note: columns may have been renamed by .name_repair (e.g., BL_HRQL_ODI__1)
# ODI: grep("^BL_HRQL_ODI") wrongly picks BL_HRQL_ODI_Date first — use total score only.
find_odi_total_col <- function(prefix, dat) {
  hit <- grep(paste0("^", prefix, "_HRQL_ODI"), names(dat), value = TRUE, ignore.case = TRUE)
  hit <- hit[grepl(paste0("^", prefix, "_HRQL_ODI(__[0-9]+)?$"), hit, ignore.case = TRUE)]
  exact <- paste0(prefix, "_HRQL_ODI")
  hit <- hit[order(hit != exact)]
  if (length(hit)) hit[1] else NA_character_
}

# Preoperative (baseline) - uses BL_ prefix
pre_srs_col <- grep("^BL_HRQL_SRS_TOTAL", names(df), value = TRUE, ignore.case = TRUE)[1]
pre_odi_col <- find_odi_total_col("BL", df)

# 6 weeks - uses W6_ prefix
w6_srs_col <- grep("^W6_HRQL_SRS_TOTAL", names(df), value = TRUE, ignore.case = TRUE)[1]
w6_odi_col <- find_odi_total_col("W6", df)

# 2 years - uses Y2_ prefix
y2_srs_col <- grep("^Y2_HRQL_SRS_TOTAL", names(df), value = TRUE, ignore.case = TRUE)[1]
y2_odi_col <- find_odi_total_col("Y2", df)

# Check which columns were found
cat("=== PROM Column Identification ===\n")
cat(sprintf("Preoperative SRS: %s\n", ifelse(is.na(pre_srs_col), "NOT FOUND", pre_srs_col)))
cat(sprintf("Preoperative ODI: %s\n", ifelse(is.na(pre_odi_col), "NOT FOUND", pre_odi_col)))
cat(sprintf("6-week SRS: %s\n", ifelse(is.na(w6_srs_col), "NOT FOUND", w6_srs_col)))
cat(sprintf("6-week ODI: %s\n", ifelse(is.na(w6_odi_col), "NOT FOUND", w6_odi_col)))
cat(sprintf("2-year SRS: %s\n", ifelse(is.na(y2_srs_col), "NOT FOUND", y2_srs_col)))
cat(sprintf("2-year ODI: %s\n", ifelse(is.na(y2_odi_col), "NOT FOUND", y2_odi_col)))
cat("\n")

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
  cat("WARNING: Age variable not found. Will exclude from PCA.\n")
}

# Calculate T4-L1 PA if needed
if ("LATpre_L1PA" %in% names(df) && "LATpre_T4PA" %in% names(df)) {
  df$LATpre_T4_L1_PA <- df$LATpre_L1PA - df$LATpre_T4PA
  cat("Calculated LATpre_T4_L1_PA = LATpre_L1PA - LATpre_T4PA\n")
} else {
  cat("WARNING: Could not find LATpre_L1PA or LATpre_T4PA columns.\n")
  df$LATpre_T4_L1_PA <- NA
}

# Calculate T4-L1 PA for other timepoints if needed
if ("LAT6W_L1PA" %in% names(df) && "LAT6W_T4PA" %in% names(df)) {
  df$LAT6W_T4_L1_PA <- df$LAT6W_L1PA - df$LAT6W_T4PA
  cat("Calculated LAT6W_T4_L1_PA = LAT6W_L1PA - LAT6W_T4PA\n")
} else {
  df$LAT6W_T4_L1_PA <- NA
}

if ("LAT1Y_L1PA" %in% names(df) && "LAT1Y_T4PA" %in% names(df)) {
  df$LAT1Y_T4_L1_PA <- df$LAT1Y_L1PA - df$LAT1Y_T4PA
  cat("Calculated LAT1Y_T4_L1_PA = LAT1Y_L1PA - LAT1Y_T4PA\n")
} else {
  df$LAT1Y_T4_L1_PA <- NA
}

if ("LAT2Y_L1PA" %in% names(df) && "LAT2Y_T4PA" %in% names(df)) {
  df$LAT2Y_T4_L1_PA <- df$LAT2Y_L1PA - df$LAT2Y_T4PA
  cat("Calculated LAT2Y_T4_L1_PA = LAT2Y_L1PA - LAT2Y_T4PA\n")
} else {
  df$LAT2Y_T4_L1_PA <- NA
}

# Function to get confounder variable names for a given timepoint prefix
get_confounder_vars <- function(timepoint_prefix) {
  # Base confounder variables
  confounders <- c(
    paste0(timepoint_prefix, "_S1PI"),
    paste0(timepoint_prefix, "_L1_S1"),
    paste0(timepoint_prefix, "_L4_S1"),
    paste0(timepoint_prefix, "_T2_T12"),
    paste0(timepoint_prefix, "_S1PT"),
    paste0(timepoint_prefix, "_SVA_C2_S1")
  )
  
  # Add T4-L1 PA (use calculated version)
  if (timepoint_prefix == "LATpre") {
    confounders <- c(confounders, "LATpre_T4_L1_PA")
  } else if (timepoint_prefix == "LAT6W") {
    confounders <- c(confounders, "LAT6W_T4_L1_PA")
  } else if (timepoint_prefix == "LAT1Y") {
    confounders <- c(confounders, "LAT1Y_T4_L1_PA")
  } else if (timepoint_prefix == "LAT2Y") {
    confounders <- c(confounders, "LAT2Y_T4_L1_PA")
  } else {
    confounders <- c(confounders, paste0(timepoint_prefix, "_T4_L1_PA"))
  }
  
  # Add age if available (age doesn't have timepoint prefix)
  if (!is.null(age_var)) {
    confounders <- c(confounders, age_var)
  }
  
  return(confounders)
}

# Function to perform PCA-adjusted regression
perform_pca_adjusted_regression <- function(data, x_var, y_var, confounder_vars, x_label, y_label) {
  # Extract variables
  x_vec <- as.numeric(data[[x_var]])
  y_vec <- as.numeric(data[[y_var]])
  
  # Check which confounders are available
  available_confounders <- confounder_vars[confounder_vars %in% names(data)]
  
  # Filter to complete cases for all variables
  complete_cases <- !is.na(x_vec) & !is.na(y_vec)
  for (var in available_confounders) {
    complete_cases <- complete_cases & !is.na(data[[var]])
  }
  
  x_clean <- abs(x_vec[complete_cases])  # Absolute value
  y_clean <- y_vec[complete_cases]
  
  if (length(x_clean) < 3) {
    return(NULL)
  }
  
  # Prepare confounder matrix
  X_confounders <- as.matrix(data[complete_cases, available_confounders, drop = FALSE])
  
  # Need at least as many observations as confounders for PCA
  if (nrow(X_confounders) < ncol(X_confounders) + 2) {
    return(NULL)  # Insufficient data for PCA
  }
  
  # Standardize confounders
  X_scaled <- scale(X_confounders)
  
  # Perform PCA
  pca_result <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  pca_components <- pca_result$x
  eigenvalues <- pca_result$sdev^2
  
  # Determine number of components using BIC
  max_components <- min(ncol(X_confounders), nrow(X_confounders) - 2)
  if (max_components < 1) {
    return(NULL)
  }
  
  # Calculate BIC for different numbers of components
  best_bic <- Inf
  best_n_comp <- 1
  best_model <- NULL
  
  for (n_comp in 1:max_components) {
    # Create data frame with PCA components
    data_temp <- data.frame(
      y = y_clean,
      x = x_clean,
      stringsAsFactors = FALSE
    )
    for (i in 1:n_comp) {
      data_temp[[paste0("PC", i)]] <- pca_components[, i]
    }
    
    # Fit model with PCA components
    formula_str <- paste("y ~ x +", paste(paste0("PC", 1:n_comp), collapse = " + "))
    model_temp <- lm(as.formula(formula_str), data = data_temp)
    bic_temp <- BIC(model_temp)
    
    if (bic_temp < best_bic) {
      best_bic <- bic_temp
      best_n_comp <- n_comp
      best_model <- model_temp
    }
  }
  
  if (is.null(best_model)) {
    return(NULL)
  }
  
  summary_best <- summary(best_model)
  
  # Extract knee flexion coefficient (second coefficient, after intercept)
  kf_coef_idx <- 2
  slope_adj <- as.numeric(summary_best$coefficients[kf_coef_idx, 1])
  p_value_adj <- as.numeric(summary_best$coefficients[kf_coef_idx, 4])
  
  # Calculate partial correlation
  # Fit model without knee flexion
  if (best_n_comp > 0) {
    formula_no_kf <- paste("y ~", paste(paste0("PC", 1:best_n_comp), collapse = " + "))
    model_no_kf <- lm(as.formula(formula_no_kf), data = data_temp)
    r2_no_kf <- summary(model_no_kf)$r.squared
    r2_full <- summary_best$r.squared
    r2_partial <- r2_full - r2_no_kf
    
    # Partial correlation
    t_stat <- summary_best$coefficients[kf_coef_idx, 3]
    df_residual <- summary_best$df[2]
    r_partial <- t_stat / sqrt(t_stat^2 + df_residual)
  } else {
    r2_partial <- summary_best$r.squared
    r_partial <- sqrt(r2_partial) * sign(slope_adj)
  }
  
  return(list(
    n = length(x_clean),
    r_partial = r_partial,
    r2_partial = r2_partial,
    p_value = p_value_adj,
    slope = slope_adj,
    intercept = as.numeric(summary_best$coefficients[1, 1]),
    n_components = best_n_comp
  ))
}

# Function to perform linear regression and create plot
perform_regression <- function(data, x_var, y_var, x_label, y_label, title, filename, confounder_vars = NULL) {
  # Extract and convert to numeric vectors explicitly to avoid difftime issues
  x_vec <- as.numeric(data[[x_var]])
  y_vec <- as.numeric(data[[y_var]])
  
  # Filter to complete cases
  complete_cases <- !is.na(x_vec) & !is.na(y_vec)
  # Use absolute value of knee flexion
  x_clean <- abs(x_vec[complete_cases])
  y_clean <- y_vec[complete_cases]
  
  if (length(x_clean) < 3) {
    cat(sprintf("WARNING: Insufficient data for %s vs %s (n = %d)\n", x_label, y_label, length(x_clean)))
    return(NULL)
  }
  
  # Create a clean data frame with only numeric columns
  data_clean <- data.frame(x = x_clean, y = y_clean)
  names(data_clean) <- c(x_var, y_var)
  
  # Perform linear regression
  model <- lm(as.formula(paste(y_var, "~", x_var)), data = data_clean)
  summary_model <- summary(model)
  r_squared <- as.numeric(summary_model$r.squared)
  p_value <- as.numeric(summary_model$coefficients[2, 4])
  intercept <- as.numeric(summary_model$coefficients[1, 1])
  slope <- as.numeric(summary_model$coefficients[2, 1])
  
  # Calculate Pearson correlation using numeric vectors
  pearson_r <- as.numeric(cor(x_clean, y_clean, use = "complete.obs"))
  
  # Print results
  cat(sprintf("\n=== %s ===\n", title))
  cat(sprintf("Sample size: %d\n", length(x_clean)))
  cat(sprintf("Pearson r: %.4f\n", pearson_r))
  cat(sprintf("R²: %.4f\n", r_squared))
  cat(sprintf("p-value: %.4e\n", p_value))
  cat(sprintf("Slope: %.4f\n", slope))
  cat(sprintf("Intercept: %.4f\n", intercept))
  
  # Update x_label to indicate absolute value
  x_label_abs <- paste0("Absolute ", x_label)
  
  # Create scatter plot with regression line
  p <- ggplot(data_clean, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = x_label_abs,
      y = y_label,
      title = title
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Add statistics annotation
  r_text <- paste0("r = ", round(pearson_r, 3))
  r2_text <- paste0("R² = ", round(r_squared, 3))
  p_text <- paste0("p = ", formatC(p_value, format = "e", digits = 2))
  n_text <- paste0("n = ", length(x_clean))
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
  
  # Perform PCA-adjusted regression if confounders provided
  pca_result <- NULL
  if (!is.null(confounder_vars)) {
    pca_result <- perform_pca_adjusted_regression(data, x_var, y_var, confounder_vars, x_label, y_label)
    if (!is.null(pca_result)) {
      cat(sprintf("PCA-Adjusted Results:\n"))
      cat(sprintf("  Partial r: %.4f\n", pca_result$r_partial))
      cat(sprintf("  Partial R²: %.4f\n", pca_result$r2_partial))
      cat(sprintf("  p-value (adjusted): %.4e\n", pca_result$p_value))
      cat(sprintf("  Slope (adjusted): %.4f\n", pca_result$slope))
      cat(sprintf("  PCA components used: %d\n", pca_result$n_components))
    } else {
      cat("WARNING: Could not perform PCA adjustment (insufficient data or confounders)\n")
    }
  }
  
  # Return summary (ensure all values are numeric)
  result <- list(
    n = as.integer(length(x_clean)),
    r = as.numeric(pearson_r),
    r_squared = as.numeric(r_squared),
    p_value = as.numeric(p_value),
    slope = as.numeric(slope),
    intercept = as.numeric(intercept)
  )
  
  # Add PCA-adjusted results if available
  if (!is.null(pca_result)) {
    result$r_partial <- as.numeric(pca_result$r_partial)
    result$r2_partial <- as.numeric(pca_result$r2_partial)
    result$p_value_adj <- as.numeric(pca_result$p_value)
    result$slope_adj <- as.numeric(pca_result$slope)
    result$n_components <- as.integer(pca_result$n_components)
  }
  
  return(result)
}

# Perform regressions for each PROM at each timepoint
results <- list()

# Preoperative PROMs - use preoperative confounders
pre_confounders <- get_confounder_vars("LATpre")
if (!is.na(pre_srs_col) && pre_srs_col %in% names(df)) {
  results[["pre_SRS"]] <- perform_regression(
    df, 
    "LATpre_LL_KneeAngle", 
    pre_srs_col,
    "Preoperative Knee Flexion (degrees)",
    "Preoperative SRS-22 Total Score",
    "Absolute Preoperative Knee Flexion vs Preoperative SRS-22 Total Score",
    "results/analysis9_preop_KF_vs_preop_SRS.png",
    confounder_vars = pre_confounders
  )
}

if (!is.na(pre_odi_col) && pre_odi_col %in% names(df)) {
  results[["pre_ODI"]] <- perform_regression(
    df, 
    "LATpre_LL_KneeAngle", 
    pre_odi_col,
    "Preoperative Knee Flexion (degrees)",
    "Preoperative ODI Score",
    "Absolute Preoperative Knee Flexion vs Preoperative ODI Score",
    "results/analysis9_preop_KF_vs_preop_ODI.png",
    confounder_vars = pre_confounders
  )
}

# 6-week PROMs - use 6-week knee flexion and 6-week confounders
w6_confounders <- get_confounder_vars("LAT6W")
if (!is.na(w6_srs_col) && w6_srs_col %in% names(df)) {
  results[["6W_SRS"]] <- perform_regression(
    df, 
    "LAT6W_LL_KneeAngle", 
    w6_srs_col,
    "6-Week Knee Flexion (degrees)",
    "6-Week SRS-22 Total Score",
    "Absolute 6-Week Knee Flexion vs 6-Week SRS-22 Total Score",
    "results/analysis9_6W_KF_vs_6W_SRS.png",
    confounder_vars = w6_confounders
  )
}

if (!is.na(w6_odi_col) && w6_odi_col %in% names(df)) {
  results[["6W_ODI"]] <- perform_regression(
    df, 
    "LAT6W_LL_KneeAngle", 
    w6_odi_col,
    "6-Week Knee Flexion (degrees)",
    "6-Week ODI Score",
    "Absolute 6-Week Knee Flexion vs 6-Week ODI Score",
    "results/analysis9_6W_KF_vs_6W_ODI.png",
    confounder_vars = w6_confounders
  )
}

# 2-year PROMs - use 2-year knee flexion (check for LAT2Y_LL_KneeAngle or LAT1Y_LL_KneeAngle)
# First check which 2-year knee flexion column exists
y2_kf_col <- NULL
if ("LAT2Y_LL_KneeAngle" %in% names(df)) {
  y2_kf_col <- "LAT2Y_LL_KneeAngle"
  y2_kf_label <- "2-Year Knee Flexion (degrees)"
} else if ("LAT1Y_LL_KneeAngle" %in% names(df)) {
  y2_kf_col <- "LAT1Y_LL_KneeAngle"
  y2_kf_label <- "1-Year Knee Flexion (degrees)"
} else {
  cat("WARNING: Could not find 2-year or 1-year knee flexion column. Skipping 2-year analyses.\n")
}

# 2-year PROMs - use 2-year confounders (or 1-year if 2-year not available)
if (!is.null(y2_kf_col)) {
  if (y2_kf_col == "LAT2Y_LL_KneeAngle") {
    y2_confounders <- get_confounder_vars("LAT2Y")
  } else {
    y2_confounders <- get_confounder_vars("LAT1Y")
  }
  
  if (!is.na(y2_srs_col) && y2_srs_col %in% names(df)) {
    results[["2Y_SRS"]] <- perform_regression(
      df, 
      y2_kf_col, 
      y2_srs_col,
      y2_kf_label,
      "2-Year SRS-22 Total Score",
      paste("Absolute", y2_kf_label, "vs 2-Year SRS-22 Total Score"),
      "results/analysis9_2Y_KF_vs_2Y_SRS.png",
      confounder_vars = y2_confounders
    )
  }
  
  if (!is.na(y2_odi_col) && y2_odi_col %in% names(df)) {
    results[["2Y_ODI"]] <- perform_regression(
      df, 
      y2_kf_col, 
      y2_odi_col,
      y2_kf_label,
      "2-Year ODI Score",
      paste("Absolute", y2_kf_label, "vs 2-Year ODI Score"),
      "results/analysis9_2Y_KF_vs_2Y_ODI.png",
      confounder_vars = y2_confounders
    )
  }
}

# Create summary table
cat("\n\n=== SUMMARY OF ALL REGRESSIONS ===\n")
summary_rows <- list()

for (name in names(results)) {
  if (!is.null(results[[name]])) {
    # Parse name to get timepoint and PROM type
    if (grepl("^pre_", name)) {
      timepoint <- "Preoperative"
      prom_type <- ifelse(grepl("SRS", name), "SRS-22", "ODI")
    } else if (grepl("^6W_", name)) {
      timepoint <- "6 Weeks"
      prom_type <- ifelse(grepl("SRS", name), "SRS-22", "ODI")
    } else if (grepl("^2Y_", name)) {
      timepoint <- "2 Years"
      prom_type <- ifelse(grepl("SRS", name), "SRS-22", "ODI")
    } else {
      timepoint <- "Unknown"
      prom_type <- "Unknown"
    }
    
    # Determine outcome label based on timepoint (using absolute value)
    if (grepl("^pre_", name)) {
      outcome_label <- "Absolute Preoperative Knee Flexion"
    } else if (grepl("^6W_", name)) {
      outcome_label <- "Absolute 6-Week Knee Flexion"
    } else if (grepl("^2Y_", name)) {
      outcome_label <- ifelse(!is.null(y2_kf_col) && y2_kf_col == "LAT2Y_LL_KneeAngle", 
                              "Absolute 2-Year Knee Flexion", "Absolute 1-Year Knee Flexion")
    } else {
      outcome_label <- "Absolute Knee Flexion"
    }
    
    # Build row as a list to avoid rbind issues
    row_data <- list(
      Outcome = outcome_label,
      Timepoint = timepoint,
      PROM = prom_type,
      n = as.integer(results[[name]]$n),
      r_unadjusted = as.numeric(results[[name]]$r),
      R_squared_unadjusted = as.numeric(results[[name]]$r_squared),
      p_value_unadjusted = as.numeric(results[[name]]$p_value),
      Slope_unadjusted = as.numeric(results[[name]]$slope),
      Intercept_unadjusted = as.numeric(results[[name]]$intercept)
    )
    
    # Add PCA-adjusted results if available
    if ("r_partial" %in% names(results[[name]])) {
      row_data$r_partial <- as.numeric(results[[name]]$r_partial)
      row_data$R_squared_partial <- as.numeric(results[[name]]$r2_partial)
      row_data$p_value_adjusted <- as.numeric(results[[name]]$p_value_adj)
      row_data$Slope_adjusted <- as.numeric(results[[name]]$slope_adj)
      row_data$n_PCA_components <- as.integer(results[[name]]$n_components)
    } else {
      row_data$r_partial <- NA
      row_data$R_squared_partial <- NA
      row_data$p_value_adjusted <- NA
      row_data$Slope_adjusted <- NA
      row_data$n_PCA_components <- NA
    }
    
    summary_rows[[length(summary_rows) + 1]] <- row_data
  }
}

# Convert list of rows to data frame
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, lapply(summary_rows, function(x) {
    data.frame(x, stringsAsFactors = FALSE)
  }))
} else {
  summary_df <- data.frame(
    Outcome = character(),
    Timepoint = character(),
    PROM = character(),
    n = integer(),
    r_unadjusted = numeric(),
    R_squared_unadjusted = numeric(),
    p_value_unadjusted = numeric(),
    Slope_unadjusted = numeric(),
    Intercept_unadjusted = numeric(),
    r_partial = numeric(),
    R_squared_partial = numeric(),
    p_value_adjusted = numeric(),
    Slope_adjusted = numeric(),
    n_PCA_components = integer(),
    stringsAsFactors = FALSE
  )
}

# Print summary table (convert to character to avoid difftime issues)
if (nrow(summary_df) > 0) {
  # Convert all columns to appropriate types to avoid difftime issues
  summary_df$n <- as.integer(summary_df$n)
  summary_df$r_unadjusted <- as.numeric(summary_df$r_unadjusted)
  summary_df$R_squared_unadjusted <- as.numeric(summary_df$R_squared_unadjusted)
  summary_df$p_value_unadjusted <- as.numeric(summary_df$p_value_unadjusted)
  summary_df$Slope_unadjusted <- as.numeric(summary_df$Slope_unadjusted)
  summary_df$Intercept_unadjusted <- as.numeric(summary_df$Intercept_unadjusted)
  summary_df$r_partial <- as.numeric(summary_df$r_partial)
  summary_df$R_squared_partial <- as.numeric(summary_df$R_squared_partial)
  summary_df$p_value_adjusted <- as.numeric(summary_df$p_value_adjusted)
  summary_df$Slope_adjusted <- as.numeric(summary_df$Slope_adjusted)
  summary_df$n_PCA_components <- as.integer(summary_df$n_PCA_components)
  
  # Apply Bonferroni correction for multiple comparisons (both unadjusted and adjusted)
  n_tests <- nrow(summary_df)
  summary_df$p_value_unadjusted_corrected <- pmin(summary_df$p_value_unadjusted * n_tests, 1.0)
  summary_df$p_value_adjusted_corrected <- pmin(summary_df$p_value_adjusted * n_tests, 1.0)
  
  cat(sprintf("\n=== Multiple Comparison Correction ===\n"))
  cat(sprintf("Number of tests: %d\n", n_tests))
  cat(sprintf("Bonferroni corrected alpha: %.4f (0.05 / %d)\n", 0.05 / n_tests, n_tests))
  cat(sprintf("Uncorrected alpha: 0.05\n"))
  cat(sprintf("Note: Plots show uncorrected p-values. Summary table shows both uncorrected and Bonferroni-corrected p-values.\n\n"))
  
  # Save summary to CSV first (safer than print)
  write.csv(summary_df, "results/analysis9_PROM_KF_correlations_summary.csv", row.names = FALSE)
  cat("Saved summary table to results/analysis9_PROM_KF_correlations_summary.csv\n")
  
  # Print summary table manually to avoid difftime issues
  cat("\nSummary Table (Unadjusted and PCA-Adjusted with Bonferroni correction):\n")
  cat("\n=== Unadjusted Results ===\n")
  cat(sprintf("%-30s %-15s %-10s %-6s %-10s %-12s %-12s %-12s\n", 
              "Outcome", "Timepoint", "PROM", "n", "r", "R²", "p_uncorr", "p_corr"))
  cat(paste(rep("-", 120), collapse = ""), "\n")
  for (i in 1:nrow(summary_df)) {
    sig_uncorr <- ifelse(!is.na(summary_df$p_value_unadjusted[i]) && summary_df$p_value_unadjusted[i] < 0.05, "*", "")
    sig_corr <- ifelse(!is.na(summary_df$p_value_unadjusted_corrected[i]) && summary_df$p_value_unadjusted_corrected[i] < 0.05, "*", "")
    
    cat(sprintf("%-30s %-15s %-10s %-6d %-10.4f %-12.4f %-12.4e%s %-12.4e%s\n",
                summary_df$Outcome[i],
                summary_df$Timepoint[i],
                summary_df$PROM[i],
                summary_df$n[i],
                summary_df$r_unadjusted[i],
                summary_df$R_squared_unadjusted[i],
                ifelse(is.na(summary_df$p_value_unadjusted[i]), NA, summary_df$p_value_unadjusted[i]),
                sig_uncorr,
                ifelse(is.na(summary_df$p_value_unadjusted_corrected[i]), NA, summary_df$p_value_unadjusted_corrected[i]),
                sig_corr))
  }
  
  cat("\n=== PCA-Adjusted Results ===\n")
  cat(sprintf("%-30s %-15s %-10s %-6s %-10s %-12s %-12s %-12s %-6s\n", 
              "Outcome", "Timepoint", "PROM", "n", "r_partial", "R²_partial", "p_adj_uncorr", "p_adj_corr", "n_PC"))
  cat(paste(rep("-", 130), collapse = ""), "\n")
  for (i in 1:nrow(summary_df)) {
    sig_adj_uncorr <- ifelse(!is.na(summary_df$p_value_adjusted[i]) && summary_df$p_value_adjusted[i] < 0.05, "*", "")
    sig_adj_corr <- ifelse(!is.na(summary_df$p_value_adjusted_corrected[i]) && summary_df$p_value_adjusted_corrected[i] < 0.05, "*", "")
    
    cat(sprintf("%-30s %-15s %-10s %-6d %-10.4f %-12.4f %-12.4e%s %-12.4e%s %-6d\n",
                summary_df$Outcome[i],
                summary_df$Timepoint[i],
                summary_df$PROM[i],
                summary_df$n[i],
                ifelse(is.na(summary_df$r_partial[i]), NA, summary_df$r_partial[i]),
                ifelse(is.na(summary_df$R_squared_partial[i]), NA, summary_df$R_squared_partial[i]),
                ifelse(is.na(summary_df$p_value_adjusted[i]), NA, summary_df$p_value_adjusted[i]),
                sig_adj_uncorr,
                ifelse(is.na(summary_df$p_value_adjusted_corrected[i]), NA, summary_df$p_value_adjusted_corrected[i]),
                sig_adj_corr,
                ifelse(is.na(summary_df$n_PCA_components[i]), 0, summary_df$n_PCA_components[i])))
  }
  cat("\n* indicates p < 0.05\n")
} else {
  cat("\nNo results to summarize.\n")
}

cat("\nAnalysis complete!\n")
