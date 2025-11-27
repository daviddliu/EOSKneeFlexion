#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration: Toggle PJK exclusion on/off
EXCLUDE_PJK <- TRUE

# Check if cache exists - if so, skip database loading
cache_file <- "results/analysis3_PCA_significance_tracking_cache.rds"
if (file.exists(cache_file)) {
  cat("\n=== Cache found - skipping database load ===\n")
  load_database <- FALSE
} else {
  cat("\n=== No cache found - loading database ===\n")
  load_database <- TRUE
}

if (load_database) {
  # Load database
  db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective projects/Databases/CADS database - 2025.10.10.xlsx"
  df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

  # Impute missing surgeon by randomly assigning to other surgeons at the same site
  # Set random seed for reproducibility
  set.seed(12345)

  if ("demo_site_txt" %in% names(df)) {
    n_before <- sum(is.na(df$SxI_Gen_Surgeon))
    
    # Convert surgeon to character for easier manipulation
    df$SxI_Gen_Surgeon <- as.character(df$SxI_Gen_Surgeon)
    
    # Get all sites with patients that have NA surgeon
    sites_with_na <- unique(df$demo_site_txt[is.na(df$SxI_Gen_Surgeon) & !is.na(df$demo_site_txt)])
    
    n_assigned <- 0
    n_fallback <- 0
    n_unknown <- 0
    
    # For each site, assign patients with NA surgeon to other surgeons at that site
    for (site in sites_with_na) {
      # Get all surgeons at this site (excluding NA)
      surgeons_at_site <- unique(df$SxI_Gen_Surgeon[df$demo_site_txt == site & !is.na(df$SxI_Gen_Surgeon)])
      
      # Get indices of patients at this site with NA surgeon
      idx_na_at_site <- which(is.na(df$SxI_Gen_Surgeon) & df$demo_site_txt == site)
      
      if (length(surgeons_at_site) > 0) {
        # Randomly assign to one of the surgeons at this site
        assigned_surgeons <- sample(surgeons_at_site, length(idx_na_at_site), replace = TRUE)
        df$SxI_Gen_Surgeon[idx_na_at_site] <- assigned_surgeons
        n_assigned <- n_assigned + length(idx_na_at_site)
      } else {
        # No surgeons at this site - use site-based grouping as fallback
        df$SxI_Gen_Surgeon[idx_na_at_site] <- paste0("Site_", site)
        n_fallback <- n_fallback + length(idx_na_at_site)
      }
    }
    
    # Handle patients with NA surgeon and NA site (use "Unknown")
    idx_na_no_site <- which(is.na(df$SxI_Gen_Surgeon))
    if (length(idx_na_no_site) > 0) {
      df$SxI_Gen_Surgeon[idx_na_no_site] <- "Unknown"
      n_unknown <- length(idx_na_no_site)
    }
    
    cat(paste("Imputed", n_before, "patients with missing surgeon:\n"))
    cat(paste("  - Assigned to other surgeons at same site:", n_assigned, "\n"))
    cat(paste("  - Assigned to site-based group (no other surgeons):", n_fallback, "\n"))
    if (n_unknown > 0) {
      cat(paste("  - Assigned to 'Unknown' (no site data):", n_unknown, "\n"))
    }
  } else {
    cat("Warning: demo_site_txt column not found. Patients with missing surgeon will be excluded.\n")
  }

  # Get unique surgeons (now includes imputed site-based groups)
  surgeons <- unique(df$SxI_Gen_Surgeon)
  surgeons <- surgeons[!is.na(surgeons)]
  cat(paste("Found", length(surgeons), "surgeon/site groups\n"))
}

# Function for Analysis 1: Preop knee flexion vs change in lordosis
analysis1_by_surgeon <- function(data, surgeon_name) {
  # Calculate change in lordosis using 6-week data
  data$change_lordosis <- data$LAT6W_L1_S1 - data$LATpre_L1_S1
  
  # Drop patients with missing data (listwise deletion)
  # No imputation - patients with missing LATpre_LL_KneeAngle or change_lordosis are excluded
  data_clean <- data %>%
    filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis))
  
  if (nrow(data_clean) == 0) {
    cat(paste("Warning: No valid data for surgeon", surgeon_name, "in Analysis 1\n"))
    return(NULL)
  }
  
  # Check if we have enough data points for regression (need at least 2 points)
  if (nrow(data_clean) < 2) {
    cat(paste("Warning: Insufficient data (n < 2) for surgeon", surgeon_name, "in Analysis 1\n"))
    return(NULL)
  }
  
  # Perform linear regression
  model <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = data_clean)
  summary_model <- summary(model)
  
  # Check if regression coefficients are available
  if (nrow(summary_model$coefficients) < 2) {
    cat(paste("Warning: Regression failed for surgeon", surgeon_name, "in Analysis 1\n"))
    return(NULL)
  }
  
  r_squared <- summary_model$r.squared
  p_value <- summary_model$coefficients[2, 4]
  
  # Calculate Pearson correlation coefficient
  pearson_r <- cor(data_clean$LATpre_LL_KneeAngle, data_clean$change_lordosis, use = "complete.obs")
  
  # Create plot
  p <- ggplot(data_clean, aes(x = LATpre_LL_KneeAngle, y = change_lordosis)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
      y = "Change in Lordosis (LAT6W_L1_S1 - LATpre_L1_S1, 6-week)",
      title = paste0("Preoperative Knee Flexion vs Change in Lordosis (6-week)\nSurgeon: ", surgeon_name),
      subtitle = paste0("n = ", nrow(data_clean))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
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
  if (!dir.exists("results/by_surgeon")) {
    dir.create("results/by_surgeon")
  }
  safe_surgeon_name <- gsub("[^A-Za-z0-9_]", "_", surgeon_name)
  filename <- paste0("analysis3_surgeon_", safe_surgeon_name, "_analysis1.png")
  filepath <- file.path("results/by_surgeon", filename)
  ggsave(filepath, plot = p, width = 10, height = 8, dpi = 300)
  cat(paste("Saved Analysis 1 for surgeon", surgeon_name, "to", filepath, "\n"))
  
  return(p)
}

# Function for Analysis 2: PI-LL mismatch vs knee flexion measures
analysis2_by_surgeon <- function(data, surgeon_name) {
  # Drop patients with missing PI-LL data (using 6-week data)
  data_filtered <- data %>%
    filter(!is.na(LAT6W_PI_LL))
  
  if (nrow(data_filtered) == 0) {
    cat(paste("Warning: No patients with PI-LL data for surgeon", surgeon_name, "\n"))
    return(NULL)
  }
  
  # Calculate difference in knee angle (using 6-week data)
  # Note: If either value is missing, knee_angle_diff will be NA (dropped in each plot)
  data_filtered$knee_angle_diff <- data_filtered$LAT6W_LL_KneeAngle - data_filtered$LATpre_LL_KneeAngle
  
  safe_surgeon_name <- gsub("[^A-Za-z0-9_]", "_", surgeon_name)
  
  # Helper function to create plot
  create_plot <- function(data, x_var, y_var, x_label, y_label, title_suffix, plot_num) {
    # Drop patients with missing data for this specific variable pair (listwise deletion)
    # No imputation - patients with missing values are excluded
    data_clean <- data %>%
      filter(!is.na(!!sym(x_var)) & !is.na(!!sym(y_var)))
    
    if (nrow(data_clean) == 0) {
      cat(paste("Warning: No valid data for surgeon", surgeon_name, "in Analysis 2, Plot", plot_num, "\n"))
      return(NULL)
    }
    
    # Check if we have enough data points for regression (need at least 2 points)
    if (nrow(data_clean) < 2) {
      cat(paste("Warning: Insufficient data (n < 2) for surgeon", surgeon_name, "in Analysis 2, Plot", plot_num, "\n"))
      return(NULL)
    }
    
    # Perform linear regression
    formula <- as.formula(paste(y_var, "~", x_var))
    model <- lm(formula, data = data_clean)
    summary_model <- summary(model)
    
    # Check if regression coefficients are available
    if (nrow(summary_model$coefficients) < 2) {
      cat(paste("Warning: Regression failed for surgeon", surgeon_name, "in Analysis 2, Plot", plot_num, "\n"))
      return(NULL)
    }
    
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
        title = paste0(title_suffix, " vs PI-LL Mismatch (6-week)\nSurgeon: ", surgeon_name),
        subtitle = paste0("n = ", nrow(data_clean))
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
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
    if (!dir.exists("results/by_surgeon")) {
      dir.create("results/by_surgeon")
    }
    filename <- paste0("analysis3_surgeon_", safe_surgeon_name, "_analysis2_", plot_num, ".png")
    filepath <- file.path("results/by_surgeon", filename)
    ggsave(filepath, plot = p, width = 10, height = 8, dpi = 300)
    cat(paste("Saved Analysis 2, Plot", plot_num, "for surgeon", surgeon_name, "to", filepath, "\n"))
    
    return(p)
  }
  
  # Plot 1: Preop Knee Angle vs PI-LL
  create_plot(
    data_filtered,
    "LATpre_LL_KneeAngle",
    "LAT6W_PI_LL",
    "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
    "PI-LL Mismatch (6-week)",
    "Preoperative Knee Flexion",
    "1"
  )
  
  # Plot 2: Postop Knee Angle vs PI-LL
  create_plot(
    data_filtered,
    "LAT6W_LL_KneeAngle",
    "LAT6W_PI_LL",
    "Postoperative Knee Flexion (6-week)",
    "PI-LL Mismatch (6-week)",
    "Postoperative Knee Flexion (6-week)",
    "2"
  )
  
  # Plot 3: Change in Knee Angle vs PI-LL
  create_plot(
    data_filtered,
    "knee_angle_diff",
    "LAT6W_PI_LL",
    "Change in Knee Flexion (6-week)",
    "PI-LL Mismatch (6-week)",
    "Change in Knee Flexion (6-week)",
    "3"
  )
}

# Function for Analysis 4: PCA-based confounder analysis (preop knee flexion adjusted for PCA components)
analysis4_by_surgeon <- function(data, surgeon_name) {
  # Calculate change in lordosis using 6-week data
  data$change_lordosis <- data$LAT6W_L1_S1 - data$LATpre_L1_S1
  
  # Calculate T4-L1 PA
  if ("LATpre_L1PA" %in% names(data) && "LATpre_T4PA" %in% names(data)) {
    data$LATpre_T4_L1_PA <- data$LATpre_L1PA - data$LATpre_T4PA
  } else {
    data$LATpre_T4_L1_PA <- NA
  }
  
  # Check for age variable (prioritize demo_Age)
  age_var <- NULL
  if ("demo_Age" %in% names(data)) {
    age_var <- "demo_Age"
  } else if ("demo_age" %in% names(data)) {
    age_var <- "demo_age"
  } else if ("age" %in% names(data)) {
    age_var <- "age"
  } else if ("Age" %in% names(data)) {
    age_var <- "Age"
  }
  
  # Drop patients with missing data for confounder analysis
  # Confounder set: S1PI (preop), Preop Lordosis (L1-S1), L4-S1, Thoracic Kyphosis, S1PT, SVA, T4-L1 PA, Age
  if (!is.null(age_var)) {
    data_clean <- data %>%
      filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
             !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) & 
             !is.na(LATpre_S1PT) & !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA) & 
             !is.na(.data[[age_var]]))
  } else {
    data_clean <- data %>%
      filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis) & !is.na(LATpre_S1PI) &
             !is.na(LATpre_L1_S1) & !is.na(LATpre_L4_S1) & !is.na(LATpre_T2_T12) & 
             !is.na(LATpre_S1PT) & !is.na(LATpre_SVA_C2_S1) & !is.na(LATpre_T4_L1_PA))
  }
  
  n_cases <- nrow(data_clean)
  
  if (n_cases == 0) {
    cat(paste("Warning: No valid data for surgeon", surgeon_name, "in Analysis 4\n"))
    return(list(surgeon = surgeon_name, n_cases = 0, pval_unadjusted = NA, pval_adjusted = NA, r2_adjusted = NA))
  }
  
  # Check if we have enough data points for regression
  if (n_cases < 2) {
    cat(paste("Warning: Insufficient data (n < 2) for surgeon", surgeon_name, "in Analysis 4\n"))
    return(list(surgeon = surgeon_name, n_cases = n_cases, pval_unadjusted = NA, pval_adjusted = NA, r2_adjusted = NA))
  }
  
  # Model 1: Simple regression (preop knee flexion only)
  model1 <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = data_clean)
  summary1 <- summary(model1)
  
  if (nrow(summary1$coefficients) < 2) {
    cat(paste("Warning: Regression failed for surgeon", surgeon_name, "in Analysis 4, Model 1\n"))
    return(list(surgeon = surgeon_name, n_cases = n_cases, pval_unadjusted = NA, pval_adjusted = NA, r2_adjusted = NA))
  }
  
  coef1_kf <- summary1$coefficients[2, 1]
  pval1_kf <- summary1$coefficients[2, 4]
  r2_model1 <- summary1$r.squared
  
  # Prepare confounder variables for PCA
  confounder_vars <- c("LATpre_S1PI", "LATpre_L1_S1", "LATpre_L4_S1", "LATpre_T2_T12", 
                       "LATpre_S1PT", "LATpre_SVA_C2_S1", "LATpre_T4_L1_PA")
  
  # Add age if available
  if (!is.null(age_var)) {
    confounder_vars <- c(confounder_vars, age_var)
  }
  
  X_confounders <- data_clean[, confounder_vars]
  
  # Check if we have enough data for PCA (need at least as many observations as variables)
  n_confounders <- length(confounder_vars)
  if (n_cases < n_confounders) {
    cat(paste("Warning: Insufficient data for PCA (n =", n_cases, "<", n_confounders, "variables) for surgeon", surgeon_name, "\n"))
    cat(paste("  Falling back to original confounder model\n"))
    # Fall back to original model
    if (!is.null(age_var)) {
      formula_str <- paste("change_lordosis ~ LATpre_LL_KneeAngle + LATpre_S1PI + LATpre_L1_S1 +",
                          "LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 +",
                          "LATpre_T4_L1_PA +", age_var)
      model2 <- lm(as.formula(formula_str), data = data_clean)
      covariates_formula_str <- paste("change_lordosis ~ LATpre_S1PI + LATpre_L1_S1 +",
                                       "LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 +",
                                       "LATpre_T4_L1_PA +", age_var)
      covariates_model <- lm(as.formula(covariates_formula_str), data = data_clean)
    } else {
      model2 <- lm(change_lordosis ~ LATpre_LL_KneeAngle + LATpre_S1PI + LATpre_L1_S1 + 
                   LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA, data = data_clean)
      covariates_model <- lm(change_lordosis ~ LATpre_S1PI + LATpre_L1_S1 + 
                             LATpre_L4_S1 + LATpre_T2_T12 + LATpre_S1PT + LATpre_SVA_C2_S1 + LATpre_T4_L1_PA, data = data_clean)
    }
    summary2 <- summary(model2)
    if (nrow(summary2$coefficients) < 2) {
      return(list(surgeon = surgeon_name, n_cases = n_cases, pval_unadjusted = pval1_kf, pval_adjusted = NA, r2_adjusted = NA))
    }
    coef2_kf <- summary2$coefficients[2, 1]
    pval2_kf <- summary2$coefficients[2, 4]
    r2_model2 <- summary2$r.squared
    data_clean$resid_from_covariates <- residuals(covariates_model)
    r2_covariates_only <- summary(covariates_model)$r.squared
    r2_kf_change <- r2_model2 - r2_covariates_only
    t_stat_kf <- summary2$coefficients[2, 3]
    df_residual <- summary2$df[2]
    r_partial <- t_stat_kf / sqrt(t_stat_kf^2 + df_residual)
    n_components_used <- "Original (insufficient data for PCA)"
  } else {
    # Perform PCA
    X_scaled <- scale(X_confounders)
    pca_result <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
    pca_components <- pca_result$x
    pca_summary <- summary(pca_result)
    variance_explained <- pca_summary$importance[2, ]
    cumulative_variance <- pca_summary$importance[3, ]
    eigenvalues <- pca_result$sdev^2
    
    # Determine number of components using BIC (similar to analysis5)
    # For small sample sizes, limit to reasonable number
    max_components_to_test <- min(n_confounders, n_cases - 2, length(cumulative_variance))  # Need at least 2 more observations than parameters
    
    if (max_components_to_test < 1) {
      cat(paste("Warning: Insufficient data for PCA model (n =", n_cases, ") for surgeon", surgeon_name, "\n"))
      return(list(surgeon = surgeon_name, n_cases = n_cases, pval_unadjusted = pval1_kf, pval_adjusted = NA, r2_adjusted = NA))
    }
    
    # Calculate BIC for different numbers of components
    ic_results <- data.frame(
      N_Components = 1:max_components_to_test,
      BIC = numeric(max_components_to_test),
      Knee_Flexion_P = numeric(max_components_to_test)
    )
    
    for (n_comp in 1:max_components_to_test) {
      # Prepare data with n_comp components
      data_temp <- data_clean
      for (i in 1:n_comp) {
        data_temp[[paste0("PC", i)]] <- pca_components[, i]
      }
      
      # Fit model
      pca_formula_temp <- as.formula(paste("change_lordosis ~ LATpre_LL_KneeAngle +", 
                                           paste(paste0("PC", 1:n_comp), collapse = " + ")))
      model_temp <- lm(pca_formula_temp, data = data_temp)
      
      ic_results$BIC[n_comp] <- BIC(model_temp)
      ic_results$Knee_Flexion_P[n_comp] <- summary(model_temp)$coefficients[2, 4]
    }
    
    # Use BIC optimal
    optimal_bic_idx <- which.min(ic_results$BIC)
    n_components_use <- ic_results$N_Components[optimal_bic_idx]
    
    # Prepare data with optimal number of components
    for (i in 1:n_components_use) {
      data_clean[[paste0("PC", i)]] <- pca_components[, i]
    }
    
    # Model 2: Multiple regression with PCA components
    pca_formula <- as.formula(paste("change_lordosis ~ LATpre_LL_KneeAngle +", 
                                    paste(paste0("PC", 1:n_components_use), collapse = " + ")))
    model2 <- lm(pca_formula, data = data_clean)
    summary2 <- summary(model2)
    
    if (nrow(summary2$coefficients) < 2) {
      cat(paste("Warning: Regression failed for surgeon", surgeon_name, "in Analysis 4, Model 2 (PCA)\n"))
      return(list(surgeon = surgeon_name, n_cases = n_cases, pval_unadjusted = pval1_kf, pval_adjusted = NA, r2_adjusted = NA))
    }
    
    coef2_kf <- summary2$coefficients[2, 1]
    pval2_kf <- summary2$coefficients[2, 4]
    r2_model2 <- summary2$r.squared
    
    # Calculate residuals from PCA covariates model (for plotting)
    pca_formula_no_kf <- as.formula(paste("change_lordosis ~", 
                                          paste(paste0("PC", 1:n_components_use), collapse = " + ")))
    covariates_model <- lm(pca_formula_no_kf, data = data_clean)
    data_clean$resid_from_covariates <- residuals(covariates_model)
    
    # Calculate R² change for knee flexion
    r2_covariates_only <- summary(covariates_model)$r.squared
    r2_kf_change <- r2_model2 - r2_covariates_only
    
    # Calculate partial correlation
    t_stat_kf <- summary2$coefficients[2, 3]
    df_residual <- summary2$df[2]
    r_partial <- t_stat_kf / sqrt(t_stat_kf^2 + df_residual)
    
    n_components_used <- paste0(n_components_use, " PCs (", round(cumulative_variance[n_components_use] * 100, 1), "% variance)")
  }
  
  # Calculate zero-order correlation
  r_kf_lordosis <- cor(data_clean$LATpre_LL_KneeAngle, data_clean$change_lordosis, use = "complete.obs")
  
  # Print summary
  cat(paste("  Analysis 4 (PCA) - Unadjusted coefficient:", round(coef1_kf, 4), 
            "(p =", formatC(pval1_kf, format = "e", digits = 2), ")\n"))
  cat(paste("  Analysis 4 (PCA) - Adjusted coefficient:", round(coef2_kf, 4), 
            "(p =", formatC(pval2_kf, format = "e", digits = 2), ")\n"))
  cat(paste("  Analysis 4 (PCA) - R² Model 1:", round(r2_model1, 3), 
            "R² Model 2:", round(r2_model2, 3), "\n"))
  cat(paste("  Analysis 4 (PCA) - Components used:", n_components_used, "\n"))
  
  # Create plot: Residuals from PCA covariates model vs Knee Flexion
  p <- ggplot(data_clean, aes(x = LATpre_LL_KneeAngle, y = resid_from_covariates)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
      x = "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
      y = paste0("Residuals from PCA Covariates Model\n(Change in Lordosis after removing effects of ", n_components_used, ")"),
      title = paste0("PCA-Based Confounder Analysis (6-week follow-up)\nSurgeon: ", surgeon_name),
      subtitle = paste0("Adjusted relationship (n = ", nrow(data_clean), ")")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("Partial r = ", round(r_partial, 3),
                     "\nAdjusted coeff = ", round(coef2_kf, 4),
                     "\np = ", formatC(pval2_kf, format = "e", digits = 2),
                     "\nR² adjusted (KF) = ", round(r2_kf_change, 3)),
      hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold"
    )
  
  # Save plot
  if (!dir.exists("results")) {
    dir.create("results")
  }
  if (!dir.exists("results/by_surgeon")) {
    dir.create("results/by_surgeon")
  }
  safe_surgeon_name <- gsub("[^A-Za-z0-9_]", "_", surgeon_name)
  filename <- paste0("analysis3_surgeon_", safe_surgeon_name, "_analysis4_PCA.png")
  filepath <- file.path("results/by_surgeon", filename)
  ggsave(filepath, plot = p, width = 10, height = 8, dpi = 300)
  cat(paste("Saved Analysis 4 (PCA) for surgeon", surgeon_name, "to", filepath, "\n"))
  
  # Return information for tracking significance changes
  return(list(surgeon = surgeon_name, n_cases = n_cases, 
              pval_unadjusted = pval1_kf, pval_adjusted = pval2_kf, r2_adjusted = r2_kf_change))
}

# Track significance changes for Analysis 4
significance_tracking <- list()

# Check if cache exists - if so, load it and skip analyses
if (file.exists(cache_file)) {
  cat("\n=== Loading cached results ===\n")
  significance_tracking_cached <- readRDS(cache_file)
  cat(paste("Loaded", length(significance_tracking_cached), "cached surgeon results\n"))
  cat("Skipping individual analyses - using cached data\n\n")
  significance_tracking <- significance_tracking_cached
  run_analyses <- FALSE
} else {
  cat("\n=== Running analyses (no cache found) ===\n")
  run_analyses <- TRUE
}

if (run_analyses) {
  if (!exists("df")) {
    stop("Error: Database not loaded but cache not found. Cannot proceed.")
  }
  # Loop through each surgeon and perform all analyses
  for (surgeon in surgeons) {
    cat(paste("\n=== Processing surgeon:", surgeon, "===\n"))
    
    # Filter data for current surgeon
    df_surgeon <- df %>%
      filter(SxI_Gen_Surgeon == surgeon)
    
    cat(paste("Surgeon", surgeon, "has", nrow(df_surgeon), "patients\n"))
    
    # Perform Analysis 1
    analysis1_by_surgeon(df_surgeon, surgeon)
    
    # Perform Analysis 2
    analysis2_by_surgeon(df_surgeon, surgeon)
    
    # Perform Analysis 4 and track results
    result_analysis4 <- analysis4_by_surgeon(df_surgeon, surgeon)
    if (!is.null(result_analysis4)) {
      significance_tracking[[length(significance_tracking) + 1]] <- result_analysis4
    }
  }
  
  # Save results to cache
  if (!dir.exists("results")) {
    dir.create("results")
  }
  saveRDS(significance_tracking, cache_file)
  cat(paste("\n=== Saved results to cache:", cache_file, "===\n"))
}

# Analyze significance changes
cat("\n=== Significance Change Analysis (Surgeons with >= 20 cases) ===\n")
if (length(significance_tracking) > 0) {
  # Convert to data frame
  tracking_df <- do.call(rbind, lapply(significance_tracking, function(x) {
    data.frame(
      surgeon = x$surgeon,
      n_cases = x$n_cases,
      pval_unadjusted = x$pval_unadjusted,
      pval_adjusted = x$pval_adjusted,
      r2_adjusted = ifelse(is.null(x$r2_adjusted), NA, x$r2_adjusted),
      stringsAsFactors = FALSE
    )
  }))
  
  # Filter for surgeons with >= 20 cases and valid p-values
  surgeons_gt20 <- tracking_df %>%
    filter(n_cases >= 20) %>%
    filter(!is.na(pval_unadjusted) & !is.na(pval_adjusted))
  
  cat(paste("Total surgeons with >= 20 cases and valid p-values:", nrow(surgeons_gt20), "\n"))
  
  # Calculate multiple comparisons correction
  n_surgeons <- nrow(surgeons_gt20)
  if (n_surgeons > 0) {
    cat(sprintf("\nSignificance Change Analysis:\n"))
    cat(sprintf("  Number of comparisons: %d\n", n_surgeons))
    cat(sprintf("  Alpha level: 0.05\n\n"))
  }
  
  if (nrow(surgeons_gt20) > 0) {
    # Categorize surgeons by significance change pattern
    # Compare unadjusted vs adjusted p-values (from multivariate model)
    stayed_sig <- surgeons_gt20 %>%
      filter(pval_unadjusted < 0.05 & pval_adjusted < 0.05)
    
    stayed_nonsig <- surgeons_gt20 %>%
      filter(pval_unadjusted >= 0.05 & pval_adjusted >= 0.05)
    
    changed_to_nonsig <- surgeons_gt20 %>%
      filter(pval_unadjusted < 0.05 & pval_adjusted >= 0.05)
    
    changed_to_sig <- surgeons_gt20 %>%
      filter(pval_unadjusted >= 0.05 & pval_adjusted < 0.05)
    
    # Print results for each category
    cat("1. Stayed Significant (both unadjusted and adjusted p < 0.05):", nrow(stayed_sig), "\n")
    if (nrow(stayed_sig) > 0) {
      for (i in seq_len(nrow(stayed_sig))) {
        cat(sprintf("   - %s: n = %d, unadjusted p = %.4e, adjusted p = %.4e\n",
                    stayed_sig$surgeon[i],
                    stayed_sig$n_cases[i],
                    stayed_sig$pval_unadjusted[i],
                    stayed_sig$pval_adjusted[i]))
      }
    }
    
    cat("\n2. Stayed Non-Significant (both unadjusted and adjusted p >= 0.05):", nrow(stayed_nonsig), "\n")
    if (nrow(stayed_nonsig) > 0) {
      for (i in seq_len(nrow(stayed_nonsig))) {
        cat(sprintf("   - %s: n = %d, unadjusted p = %.4e, adjusted p = %.4e\n",
                    stayed_nonsig$surgeon[i],
                    stayed_nonsig$n_cases[i],
                    stayed_nonsig$pval_unadjusted[i],
                    stayed_nonsig$pval_adjusted[i]))
      }
    }
    
    cat("\n3. Changed from Significant to Non-Significant (unadjusted p < 0.05, adjusted p >= 0.05):", nrow(changed_to_nonsig), "\n")
    if (nrow(changed_to_nonsig) > 0) {
      for (i in seq_len(nrow(changed_to_nonsig))) {
        cat(sprintf("   - %s: n = %d, unadjusted p = %.4e, adjusted p = %.4e\n",
                    changed_to_nonsig$surgeon[i],
                    changed_to_nonsig$n_cases[i],
                    changed_to_nonsig$pval_unadjusted[i],
                    changed_to_nonsig$pval_adjusted[i]))
      }
    }
    
    cat("\n4. Changed from Non-Significant to Significant (unadjusted p >= 0.05, adjusted p < 0.05):", nrow(changed_to_sig), "\n")
    if (nrow(changed_to_sig) > 0) {
      for (i in seq_len(nrow(changed_to_sig))) {
        cat(sprintf("   - %s: n = %d, unadjusted p = %.4e, adjusted p = %.4e\n",
                    changed_to_sig$surgeon[i],
                    changed_to_sig$n_cases[i],
                    changed_to_sig$pval_unadjusted[i],
                    changed_to_sig$pval_adjusted[i]))
      }
    }
    
    # Create formatted table for figure
    cat("\n=== Creating formatted table for figure ===\n")
    
    # Format p-values as 0.xxx (not scientific notation)
    format_pval <- function(p) {
      if (is.na(p) || is.nan(p)) {
        return("NA")
      } else if (p < 0.001) {
        return("< 0.001")
      } else {
        return(sprintf("%.3f", p))
      }
    }
    
    # Prepare data for each category
    if (!dir.exists("results")) {
      dir.create("results")
    }
    
    # Calculate total rows for height calculation
    # For empty categories, count as 1 row (for "None")
    total_rows <- sum(
      max(1, nrow(stayed_sig)), 
      max(1, nrow(stayed_nonsig)), 
      max(1, nrow(changed_to_nonsig)), 
      max(1, nrow(changed_to_sig))
    )
    # Calculate height: title + subheadings + tables + margins
    title_height <- 0.6
    subheading_height <- 0.25  # Per subheading
    header_height <- 0.3  # Per table header
    row_height <- 0.2  # Per data row
    spacing <- 0.15  # Between sections
    margin <- 0.4  # Bottom margin
    
    # Always show all 4 categories
    n_subheadings <- 4
    n_tables <- 4
    
    img_height <- title_height + 
                  (n_subheadings * subheading_height) + 
                  (n_tables * header_height) + 
                  (total_rows * row_height) + 
                  ((n_tables - 1) * spacing) + 
                  margin
    img_height <- max(14, img_height)  # Minimum 14 inches
    img_width <- 16  # Increased width for better visibility
    
    png("results/analysis3_PCA_bonferroni_table.png", width = img_width, height = img_height, units = "in", res = 300)
    grid.newpage()
    
    # Create title
    title_text <- paste0("Adjusted p-values after PCA-based confounder adjustment\nPreop knee flexion vs change in lordosis at 6 weeks (Surgeons with >=20 cases with available data)")
    
    title_grob <- textGrob(title_text, gp = gpar(fontsize = 10, fontface = "bold"),
                          x = 0.5, y = unit(1, "npc") - unit(0.2, "in"), just = "top")
    grid.draw(title_grob)
    
    # Helper function to create a category table with subheading
    create_category_table <- function(category_df, category_name, y_pos) {
      # If empty, create a table with "None"
      if (nrow(category_df) == 0) {
        table_data <- data.frame(
          "Surgeon/Site" = "None",
          "n" = "",
          "Unadjusted p" = "",
          "Adjusted p" = "",
          "R² adjusted (KF)" = "",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      } else {
        # Format the data - use adjusted p-values from multivariate model
        table_data <- category_df %>%
          mutate(
            pval_unadjusted_fmt = sapply(pval_unadjusted, format_pval),
            pval_adjusted_fmt = sapply(pval_adjusted, format_pval),
            r2_adjusted_fmt = ifelse(is.na(r2_adjusted), "", sprintf("%.3f", r2_adjusted))
          ) %>%
          select(surgeon, n_cases, pval_unadjusted_fmt, pval_adjusted_fmt, r2_adjusted_fmt) %>%
          rename(
            "Surgeon/Site" = surgeon,
            "n" = n_cases,
            "Unadjusted p" = pval_unadjusted_fmt,
            "Adjusted p" = pval_adjusted_fmt,
            "R² adjusted (KF)" = r2_adjusted_fmt
          )
      }
      
      # Create subheading - position even further to the right
      subheading_grob <- textGrob(category_name, 
                                  gp = gpar(fontsize = 9, fontface = "bold"),
                                  x = unit(6.5, "in"), y = y_pos, just = c("left", "top"))
      grid.draw(subheading_grob)
      y_pos <- y_pos - unit(0.25, "in")
      
      # Create table for this category
      table_grob <- tableGrob(table_data,
                             theme = ttheme_minimal(
                               core = list(
                                 fg_params = list(fontsize = 7.5, hjust = 0, x = 0.01),
                                 bg_params = list(fill = rep("white", nrow(table_data)))
                               ),
                               colhead = list(
                                 fg_params = list(fontsize = 8, fontface = "bold", hjust = 0, x = 0.01),
                                 bg_params = list(fill = "grey80")
                               ),
                               rowhead = list(fg_params = list(fontsize = 7.5))
                             ),
                             rows = NULL)
      
      # Adjust column widths - include R² column
      table_grob$widths <- unit(c(1.3, 0.6, 1.0, 1.0, 0.8), "in")
      
      # Adjust row heights
      table_grob$heights[1] <- unit(0.3, "in")  # Header row
      for (i in 2:length(table_grob$heights)) {
        table_grob$heights[i] <- unit(0.2, "in")  # Data rows
      }
      
      # Calculate table height
      table_height <- sum(table_grob$heights)
      
      # Draw table - position even further to the right
      pushViewport(viewport(y = y_pos, 
                           height = table_height,
                           just = "top",
                           x = unit(6.5, "in"), 
                           width = unit(img_width - 7, "in"),
                           clip = "off"))
      grid.draw(table_grob)
      popViewport()
      
      # Return new y position (move down by subheading + table + spacing)
      y_pos - table_height - unit(0.15, "in")
    }
    
    # Create tables for each category, starting below title
    # Always show all 4 categories, even if empty (will show "None")
    y_pos <- unit(1, "npc") - unit(0.8, "in")
    
    y_pos <- create_category_table(stayed_sig, 
                                   "1. Stayed Significant (both unadjusted and adjusted p < 0.05)", 
                                   y_pos)
    
    y_pos <- create_category_table(stayed_nonsig,
                                   "2. Stayed Non-Significant (both unadjusted and adjusted p >= 0.05)",
                                   y_pos)
    
    y_pos <- create_category_table(changed_to_nonsig,
                                   "3. Changed from Significant to Non-Significant (unadjusted p < 0.05, adjusted p >= 0.05)",
                                   y_pos)
    
    y_pos <- create_category_table(changed_to_sig,
                                   "4. Changed from Non-Significant to Significant (unadjusted p >= 0.05, adjusted p < 0.05)",
                                   y_pos)
    
    dev.off()
    
    cat("Saved formatted table image to results/analysis3_PCA_bonferroni_table.png\n")
  } else {
    cat("No surgeons with >= 20 cases and valid p-values found.\n")
  }
} else {
  cat("No valid Analysis 4 results to analyze.\n")
}

cat("\n=== All analyses completed ===\n")

