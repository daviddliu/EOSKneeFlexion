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

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
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

# Filter patients with satisfactory PI-LL mismatch (between -10 and 10)
# Drop patients with missing PI-LL data (using 6-week data)
df_filtered <- df %>%
  filter(LAT6W_PI_LL >= -10 & LAT6W_PI_LL <= 10, !is.na(LAT6W_PI_LL))

cat(paste("\nFiltered to", nrow(df_filtered), "patients with PI-LL between -10 and 10 (6-week)\n"))

# Calculate difference in knee angle (using 6-week data)
df_filtered$knee_angle_diff <- df_filtered$LAT6W_LL_KneeAngle - df_filtered$LATpre_LL_KneeAngle

# Get unique surgeons
surgeons <- unique(df_filtered$SxI_Gen_Surgeon)
surgeons <- surgeons[!is.na(surgeons)]
cat(paste("Found", length(surgeons), "surgeon/site groups\n\n"))

# Function to analyze corrected LL vs knee flexion for a surgeon
analysis_correctedLL_kf_by_surgeon <- function(data, surgeon_name) {
  # Filter to this surgeon's patients
  df_surgeon <- data %>%
    filter(SxI_Gen_Surgeon == surgeon_name)
  
  n_total <- nrow(df_surgeon)
  cat(paste("=== Processing surgeon:", surgeon_name, "===\n"))
  cat(paste("Surgeon", surgeon_name, "has", n_total, "patients\n"))
  
  # Analysis 1: PI-LL vs Preop Knee Angle
  data_preop <- df_surgeon %>%
    filter(!is.na(LAT6W_PI_LL) & !is.na(LATpre_LL_KneeAngle))
  
  n_preop <- nrow(data_preop)
  
  if (n_preop < 2) {
    cat(paste("Warning: Insufficient data (n < 2) for surgeon", surgeon_name, "in Preop analysis\n"))
    return(list(
      surgeon = surgeon_name,
      n_total = n_total,
      n_preop = n_preop,
      n_change = 0,
      preop_r2 = NA,
      preop_p = NA,
      change_r2 = NA,
      change_p = NA
    ))
  }
  
  # Regression: PI-LL vs Preop Knee Angle
  model_preop <- lm(LATpre_LL_KneeAngle ~ LAT6W_PI_LL, data = data_preop)
  summary_preop <- summary(model_preop)
  
  if (nrow(summary_preop$coefficients) < 2) {
    cat(paste("Warning: Regression failed for surgeon", surgeon_name, "in Preop analysis\n"))
    return(list(
      surgeon = surgeon_name,
      n_total = n_total,
      n_preop = n_preop,
      n_change = 0,
      preop_r2 = NA,
      preop_p = NA,
      change_r2 = NA,
      change_p = NA
    ))
  }
  
  preop_r2 <- summary_preop$r.squared
  preop_p <- summary_preop$coefficients[2, 4]
  
  # Analysis 2: PI-LL vs Change in Knee Angle
  data_change <- df_surgeon %>%
    filter(!is.na(LAT6W_PI_LL) & !is.na(knee_angle_diff))
  
  n_change <- nrow(data_change)
  
  if (n_change < 2) {
    cat(paste("Warning: Insufficient data (n < 2) for surgeon", surgeon_name, "in Change analysis\n"))
    return(list(
      surgeon = surgeon_name,
      n_total = n_total,
      n_preop = n_preop,
      n_change = n_change,
      preop_r2 = preop_r2,
      preop_p = preop_p,
      change_r2 = NA,
      change_p = NA
    ))
  }
  
  # Regression: PI-LL vs Change in Knee Angle
  model_change <- lm(knee_angle_diff ~ LAT6W_PI_LL, data = data_change)
  summary_change <- summary(model_change)
  
  if (nrow(summary_change$coefficients) < 2) {
    cat(paste("Warning: Regression failed for surgeon", surgeon_name, "in Change analysis\n"))
    return(list(
      surgeon = surgeon_name,
      n_total = n_total,
      n_preop = n_preop,
      n_change = n_change,
      preop_r2 = preop_r2,
      preop_p = preop_p,
      change_r2 = NA,
      change_p = NA
    ))
  }
  
  change_r2 <- summary_change$r.squared
  change_p <- summary_change$coefficients[2, 4]
  
  cat(paste("  Preop Knee Flexion: R² =", round(preop_r2, 3), ", p =", formatC(preop_p, format = "e", digits = 2), "\n"))
  cat(paste("  Change in Knee Flexion: R² =", round(change_r2, 3), ", p =", formatC(change_p, format = "e", digits = 2), "\n\n"))
  
  # Create plots
  if (!dir.exists("results/by_surgeon")) {
    dir.create("results/by_surgeon", recursive = TRUE)
  }
  
  # Safe filename
  safe_surgeon_name <- gsub("[^A-Za-z0-9_]", "_", surgeon_name)
  
  # Plot 1: PI-LL vs Preop Knee Angle
  p1 <- ggplot(data_preop, aes(x = LAT6W_PI_LL, y = LATpre_LL_KneeAngle)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = "PI-LL Mismatch (6-week)",
      y = "Preoperative Knee Flexion",
      title = paste0("PI-LL vs Preoperative Knee Flexion\nSurgeon: ", surgeon_name),
      subtitle = paste0("n = ", n_preop, " (PI-LL between -10 and 10)")
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
      label = paste0("R² = ", round(preop_r2, 3),
                     "\np = ", formatC(preop_p, format = "e", digits = 2)),
      hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
    )
  
  filename1 <- paste0("analysis6_surgeon_", safe_surgeon_name, "_preop.png")
  filepath1 <- file.path("results/by_surgeon", filename1)
  ggsave(filepath1, plot = p1, width = 10, height = 8, dpi = 300)
  cat(paste("Saved Preop plot to", filepath1, "\n"))
  
  # Plot 2: PI-LL vs Change in Knee Angle
  p2 <- ggplot(data_change, aes(x = LAT6W_PI_LL, y = knee_angle_diff)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid") +
    labs(
      x = "PI-LL Mismatch (6-week)",
      y = "Change in Knee Flexion (6-week)",
      title = paste0("PI-LL vs Change in Knee Flexion\nSurgeon: ", surgeon_name),
      subtitle = paste0("n = ", n_change, " (PI-LL between -10 and 10)")
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
      label = paste0("R² = ", round(change_r2, 3),
                     "\np = ", formatC(change_p, format = "e", digits = 2)),
      hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
    )
  
  filename2 <- paste0("analysis6_surgeon_", safe_surgeon_name, "_change.png")
  filepath2 <- file.path("results/by_surgeon", filename2)
  ggsave(filepath2, plot = p2, width = 10, height = 8, dpi = 300)
  cat(paste("Saved Change plot to", filepath2, "\n\n"))
  
  return(list(
    surgeon = surgeon_name,
    n_total = n_total,
    n_preop = n_preop,
    n_change = n_change,
    preop_r2 = preop_r2,
    preop_p = preop_p,
    change_r2 = change_r2,
    change_p = change_p
  ))
}

# Run analysis for each surgeon
results_list <- list()

for (surgeon in surgeons) {
  result <- analysis_correctedLL_kf_by_surgeon(df_filtered, surgeon)
  results_list[[length(results_list) + 1]] <- result
}

# Convert results to data frame
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    Surgeon = x$surgeon,
    N_Total = x$n_total,
    N_Preop = x$n_preop,
    N_Change = x$n_change,
    Preop_R2 = x$preop_r2,
    Preop_P = x$preop_p,
    Change_R2 = x$change_r2,
    Change_P = x$change_p,
    stringsAsFactors = FALSE
  )
}))

# Load analysis3 cache to get the same surgeons
cache_file_analysis3 <- "results/analysis3_PCA_significance_tracking_cache.rds"
if (file.exists(cache_file_analysis3)) {
  cat("\n=== Loading analysis3 cache to get surgeon list ===\n")
  significance_tracking <- readRDS(cache_file_analysis3)
  
  # Convert to data frame and get surgeons with >= 20 cases
  tracking_df <- do.call(rbind, lapply(significance_tracking, function(x) {
    data.frame(
      surgeon = x$surgeon,
      n_cases = x$n_cases,
      stringsAsFactors = FALSE
    )
  }))
  
  surgeons_gt20 <- tracking_df %>%
    filter(n_cases >= 20) %>%
    filter(!is.na(n_cases))
  
  surgeons_to_include <- surgeons_gt20$surgeon
  cat(paste("Found", length(surgeons_to_include), "surgeons with >= 20 cases in analysis3\n"))
  
  # Filter to only include these surgeons
  results_df_filtered <- results_df %>%
    filter(Surgeon %in% surgeons_to_include)
} else {
  cat("\n=== Warning: analysis3 cache not found, using >= 20 cases filter ===\n")
  # Fallback: filter to surgeons with >= 20 cases
  results_df_filtered <- results_df %>%
    filter(N_Total >= 20)
}

cat("\n=== Results for Surgeons (same as analysis3) ===\n")
cat(paste("Total surgeons included:", nrow(results_df_filtered), "\n\n"))

# Print table (without N columns)
print_table <- results_df_filtered %>%
  select(Surgeon, Preop_R2, Preop_P, Change_R2, Change_P)
print(print_table, row.names = FALSE)

# Create formatted table for figure
if (nrow(results_df_filtered) > 0) {
  cat("\n=== Creating formatted table for figure ===\n")
  
  # Format p-values
  format_p_value <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("< 0.001")
    return(sprintf("%.3f", p))
  }
  
  # Format R²
  format_r2 <- function(r2) {
    if (is.na(r2)) return("NA")
    return(sprintf("%.3f", r2))
  }
  
  # Create formatted data frame
  table_df <- results_df_filtered %>%
    mutate(
      Preop_R2_Formatted = sapply(Preop_R2, format_r2),
      Preop_P_Formatted = sapply(Preop_P, format_p_value),
      Change_R2_Formatted = sapply(Change_R2, format_r2),
      Change_P_Formatted = sapply(Change_P, format_p_value)
    ) %>%
    select(Surgeon, 
           Preop_R2_Formatted, Preop_P_Formatted,
           Change_R2_Formatted, Change_P_Formatted)
  
  # Rename columns for display
  colnames(table_df) <- c("Surgeon",
                          "R² (Preop)", "p-value (Preop)",
                          "R² (Change)", "p-value (Change)")
  
  # Calculate image dimensions
  n_rows <- nrow(table_df)
  img_width <- 14
  img_height <- max(6, 2 + n_rows * 0.4)
  
  # Create table plot
  png("results/analysis6_by_surgeon_table.png", width = img_width, height = img_height, units = "in", res = 300)
  
  # Create table grob
  table_grob <- tableGrob(
    table_df,
    rows = NULL,
    theme = ttheme_default(
      base_size = 10,
      core = list(
        fg_params = list(hjust = 0, x = 0.05),
        bg_params = list(fill = c("white", "gray95"))
      ),
      colhead = list(
        fg_params = list(fontface = "bold", hjust = 0, x = 0.05),
        bg_params = list(fill = "gray80")
      )
    )
  )
  
  # Add title
  title_text <- paste0("PI-LL vs Knee Flexion Analysis by Surgeon\n",
                       "Same surgeons as analysis3 (PI-LL between -10 and 10)")
  title_grob <- textGrob(title_text, gp = gpar(fontsize = 14, fontface = "bold"))
  
  # Combine title and table
  final_grob <- arrangeGrob(
    title_grob,
    table_grob,
    ncol = 1,
    heights = c(0.1, 0.9),
    padding = unit(0.5, "in")
  )
  
  grid.draw(final_grob)
  dev.off()
  
  cat("Saved formatted table image to results/analysis6_by_surgeon_table.png\n")
}

cat("\n=== All analyses completed ===\n")

