# Utility functions for data loading and preprocessing

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

#' Load and preprocess COMBINE data with optional PJK exclusion
#'
#' @param file_path Path to the Excel database file
#' @param exclude_pjk Logical. If TRUE, excludes patients with PJK (AE_Spine_Radio_PJK > 0). Default is TRUE.
#' @return A data frame with the COMBINE data, optionally filtered to exclude PJK patients
#'
#' @export
load_combine_data <- function(file_path, exclude_pjk = TRUE) {
  # Load COMBINE sheet
  df <- readxl::read_excel(
    file_path,
    sheet = "COMBINE",
    .name_repair = ~make.unique(.x, sep = "__")
  )
  cat("Loaded CADS database\n")
  
  if (exclude_pjk) {
    # Load AE_flat sheet and exclude patients with PJK
    ae_flat <- readxl::read_excel(
      file_path,
      sheet = "AE_flat",
      .name_repair = ~make.unique(.x, sep = "__")
    )
    cat("Loaded AE_flat sheet\n")
    
    # Merge with COMBINE data to get PJK information
    # Find the demo_id column in COMBINE (could be demo_id or demo_id...1 after name repair)
    demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
    if (is.na(demo_id_col) || length(demo_id_col) == 0) {
      demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
      if (is.na(demo_id_col) || length(demo_id_col) == 0) {
        demo_id_col <- "demo_id"
      }
    }
    
    # Merge datasets
    df <- df %>%
      left_join(ae_flat %>% select(demo_id, AE_Spine_Radio_PJK), 
                by = setNames("demo_id", demo_id_col))
    
    # Exclude patients with PJK (AE_Spine_Radio_PJK > 0)
    n_before_pjk <- nrow(df)
    df <- df %>%
      filter(is.na(AE_Spine_Radio_PJK) | AE_Spine_Radio_PJK <= 0)
    n_after_pjk <- nrow(df)
    n_excluded_pjk <- n_before_pjk - n_after_pjk
    
    cat(paste("Excluded", n_excluded_pjk, "patients with PJK (AE_Spine_Radio_PJK > 0)\n"))
    cat(paste("Patients remaining after PJK exclusion:", n_after_pjk, "\n"))
  } else {
    cat("PJK exclusion is disabled\n")
  }
  
  return(df)
}

#' Filter data to exclude patients from surgeons with low case volumes
#'
#' @param df Data frame with the COMBINE data
#' @param min_surgeon_cases Integer. Minimum number of cases per surgeon to keep. Default is 6 (excludes surgeons with 5 or fewer).
#' @param surgeon_col Character. Name of the surgeon column. Default is "SxI_Gen_Surgeon".
#' @return A data frame with patients from low-volume surgeons excluded
#'
#' @export
filter_low_volume_surgeons <- function(df, min_surgeon_cases = 6, surgeon_col = "SxI_Gen_Surgeon") {
  # Count cases per surgeon
  surgeon_counts <- df %>%
    filter(!is.na(!!sym(surgeon_col))) %>%
    count(!!sym(surgeon_col), name = "case_count")
  
  # Identify surgeons with sufficient cases
  high_volume_surgeons <- surgeon_counts %>%
    filter(case_count >= min_surgeon_cases) %>%
    pull(!!sym(surgeon_col))
  
  n_before <- nrow(df)
  
  # Filter to keep only patients from high-volume surgeons
  df_filtered <- df %>%
    filter(is.na(!!sym(surgeon_col)) | !!sym(surgeon_col) %in% high_volume_surgeons)
  
  n_after <- nrow(df_filtered)
  n_excluded <- n_before - n_after
  
  cat(paste("Excluded", n_excluded, "patients from surgeons with", min_surgeon_cases - 1, "or fewer cases\n"))
  cat(paste("Patients remaining after surgeon volume filter:", n_after, "\n"))
  cat(paste("Number of surgeons retained:", length(high_volume_surgeons), "\n"))
  
  return(df_filtered)
}

