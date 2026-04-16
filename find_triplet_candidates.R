#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(parallel)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration: Toggle PJK exclusion on/off
EXCLUDE_PJK <- TRUE

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Get patient ID column (demo_id)
demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || length(demo_id_col) == 0) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
  if (is.na(demo_id_col) || length(demo_id_col) == 0) {
    demo_id_col <- "demo_id"
  }
}

# Filter for patients with satisfactory 6-week alignment (-10 < PI_LL < 10)
# and complete knee flexion data
df_clean <- df %>%
  filter(!is.na(LAT6W_PI_LL) &
         LAT6W_PI_LL > -10 & 
         LAT6W_PI_LL < 10 &
         !is.na(LATpre_LL_KneeAngle) &
         !is.na(LAT6W_LL_KneeAngle))

# Calculate change in knee flexion (preop to 6-week)
df_clean$change_kf <- df_clean$LAT6W_LL_KneeAngle - df_clean$LATpre_LL_KneeAngle

cat(paste("Found", nrow(df_clean), "patients with satisfactory 6-week alignment (-10 < PI_LL < 10)\n"))
cat(paste("and complete knee flexion data (preop and 6-week)\n"))
cat(paste("Patient ID column:", demo_id_col, "\n\n"))

# Function to calculate triplet score (spread of knee flexion changes)
calculate_triplet_score <- function(ids, data, demo_id_col) {
  triplet_data <- data[data[[demo_id_col]] %in% ids, ]
  
  if (nrow(triplet_data) != 3) return(NULL)
  
  # Calculate spread: max change - min change (larger is better)
  kf_changes <- triplet_data$change_kf
  kf_spread <- max(kf_changes) - min(kf_changes)
  
  return(list(
    ids = ids,
    kf_spread = kf_spread,
    kf_changes = sort(kf_changes),
    data = triplet_data
  ))
}

# Function to process a batch of combinations (for parallel processing)
process_combination_batch <- function(batch_indices, all_ids, df_clean, demo_id_col, combs_matrix) {
  batch_results <- list()
  for (idx in batch_indices) {
    ids <- all_ids[combs_matrix[, idx]]
    result <- calculate_triplet_score(ids, df_clean, demo_id_col)
    if (!is.null(result)) {
      batch_results[[length(batch_results) + 1]] <- result
    }
  }
  return(batch_results)
}

# Find all triplets and calculate their scores
cat("Finding triplets with maximum spread in knee flexion change...\n")
all_ids <- df_clean[[demo_id_col]]
n_patients <- length(all_ids)

cat(paste("Total eligible patients:", n_patients, "\n"))
n_combinations <- choose(n_patients, 3)
cat(paste("Total combinations to check:", n_combinations, "\n\n"))

if (n_combinations > 1000000) {
  cat("WARNING: Large number of combinations. This may take a while.\n\n")
}

cat("Generating combinations and calculating triplet scores...\n")
cat("Using parallel processing with 4 threads...\n")

# Use simplify = TRUE for better memory efficiency
# Generate combinations as a matrix (more memory efficient)
combs_matrix <- combn(n_patients, 3, simplify = TRUE)
n_combs <- ncol(combs_matrix)

cat(sprintf("Processing %d combinations with 4 parallel threads...\n", n_combs))

# Set up parallel processing
n_cores <- 4
cl <- makeCluster(n_cores, type = "PSOCK")

# Load required libraries on workers
clusterEvalQ(cl, {
  library(dplyr)
})

# Export necessary variables and functions to workers
clusterExport(cl, c("all_ids", "df_clean", "demo_id_col", "combs_matrix", 
                    "calculate_triplet_score", "process_combination_batch"))

# Split work into batches for parallel processing
batch_size <- max(1, floor(n_combs / (n_cores * 10)))  # Create enough batches for good load balancing
batches <- split(seq_len(n_combs), ceiling(seq_len(n_combs) / batch_size))

cat(sprintf("Processing in %d batches across %d cores...\n", length(batches), n_cores))

# Process batches in parallel
start_time <- Sys.time()
all_triplets_list <- parLapply(cl, batches, function(batch_indices) {
  process_combination_batch(batch_indices, all_ids, df_clean, demo_id_col, combs_matrix)
})

# Combine results from all batches
all_triplets <- do.call(c, all_triplets_list)

# Clean up cluster
stopCluster(cl)

elapsed_time <- Sys.time() - start_time
cat(sprintf("Completed in %.2f seconds\n", as.numeric(elapsed_time, units = "secs")))
cat("\n")  # New line after progress

cat(paste("\nFound", length(all_triplets), "triplets\n\n"))

if (length(all_triplets) == 0) {
  cat("No triplets found.\n")
  quit(status = 0)
}

# Sort by KF change spread (largest first - most different changes)
all_triplets <- all_triplets[order(sapply(all_triplets, function(x) -x$kf_spread))]

# Output results
cat("=== TRIPLET CANDIDATES (sorted by knee flexion change spread, largest first) ===\n\n")

# Limit output to top 20 triplets to avoid overwhelming output
n_output <- min(20, length(all_triplets))
cat(paste("Showing top", n_output, "triplets (out of", length(all_triplets), "total):\n\n"))

for (i in seq_len(n_output)) {
  triplet <- all_triplets[[i]]
  cat(sprintf("--- Triplet %d (KF change spread: %.2f degrees) ---\n", i, triplet$kf_spread))
  cat("\nPatient details:\n")
  
  # Sort by KF change for display
  triplet_data_sorted <- triplet$data[order(triplet$data$change_kf), ]
  
  for (j in 1:3) {
    patient <- triplet_data_sorted[j, ]
    cat(sprintf("  Patient %d (ID: %s):\n", j, patient[[demo_id_col]]))
    cat(sprintf("    KF preop: %.2f degrees\n", patient$LATpre_LL_KneeAngle))
    cat(sprintf("    KF postop (6W): %.2f degrees\n", patient$LAT6W_LL_KneeAngle))
    cat(sprintf("    KF change: %.2f degrees\n", patient$change_kf))
    cat(sprintf("    PI-LL (6W): %.2f degrees\n", patient$LAT6W_PI_LL))
    cat("\n")
  }
  cat("\n")
}

# Also save to CSV for easier analysis
output_df <- data.frame(
  triplet_num = integer(),
  patient_num = integer(),
  patient_id = character(),
  kf_preop = numeric(),
  kf_postop_6w = numeric(),
  kf_change = numeric(),
  pi_ll_6w = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(min(50, length(all_triplets)))) {  # Save top 50 triplets
  triplet <- all_triplets[[i]]
  triplet_data_sorted <- triplet$data[order(triplet$data$change_kf), ]
  
  for (j in 1:3) {
    patient <- triplet_data_sorted[j, ]
    output_df <- rbind(output_df, data.frame(
      triplet_num = i,
      patient_num = j,
      patient_id = as.character(patient[[demo_id_col]]),
      kf_preop = patient$LATpre_LL_KneeAngle,
      kf_postop_6w = patient$LAT6W_LL_KneeAngle,
      kf_change = patient$change_kf,
      pi_ll_6w = patient$LAT6W_PI_LL,
      stringsAsFactors = FALSE
    ))
  }
}

# Add summary columns
summary_df <- output_df %>%
  group_by(triplet_num) %>%
  summarise(
    kf_change_spread = max(kf_change) - min(kf_change),
    kf_change_min = min(kf_change),
    kf_change_max = max(kf_change),
    kf_change_mean = mean(kf_change),
    .groups = 'drop'
  )

# Merge with main output
output_df <- output_df %>%
  left_join(summary_df, by = "triplet_num")

# Reorder columns for better readability
output_df <- output_df %>%
  select(triplet_num, patient_num, patient_id, 
         kf_preop, kf_postop_6w, kf_change, pi_ll_6w,
         kf_change_spread, kf_change_min, kf_change_max, kf_change_mean)

# Save to CSV
if (!dir.exists("planned_results")) {
  dir.create("planned_results")
}
output_file <- "planned_results/triplet_candidates.csv"
write.csv(output_df, output_file, row.names = FALSE)
cat(sprintf("Saved detailed results to %s\n", output_file))
cat(sprintf("Total triplets found: %d\n", length(all_triplets)))
