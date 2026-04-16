#!/usr/bin/env Rscript
#
# Rows with parsed planned PI–LL goal > threshold: print and write CSV of preop PI, L1–S1 LL,
# and raw PILL text for manual verification of extreme alignment goals.

suppressPackageStartupMessages({
  library(tidyverse)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

PLANNED_PI_LL_GT <- 20

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)

demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
}
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- "demo_id"
}

pi_col <- NULL
if ("PI" %in% names(df)) {
  pi_col <- "PI"
} else if ("LATpre_PI" %in% names(df)) {
  pi_col <- "LATpre_PI"
} else if ("LATpre_S1PI" %in% names(df)) {
  pi_col <- "LATpre_S1PI"
  cat("Using LATpre_S1PI as preop PI (PI column not found).\n")
}

if (is.null(pi_col)) {
  stop("No pelvic incidence column found (expected PI, LATpre_PI, or LATpre_S1PI).")
}

if (!"LATpre_L1_S1" %in% names(df)) {
  stop("COMBINE data must contain LATpre_L1_S1 for preoperative L1–S1 lordosis.")
}

raw_col <- if ("demo_AlignGoal_PILL_raw" %in% names(df)) {
  "demo_AlignGoal_PILL_raw"
} else {
  NULL
}

out <- df %>%
  dplyr::filter(!is.na(.data$planned_PI_LL), .data$planned_PI_LL > PLANNED_PI_LL_GT) %>%
  dplyr::mutate(
    demo_id = .data[[demo_id_col]],
    preop_PI = .data[[pi_col]],
    preop_LL_L1_S1 = .data$LATpre_L1_S1
  ) %>%
  dplyr::arrange(dplyr::desc(.data$planned_PI_LL))

show <- out %>%
  dplyr::transmute(
    demo_id = .data$demo_id,
    planned_PI_LL = .data$planned_PI_LL,
    preop_PI = .data$preop_PI,
    preop_LL_L1_S1 = .data$preop_LL_L1_S1,
    LATpre_PI_LL = .data$LATpre_PI_LL
  )
if (!is.null(raw_col)) {
  show$demo_AlignGoal_PILL_raw <- stringr::str_trunc(as.character(out[[raw_col]]), 70)
}

cat(sprintf(
  "\nPlanned PI–LL > %.0f°: n = %d | preop PI from column \"%s\" | preop LL = LATpre_L1_S1\n",
  PLANNED_PI_LL_GT,
  nrow(out),
  pi_col
))
cat("(demo_AlignGoal_PILL_raw truncated to 70 chars for console.)\n\n")

if (nrow(show) == 0) {
  cat("(no rows)\n")
} else {
  print(show, n = Inf, width = Inf)
}

pr_dir <- ensure_planned_results_dir()
csv_path <- file.path(pr_dir, "analysis24_verify_planned_PI_LL_gt20.csv")
csv_tbl <- out %>%
  dplyr::mutate(preop_PI_source = pi_col) %>%
  dplyr::select(
    "demo_id",
    "planned_PI_LL",
    "preop_PI",
    "preop_PI_source",
    "preop_LL_L1_S1",
    "LATpre_PI_LL",
    dplyr::any_of("demo_AlignGoal_PILL_raw")
  )
readr::write_csv(csv_tbl, csv_path, na = "")
cat(sprintf("\nWrote %s (%d rows, full raw PILL text)\n", csv_path, nrow(csv_tbl)))
