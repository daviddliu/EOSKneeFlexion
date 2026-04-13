#!/usr/bin/env Rscript
# Table: Pearson r and p-values for |KF| vs ODI / SRS-22 total / SRS-22 subdomains
# at all HRQL timepoints (baseline, 6 weeks, 1 year, 2 years).
# Multiple comparisons: Benjamini–Hochberg FDR within each stratum sheet.

suppressPackageStartupMessages({
  library(tidyverse)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# COMBINE already includes HRQL/PROM fields; joining the HRQL sheet again duplicates
# column names (e.g. BL_HRQL_SRS_TOTAL.x / .y) and breaks grep-based lookups.

time_map <- tibble::tribble(
  ~hrql_prefix, ~time_label, ~kf_col,
  "BL",  "Preoperative", "LATpre_LL_KneeAngle",
  "W6",  "6 weeks",      "LAT6W_LL_KneeAngle",
  "Y1",  "1 year",       "LAT1Y_LL_KneeAngle",
  "Y2",  "2 years",      "LAT2Y_LL_KneeAngle"
)

outcome_defs <- tribble(
  ~outcome_key, ~suffix_pattern, ~outcome_label,
  "ODI",         "ODI",          "ODI (total)",
  "SRS_TOTAL",   "SRS_TOTAL",    "SRS-22 total",
  "ACTIVITY",    "SRS_ACTIVITY", "SRS-22 Activity",
  "PAIN",        "SRS_PAIN",     "SRS-22 Pain",
  "APPEARANCE",  "SRS_APPEARANCE", "SRS-22 Appearance",
  "MENTAL",      "SRS_MENTAL",   "SRS-22 Mental health",
  "SATIS",       "SRS_SATIS",    "SRS-22 Satisfaction with management"
)

find_hrql_col <- function(dat, prefix, suffix_pattern) {
  pat <- paste0("^", prefix, "_HRQL_", suffix_pattern)
  hit <- grep(pat, names(dat), value = TRUE, ignore.case = TRUE)
  # Keep only exact PROM column (e.g. BL_HRQL_ODI) or readxl duplicate (…__1), not BL_HRQL_ODI_Date / ODI_Pain
  hit <- hit[grepl(paste0("^", prefix, "_HRQL_", suffix_pattern, "(__[0-9]+)?$"), hit, ignore.case = TRUE)]
  exact <- paste0(prefix, "_HRQL_", suffix_pattern)
  hit <- hit[order(hit != exact)]
  if (length(hit)) hit[1] else NA_character_
}

correlation_row <- function(dat, time_label, hrql_prefix, kf_col, outcome_label, ycol) {
  if (is.na(ycol) || !(ycol %in% names(dat))) {
    return(NULL)
  }
  if (!(kf_col %in% names(dat))) {
    return(NULL)
  }
  x <- abs(as.numeric(dat[[kf_col]]))
  y <- as.numeric(dat[[ycol]])
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 3) {
    return(NULL)
  }
  fit <- lm(y ~ x)
  sm <- summary(fit)
  tibble(
    Timepoint = time_label,
    `HRQL prefix` = hrql_prefix,
    Outcome = outcome_label,
    `HRQL column` = ycol,
    `KF column` = kf_col,
    n = n,
    r = unname(cor(x, y)),
    R2 = sm$r.squared,
    p = sm$coefficients[2, 4],
    slope = sm$coefficients[2, 1]
  )
}

run_all_for_dat <- function(dat) {
  rows <- list()
  for (ti in seq_len(nrow(time_map))) {
    pref <- time_map$hrql_prefix[ti]
    tlab <- time_map$time_label[ti]
    kf <- time_map$kf_col[ti]
    for (j in seq_len(nrow(outcome_defs))) {
      ycol <- find_hrql_col(dat, pref, outcome_defs$suffix_pattern[j])
      row <- correlation_row(dat, tlab, pref, kf, outcome_defs$outcome_label[j], ycol)
      if (!is.null(row)) rows[[length(rows) + 1]] <- row
    }
  }
  bind_rows(rows)
}

add_mc <- function(tbl) {
  if (nrow(tbl) == 0) {
    return(mutate(tbl, `p (BH-FDR)` = numeric(0)))
  }
  mutate(tbl, `p (BH-FDR)` = p.adjust(p, method = "fdr"))
}

out_all <- run_all_for_dat(df) %>%
  mutate(`PT stratum` = "All") %>%
  add_mc()

df_low <- df %>% filter(!is.na(LATpre_S1PT), LATpre_S1PT < 20)
df_high <- df %>% filter(!is.na(LATpre_S1PT), LATpre_S1PT >= 20)

out_low <- run_all_for_dat(df_low) %>%
  mutate(`PT stratum` = "Preop PT < 20") %>%
  add_mc()
out_high <- run_all_for_dat(df_high) %>%
  mutate(`PT stratum` = "Preop PT >= 20") %>%
  add_mc()

out_combined <- bind_rows(out_all, out_low, out_high) %>%
  select(
    `PT stratum`, Timepoint, Outcome, n, r, R2, p, `p (BH-FDR)`,
    `HRQL prefix`, `HRQL column`, `KF column`, slope
  )

notes <- tibble::tibble(Note = c(
  "Exposure: absolute knee flexion (|KF|) at the same imaging timepoint as the PROM (matched LAT*_LL_KneeAngle).",
  "Association: linear regression PROM ~ |KF|; r = Pearson correlation; p = slope t-test p-value (two-sided).",
  "p (BH-FDR): Benjamini–Hochberg FDR within each sheet (separate adjustment for All, Preop PT < 20, Preop PT >= 20).",
  "PROMs are read from the COMBINE sheet (HRQL fields embedded there for this database).",
  "PJK exclusion and database path follow other analyses in this repository.",
  "Regenerate: Rscript analysis20_KF_ODI_SRS_correlations_excel.R"
))

if (!dir.exists("results")) dir.create("results")

excel_path <- "results/KF_ODI_SRS22_correlations_by_timepoint.xlsx"

if (!requireNamespace("writexl", quietly = TRUE)) {
  stop("Install writexl: install.packages(\"writexl\")")
}

writexl::write_xlsx(
  list(
    All = select(out_all, -`PT stratum`),
    `Preop PT lt 20` = select(out_low, -`PT stratum`),
    `Preop PT ge 20` = select(out_high, -`PT stratum`),
    Combined = out_combined,
    Notes = notes
  ),
  excel_path
)

cat("\nWrote ", excel_path, "\n", sep = "")
cat("Rows (All / low PT / high PT):", nrow(out_all), nrow(out_low), nrow(out_high), "\n")
