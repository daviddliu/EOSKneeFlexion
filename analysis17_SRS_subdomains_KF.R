#!/usr/bin/env Rscript
# SRS-22 subdomain scores vs absolute knee flexion (matched timepoint).
# Subdomains: Activity, Pain, Appearance, Mental health, Satisfaction with management.
# Bonferroni: within each PT stratum (and within "All"), 5 subdomains × 3 timepoints = 15 tests.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

hrql_df <- read_excel(db_path, sheet = "HRQL", .name_repair = ~ make.unique(.x, sep = "__"))

demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || length(demo_id_col) == 0) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
  if (is.na(demo_id_col) || length(demo_id_col) == 0) demo_id_col <- "demo_id"
}
hrql_demo_id_col <- grep("^BL_demo_id$", names(hrql_df), value = TRUE)[1]
if (is.na(hrql_demo_id_col) || length(hrql_demo_id_col) == 0) {
  hrql_demo_id_col <- grep("^BL_demo_id", names(hrql_df), value = TRUE)[1]
  if (is.na(hrql_demo_id_col) || length(hrql_demo_id_col) == 0) hrql_demo_id_col <- "BL_demo_id"
}

df <- df %>%
  left_join(hrql_df, by = setNames(hrql_demo_id_col, demo_id_col))

subdomains <- c("ACTIVITY", "PAIN", "APPEARANCE", "MENTAL", "SATIS")
subdomain_labels <- c(
  ACTIVITY = "Activity",
  PAIN = "Pain",
  APPEARANCE = "Appearance",
  MENTAL = "Mental health",
  SATIS = "Satisfaction with management"
)

time_map <- tibble::tribble(
  ~hrql_prefix, ~time_label, ~kf_col,
  "BL",  "Preoperative", "LATpre_LL_KneeAngle",
  "W6",  "6 Weeks",      "LAT6W_LL_KneeAngle",
  "Y2",  "2 Years",      "LAT2Y_LL_KneeAngle"
)

find_srs_col <- function(dat, prefix, sub) {
  pat <- paste0("^", prefix, "_HRQL_SRS_", sub)
  hit <- grep(pat, names(dat), value = TRUE, ignore.case = TRUE)
  exact <- paste0(prefix, "_HRQL_SRS_", sub)
  hit <- hit[order(hit != exact)]
  if (length(hit)) hit[1] else NA_character_
}

#' One stratum: all subdomain × timepoint correlations |KF| ~ SRS subdomain
run_subdomain_correlations <- function(dat, label = "") {
  rows <- list()
  for (ti in seq_len(nrow(time_map))) {
    pref <- time_map$hrql_prefix[ti]
    tlab <- time_map$time_label[ti]
    kf <- time_map$kf_col[ti]
    if (!kf %in% names(dat)) next
    for (sub in subdomains) {
      ycol <- find_srs_col(dat, pref, sub)
      if (is.na(ycol) || !ycol %in% names(dat)) next
      x <- abs(as.numeric(dat[[kf]]))
      y <- as.numeric(dat[[ycol]])
      ok <- !is.na(x) & !is.na(y)
      x <- x[ok]
      y <- y[ok]
      n <- length(x)
      if (n < 3) next
      fit <- lm(y ~ x)
      sm <- summary(fit)
      rows[[length(rows) + 1]] <- tibble::tibble(
        Timepoint = tlab,
        Subdomain = subdomain_labels[[sub]],
        SRS_column = ycol,
        KF_column = kf,
        n = n,
        r = unname(cor(x, y)),
        R2 = sm$r.squared,
        p = sm$coefficients[2, 4],
        slope = sm$coefficients[2, 1]
      )
    }
  }
  bind_rows(rows)
}

add_bonferroni <- function(tbl) {
  k <- nrow(tbl)
  if (k == 0) return(tbl %>% mutate(p_bonferroni = numeric(0)))
  tbl %>% mutate(p_bonferroni = pmin(p * k, 1))
}

# --- Full cohort ---
out_all <- run_subdomain_correlations(df, "All") %>% mutate(PT_group = "All")
out_all <- add_bonferroni(out_all)

# --- PT strata (preoperative S1PT) ---
df_low <- df %>% filter(!is.na(LATpre_S1PT), LATpre_S1PT < 20)
df_high <- df %>% filter(!is.na(LATpre_S1PT), LATpre_S1PT >= 20)

out_low <- run_subdomain_correlations(df_low) %>% mutate(PT_group = "Preop PT < 20")
out_high <- run_subdomain_correlations(df_high) %>% mutate(PT_group = "Preop PT >= 20")
out_low <- add_bonferroni(out_low)
out_high <- add_bonferroni(out_high)

out_pt <- bind_rows(out_low, out_high)

cat("\n=== SRS-22 subdomains vs |KF| — FULL COHORT ===\n")
cat("Bonferroni: ", nrow(out_all), " tests\n\n")
print(as.data.frame(out_all), row.names = FALSE, digits = 3)

cat("\n=== SRS-22 subdomains vs |KF| — BY PREOP PT ===\n")
cat("Bonferroni applied separately within each PT group (", nrow(out_low), " and ", nrow(out_high), " tests)\n\n", sep = "")
print(as.data.frame(out_pt), row.names = FALSE, digits = 3)

if (!dir.exists("planned_results")) dir.create("planned_results")
readr::write_csv(out_all, "planned_results/analysis17_SRS_subdomains_KF_summary.csv")
readr::write_csv(out_pt, "planned_results/analysis17_SRS_subdomains_KF_by_PT_summary.csv")
cat("\nWrote planned_results/analysis17_SRS_subdomains_KF_summary.csv\n")
cat("Wrote planned_results/analysis17_SRS_subdomains_KF_by_PT_summary.csv\n")

for (g in unique(out_pt$PT_group)) {
  sig <- out_pt %>% filter(PT_group == g, p_bonferroni < 0.05)
  cat("\n", g, ": Bonferroni p < 0.05 (n = ", nrow(sig), ")\n", sep = "")
  if (nrow(sig) > 0) print(as.data.frame(sig), row.names = FALSE)
}
