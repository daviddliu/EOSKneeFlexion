#!/usr/bin/env Rscript
# analysis23_planned_vs_achieved_PILL.R
#
# (1) Planned ΔLL (planned_DLL) vs actual change in lumbar lordosis (L1–S1) at 6 weeks.
# (2) Histogram: planned_DLL − achieved ΔL1–S1 (constant-PI planned lordosis change).
# (3) 6-week PI–LL (LAT6W_PI_LL) vs planned alignment goal (planned_PI_LL from PILL).
# Cohort: |preop PI–LL| > threshold; no PREOP_SVA_C2_C7_LT filter.
#
# Optional: analysis24_manually_verify_planned_PI_LL_gt20.csv with column "Exclude?" = Y
# removes those demo_ids from all analyses below (if the file exists).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

# |LATpre_PI_LL| must exceed this (°)
PREOP_ABS_PI_LL_MIN <- 10

# Drop implausible angles (same cap as analysis22)
MAX_DEG <- 100

# Rows with Exclude? = Y are dropped from df when this file is present (project root).
MANUAL_EXCLUDE_PLANNED_PILL_CSV <- "analysis24_manually_verify_planned_PI_LL_gt20.csv"

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)

# Actual lordosis change: 6-week L1–S1 minus preop (degrees)
df$change_lordosis <- df$LAT6W_L1_S1 - df$LATpre_L1_S1

demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
}
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- "demo_id"
}

if (file.exists(MANUAL_EXCLUDE_PLANNED_PILL_CSV)) {
  mx <- readr::read_csv(MANUAL_EXCLUDE_PLANNED_PILL_CSV, show_col_types = FALSE)
  ex_nm <- names(mx)[names(mx) == "Exclude?" |
    tolower(trimws(gsub("\\?", "", names(mx), fixed = TRUE))) == "exclude"]
  if (length(ex_nm) != 1L) {
    warning(
      "Manual exclusion CSV found but no single \"Exclude?\" column; skipping exclusions.",
      call. = FALSE
    )
  } else {
    id_nm <- names(mx)[tolower(names(mx)) == "demo_id"][1]
    if (is.na(id_nm)) {
      id_nm <- grep("^demo_id", names(mx), ignore.case = TRUE, value = TRUE)[1]
    }
    if (is.na(id_nm) || !nzchar(id_nm)) {
      id_nm <- names(mx)[1]
    }
    fl <- toupper(trimws(as.character(mx[[ex_nm]]))) == "Y"
    drop_ids <- unique(as.character(mx[[id_nm]][fl]))
    drop_ids <- drop_ids[!is.na(drop_ids) & nzchar(drop_ids)]
    n0 <- nrow(df)
    df <- df %>%
      dplyr::filter(!as.character(.data[[demo_id_col]]) %in% drop_ids)
    cat(sprintf(
      "Manual exclusions (Exclude? = Y): removed %d row(s) via %s\n",
      n0 - nrow(df),
      MANUAL_EXCLUDE_PLANNED_PILL_CSV
    ))
  }
}

d <- df %>%
  dplyr::filter(
    !is.na(.data$planned_DLL),
    !is.na(.data$LATpre_PI_LL),
    abs(.data$LATpre_PI_LL) > PREOP_ABS_PI_LL_MIN,
    !is.na(.data$change_lordosis)
  )

n_before_angle <- nrow(d)
d <- d %>%
  dplyr::filter(
    abs(.data$LATpre_PI_LL) <= MAX_DEG,
    abs(.data$planned_PI_LL) <= MAX_DEG,
    abs(.data$planned_DLL) <= MAX_DEG,
    .data$planned_DLL >= PLANNED_DLL_MIN_KEEP
  )

cat(sprintf(
  paste0(
    "Cohort: |preop PI–LL| > %d°, planned_DLL + 6W L1–S1 change (no SVA C2–C7 screen)\n",
    "  After mismatch + completeness: n = %d\n"
  ),
  PREOP_ABS_PI_LL_MIN,
  n_before_angle
))
if (n_before_angle - nrow(d) > 0) {
  cat(sprintf(
    paste0(
      "  Excluded %d row(s): any |baseline|, |planned|, |Planned ΔLL| > %d° ",
      "or Planned ΔLL < %.0f°\n"
    ),
    n_before_angle - nrow(d),
    MAX_DEG,
    PLANNED_DLL_MIN_KEEP
  ))
}
cat(sprintf("  Final n = %d\n", nrow(d)))

# Planned ΔL1–S1 (°) equals planned_DLL if PI is constant (same as plot x-axis).
# Within 1 SD: |planned_DLL − achieved ΔL1–S1| ≤ SD(planned_DLL − achieved ΔL1–S1).
diff_pa <- d$planned_DLL - d$change_lordosis
sd_diff <- stats::sd(diff_pa)
if (is.na(sd_diff) || sd_diff <= 0) {
  cat("\nPlanned vs achieved Δlordosis (1 SD rule): SD(planned − achieved) not defined; skip.\n")
} else {
  within_1sd <- abs(diff_pa) <= sd_diff
  n_w <- sum(within_1sd)
  n_tot <- nrow(d)
  cat(sprintf(
    paste0(
      "\n|planned_DLL \u2212 achieved \u0394L1–S1| \u2264 SD(planned \u2212 achieved): ",
      "%d / %d (%.1f%%) | SD(planned \u2212 achieved) = %.2f\u00b0\n"
    ),
    n_w,
    n_tot,
    100 * n_w / n_tot,
    sd_diff
  ))
}

d_diff <- tibble::tibble(planned_minus_achieved = diff_pa)
cat(sprintf(
  paste0(
    "\nPlanned \u2212 achieved \u0394L1–S1 (°): mean = %.2f, SD = %.2f, median = %.2f\n"
  ),
  mean(diff_pa),
  stats::sd(diff_pa),
  stats::median(diff_pa)
))

# --- Plot 2: planned PI–LL goal vs 6-week PI–LL (separate completeness / angle screen) ---
d_pill <- df %>%
  dplyr::filter(
    !is.na(.data$planned_PI_LL),
    !is.na(.data$LAT6W_PI_LL),
    !is.na(.data$LATpre_PI_LL),
    abs(.data$LATpre_PI_LL) > PREOP_ABS_PI_LL_MIN
  )
n_pill_before <- nrow(d_pill)
d_pill <- d_pill %>%
  dplyr::filter(
    abs(.data$LATpre_PI_LL) <= MAX_DEG,
    abs(.data$planned_PI_LL) <= MAX_DEG,
    abs(.data$LAT6W_PI_LL) <= MAX_DEG
  )
cat(sprintf(
  "\nPlanned PILL vs 6W PI–LL: complete + |angles| ≤ %d°: n = %d",
  MAX_DEG,
  nrow(d_pill)
))
if (n_pill_before - nrow(d_pill) > 0) {
  cat(sprintf(" (excluded %d row(s) by angle cap)\n", n_pill_before - nrow(d_pill)))
} else {
  cat("\n")
}

model <- lm(change_lordosis ~ planned_DLL, data = d)
sm <- summary(model)
r_squared <- sm$r.squared
p_value <- sm$coefficients[2L, 4L]
pearson_r <- cor(d$planned_DLL, d$change_lordosis, use = "complete.obs")

r_text <- paste0("r = ", round(pearson_r, 3))
r2_text <- paste0("R² = ", round(r_squared, 3))
p_text <- paste0("p = ", formatC(p_value, format = "e", digits = 2))
annotation_text <- paste0(r_text, "\n", r2_text, "\n", p_text)

p <- ggplot(d, aes(x = planned_DLL, y = change_lordosis)) +
  geom_point(alpha = 0.55) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linewidth = 0.6) +
  labs(
    x = "Planned ΔLL (°) = preop PI–LL − planned PI–LL",
    y = "Actual ΔL1–S1 lordosis (°) = LAT6W − preop",
    title = "Planned ΔLL vs actual change in lordosis"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = annotation_text,
    hjust = 1.05, vjust = 1.2,
    size = 3.8,
    fontface = "bold"
  )

model_pill <- lm(planned_PI_LL ~ LAT6W_PI_LL, data = d_pill)
sm_pill <- summary(model_pill)
r2_pill <- sm_pill$r.squared
p_pill <- sm_pill$coefficients[2L, 4L]
r_pill <- cor(d_pill$LAT6W_PI_LL, d_pill$planned_PI_LL, use = "complete.obs")
ann_pill <- paste0(
  "r = ", round(r_pill, 3), "\n",
  "R² = ", round(r2_pill, 3), "\n",
  "p = ", formatC(p_pill, format = "e", digits = 2)
)

p_pill <- ggplot(d_pill, aes(x = LAT6W_PI_LL, y = planned_PI_LL)) +
  geom_point(alpha = 0.55) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linewidth = 0.6) +
  labs(
    x = "6-week PI–LL (LAT6W_PI_LL, °)",
    y = "Planned PI–LL (parsed PILL goal, °)",
    title = "6-week PI–LL vs planned alignment goal"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = ann_pill,
    hjust = 1.05, vjust = 1.2,
    size = 3.8,
    fontface = "bold"
  )

p_hist_diff <- ggplot(d_diff, aes(x = planned_minus_achieved)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 40, boundary = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray25") +
  geom_vline(
    xintercept = mean(d_diff$planned_minus_achieved),
    color = "firebrick",
    linewidth = 0.75
  ) +
  labs(
    x = "Planned \u2212 achieved \u0394L1–S1 (°) (planned_DLL \u2212 [LAT6W\u2212pre L1\u2013S1])",
    y = "Patients",
    title = "Distribution of planned minus achieved lordosis change",
    subtitle = sprintf(
      "n = %d | mean = %.2f\u00b0 | SD = %.2f\u00b0 | red = mean",
      nrow(d_diff),
      mean(d_diff$planned_minus_achieved),
      stats::sd(d_diff$planned_minus_achieved)
    )
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

pr_dir <- ensure_planned_results_dir()
out_png <- file.path(pr_dir, "analysis23_planned_dLL_vs_actual_lordosis_change.png")
ggsave(out_png, plot = p, width = 9, height = 6.5, dpi = 300)
cat(sprintf("Saved %s\n", out_png))

out_png_hist <- file.path(pr_dir, "analysis23_planned_minus_achieved_delta_LL_histogram.png")
ggsave(out_png_hist, plot = p_hist_diff, width = 9, height = 5.5, dpi = 300)
cat(sprintf("Saved %s\n", out_png_hist))

out_png_pill <- file.path(pr_dir, "analysis23_planned_PILL_vs_6W_PI_LL.png")
ggsave(out_png_pill, plot = p_pill, width = 9, height = 6.5, dpi = 300)
cat(sprintf("Saved %s\n", out_png_pill))
