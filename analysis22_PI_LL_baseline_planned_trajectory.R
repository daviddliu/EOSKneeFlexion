#!/usr/bin/env Rscript
#
# Spaghetti + mean trajectory: preoperative PI–LL mismatch (LATpre_PI_LL) → planned
# PI–LL from Patients_Attributes$demo_AlignGoal_PILL (parsed + recode), plus Planned ΔLL histogram.
# Planned ΔLL (planned_DLL) = LATpre_PI_LL − planned_PI_LL (same convention as other PILL analyses).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)
df <- restrict_planned_dll_analysis_cohort(df)
cat(sprintf(
  "\nTrajectory cohort: |preop PI–LL| > %d°, preop SVA C2–C7 < %d (non-missing): n = %d rows\n",
  PREOP_ABS_PI_LL_GT,
  PREOP_SVA_C2_C7_LT,
  nrow(df)
))

demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
}
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- "demo_id"
}

if (!"LATpre_PI_LL" %in% names(df)) {
  stop("COMBINE must contain LATpre_PI_LL")
}

# Drop implausible angles (same cap as analysis19) so plots are not dominated by typos
MAX_DEG <- 100

traj_wide <- df %>%
  filter(
    !is.na(LATpre_PI_LL),
    !is.na(planned_PI_LL),
    !is.na(planned_DLL)
  ) %>%
  mutate(demo_id = .data[[demo_id_col]])
n_before_screen <- nrow(traj_wide)
traj <- traj_wide %>%
  filter(
    abs(LATpre_PI_LL) <= MAX_DEG,
    abs(planned_PI_LL) <= MAX_DEG,
    abs(planned_DLL) <= MAX_DEG,
    planned_DLL >= PLANNED_DLL_MIN_KEEP
  )
if (n_before_screen - nrow(traj) > 0) {
  cat(sprintf(
    paste0(
      "Excluded %d patient-row(s) with any |baseline|, |planned|, |Planned \u0394LL| > %d\u00b0 ",
      "or Planned \u0394LL < %.0f\u00b0\n"
    ),
    n_before_screen - nrow(traj),
    MAX_DEG,
    PLANNED_DLL_MIN_KEEP
  ))
}

n_pat <- n_distinct(traj$demo_id)
n_neg <- sum(traj$planned_DLL < 0)
n_zero <- sum(traj$planned_DLL == 0)
n_pos <- sum(traj$planned_DLL > 0)

cat("\n=== PI–LL: baseline → planned (parsed goal, recoded) ===\n")
cat(sprintf(
  "After |baseline|, |planned|, |Planned \u0394LL| \u2264 %d\u00b0 and Planned \u0394LL \u2265 %.0f\u00b0\n",
  MAX_DEG,
  PLANNED_DLL_MIN_KEEP
))
cat(sprintf("Patients with both baseline and planned PI–LL: n = %d\n", n_pat))
cat(sprintf(
  "Baseline PI–LL:  mean = %.2f\u00b0, SD = %.2f\u00b0\n",
  mean(traj$LATpre_PI_LL),
  stats::sd(traj$LATpre_PI_LL)
))
cat(sprintf(
  "Planned PI–LL:   mean = %.2f\u00b0, SD = %.2f\u00b0\n",
  mean(traj$planned_PI_LL),
  stats::sd(traj$planned_PI_LL)
))
cat(sprintf(
  "Planned \u0394LL (preop \u2212 planned): mean = %.2f\u00b0, SD = %.2f\u00b0\n",
  mean(traj$planned_DLL),
  stats::sd(traj$planned_DLL)
))
cat(sprintf(
  "Planned \u0394LL < 0: n = %d (%.1f%%); = 0: n = %d; > 0: n = %d (%.1f%%)\n",
  n_neg,
  if (n_pat > 0) 100 * n_neg / n_pat else NA_real_,
  n_zero,
  n_pos,
  if (n_pat > 0) 100 * n_pos / n_pat else NA_real_
))

phase_levels <- c("Baseline PI–LL", "Planned PI–LL")
long <- traj %>%
  transmute(
    demo_id,
    baseline = LATpre_PI_LL,
    planned = planned_PI_LL
  ) %>%
  pivot_longer(-demo_id, names_to = "phase_key", values_to = "PI_LL") %>%
  mutate(
    phase = factor(
      ifelse(phase_key == "baseline", phase_levels[1], phase_levels[2]),
      levels = phase_levels
    )
  )

summ <- long %>%
  group_by(phase) %>%
  summarize(mean_PI_LL = mean(PI_LL), .groups = "drop")

p_traj <- ggplot(long, aes(x = phase, y = PI_LL, group = demo_id)) +
  geom_line(alpha = 0.12, color = "gray40") +
  geom_line(
    data = summ,
    aes(x = phase, y = mean_PI_LL, group = 1),
    inherit.aes = FALSE,
    color = "firebrick",
    linewidth = 1
  ) +
  geom_point(
    data = summ,
    aes(x = phase, y = mean_PI_LL, group = 1),
    inherit.aes = FALSE,
    color = "firebrick",
    size = 2.8
  ) +
  labs(
    x = NULL,
    y = "PI–LL mismatch (\u00b0)",
    title = "PI–LL trajectory: preoperative → planned",
    subtitle = sprintf(
      "|values| \u2264 %d\u00b0 | n = %d | Mean Planned \u0394LL = %.2f\u00b0 | %.1f%% with Planned \u0394LL < 0",
      MAX_DEG,
      n_pat,
      mean(traj$planned_DLL),
      if (n_pat > 0) 100 * n_neg / n_pat else NA_real_
    )
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_hist <- ggplot(traj, aes(x = planned_DLL)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 40, boundary = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(
    xintercept = mean(traj$planned_DLL),
    color = "firebrick",
    linewidth = 0.8
  ) +
  labs(
    x = "Planned \u0394LL (\u00b0) = preop \u2212 planned PI–LL",
    y = "Patients",
    title = "Distribution of Planned \u0394LL",
    subtitle = sprintf("|Planned \u0394LL| \u2264 %d\u00b0 | n = %d | Red = mean", MAX_DEG, n_pat)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

pr_dir <- ensure_planned_results_dir()

ggsave(
  file.path(pr_dir, "analysis22_PI_LL_baseline_planned_trajectory.png"),
  plot = p_traj,
  width = 9,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(pr_dir, "analysis22_PI_LL_baseline_planned_delta_histogram.png"),
  plot = p_hist,
  width = 9,
  height = 5.5,
  dpi = 300
)

cat("\nSaved planned_results/analysis22_PI_LL_baseline_planned_trajectory.png\n")
cat("Saved planned_results/analysis22_PI_LL_baseline_planned_delta_histogram.png\n")
