#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
})

options(warn = -1)

# Source utility functions
source("utils/utils.R")

# Configuration: Toggle PJK exclusion on/off
EXCLUDE_PJK <- TRUE

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)
df <- restrict_planned_dll_analysis_cohort(df)

# Outcome: Planned ΔLL (preop − planned PI–LL goal)
df_clean <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle),
    !is.na(planned_DLL),
    planned_DLL >= PLANNED_DLL_MIN_KEEP
  )

cat(sprintf(
  "Cohort: |preop PI–LL| > %d°, preop SVA C2–C7 < %d (non-missing)\n",
  PREOP_ABS_PI_LL_GT,
  PREOP_SVA_C2_C7_LT
))
cat(paste("Analysis includes", nrow(df_clean), "patients (dropped", nrow(df) - nrow(df_clean), "patients with missing KF / planned goal / truncation)\n"))

# Perform linear regression
model <- lm(planned_DLL ~ LATpre_LL_KneeAngle, data = df_clean)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

# Calculate Pearson correlation coefficient
pearson_r <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$planned_DLL, use = "complete.obs")

# Create scatter plot with regression line
p <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = planned_DLL)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative Knee Flexion (°)",
    y = "Planned \u0394LL (\u00b0)",
    title = "Preoperative Knee Flexion vs Planned \u0394LL"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
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

pr_dir <- ensure_planned_results_dir()
ggsave(
  file.path(pr_dir, "analysis1_kf_correction_PILL.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)
cat("Saved plot to planned_results/analysis1_kf_correction_PILL.png\n")
