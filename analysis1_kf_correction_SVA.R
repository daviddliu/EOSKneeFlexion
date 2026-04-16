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

# Planned ΔSVA (cm): drop non-physical / data-entry tails (main cluster within a few tens of cm)
MAX_PLANNED_DSVA_ABS_CM <- 100

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)

# Calculate change in SVA (convert from mm to cm by dividing by 10)
# Note: If either LAT1Y_SVA_C2_S1 or LATpre_SVA_C2_S1 is missing, change_SVA will be NA
df$change_SVA <- (df$LAT1Y_SVA_C2_S1 - df$LATpre_SVA_C2_S1) / 10

# Planned ΔSVA (cm): preop C2–S1 SVA (mm → cm) minus parsed planned SVA from demo_AlignGoal_SVA
# (same cm scale as change_SVA above; planned_SVA from utils/planned_pi_ll.R)
df$planned_dSVA <- (df$LATpre_SVA_C2_S1 / 10) - df$planned_SVA

# Drop patients with missing data (listwise deletion)
# Patients with missing LATpre_LL_KneeAngle or change_SVA are excluded from analysis
df_clean <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_SVA))

cat(paste("Analysis includes", nrow(df_clean), "patients (dropped", nrow(df) - nrow(df_clean), "patients with missing data)\n"))

# Perform linear regression
model <- lm(change_SVA ~ LATpre_LL_KneeAngle, data = df_clean)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

# Calculate Pearson correlation coefficient
pearson_r <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$change_SVA, use = "complete.obs")

# Create scatter plot with regression line
p <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = change_SVA)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
    y = "Change in SVA (cm) ((LAT1Y_SVA_C2_S1 - LATpre_SVA_C2_S1) / 10)",
    title = "Preoperative Knee Flexion vs Change in SVA"
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

# Save plot
if (!dir.exists("planned_results")) {
  dir.create("planned_results")
}
ggsave("planned_results/analysis1_kf_correction_SVA.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Saved plot to planned_results/analysis1_kf_correction_SVA.png\n")

# --- Knee flexion vs planned ΔSVA (preop SVA − planned SVA, cm) ---
df_planned_complete <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(planned_dSVA))

df_planned <- df_planned_complete %>%
  filter(abs(planned_dSVA) <= MAX_PLANNED_DSVA_ABS_CM)

n_planned_out <- nrow(df_planned_complete) - nrow(df_planned)
if (n_planned_out > 0L) {
  cat(sprintf(
    "Excluded %d patient(s) with |planned \u0394SVA| > %d cm\n",
    n_planned_out,
    MAX_PLANNED_DSVA_ABS_CM
  ))
}

cat(paste(
  "Planned ΔSVA analysis includes", nrow(df_planned), "patients (dropped",
  nrow(df) - nrow(df_planned_complete), "with missing KF or preop/planned SVA)\n"
))

if (nrow(df_planned) >= 3L) {
  model_p <- lm(planned_dSVA ~ LATpre_LL_KneeAngle, data = df_planned)
  summary_p <- summary(model_p)
  r_squared_p <- summary_p$r.squared
  p_value_p <- summary_p$coefficients[2, 4]
  pearson_r_p <- cor(df_planned$LATpre_LL_KneeAngle, df_planned$planned_dSVA, use = "complete.obs")

  p_planned <- ggplot(df_planned, aes(x = LATpre_LL_KneeAngle, y = planned_dSVA)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linetype = "solid") +
    labs(
      x = "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
      y = "Planned \u0394SVA (cm) (LATpre_SVA_C2_S1/10 \u2212 planned_SVA)",
      title = "Preoperative Knee Flexion vs Planned \u0394SVA"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ann_p <- paste0(
    "r = ", round(pearson_r_p, 3), "\n",
    "R² = ", round(r_squared_p, 3), "\n",
    "p = ", formatC(p_value_p, format = "e", digits = 2)
  )
  p_planned <- p_planned + annotate(
    "text",
    x = Inf, y = Inf,
    label = ann_p,
    hjust = 1.1, vjust = 1.5,
    size = 4,
    fontface = "bold"
  )

  ggsave(
    "planned_results/analysis1_kf_correction_planned_dSVA.png",
    plot = p_planned,
    width = 10,
    height = 8,
    dpi = 300
  )
  cat("Saved plot to planned_results/analysis1_kf_correction_planned_dSVA.png\n")
} else {
  cat("Skipped planned \u0394SVA plot: fewer than 3 patients with complete data.\n")
}

