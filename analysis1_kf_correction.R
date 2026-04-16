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

# Preop PI–LL mismatch (|LATpre_PI_LL|) must exceed this (°) to be included; NULL = no restriction
PREOP_PILL_MIN <- 10

# Load database
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

# Calculate change in lordosis
# Note: If either LAT1Y_L1_S1 or LATpre_L1_S1 is missing, change_lordosis will be NA
df$change_lordosis <- df$LAT1Y_L1_S1 - df$LATpre_L1_S1

# Drop patients with missing data (listwise deletion)
# Patients with missing LATpre_LL_KneeAngle or change_lordosis are excluded from analysis
df_clean <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_lordosis))

cat(sprintf(
  "After complete KF + 1Y lordosis change: n = %d (dropped %d)\n",
  nrow(df_clean),
  nrow(df) - nrow(df_clean)
))

if (!is.null(PREOP_PILL_MIN)) {
  if (!"LATpre_PI_LL" %in% names(df_clean)) {
    stop("COMBINE data must contain LATpre_PI_LL when PREOP_PILL_MIN is set")
  }
  n_before <- nrow(df_clean)
  df_clean <- df_clean %>%
    filter(!is.na(LATpre_PI_LL), abs(LATpre_PI_LL) > PREOP_PILL_MIN)
  cat(sprintf(
    "Restricted to |preop PI–LL| > %d°: n = %d (excluded %d)\n",
    PREOP_PILL_MIN,
    nrow(df_clean),
    n_before - nrow(df_clean)
  ))
}

cat(sprintf("Analysis sample: n = %d\n", nrow(df_clean)))

# Perform linear regression
model <- lm(change_lordosis ~ LATpre_LL_KneeAngle, data = df_clean)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

# Calculate Pearson correlation coefficient
pearson_r <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$change_lordosis, use = "complete.obs")

# Create scatter plot with regression line
p <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = change_lordosis)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
    y = "Change in Lordosis (LAT1Y_L1_S1 - LATpre_L1_S1)",
    title = "Preoperative Knee Flexion vs Change in Lordosis",
    subtitle = if (is.null(PREOP_PILL_MIN)) {
      NULL
    } else {
      sprintf("|preop PI–LL| > %d°", PREOP_PILL_MIN)
    }
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
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
out_png <- if (is.null(PREOP_PILL_MIN)) {
  "planned_results/analysis1_kf_correction.png"
} else {
  sprintf("planned_results/analysis1_kf_correction_preopPILLgt%d.png", PREOP_PILL_MIN)
}
ggsave(out_png, plot = p, width = 10, height = 8, dpi = 300)
cat(sprintf("Saved plot to %s\n", out_png))

