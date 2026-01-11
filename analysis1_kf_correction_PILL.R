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

# Calculate change in PI-LL mismatch
# Note: If either LAT1Y_PI_LL or LATpre_PI_LL is missing, change_PI_LL will be NA
df$change_PI_LL <- df$LAT1Y_PI_LL - df$LATpre_PI_LL

# Drop patients with missing data (listwise deletion)
# Patients with missing LATpre_LL_KneeAngle or change_PI_LL are excluded from analysis
df_clean <- df %>%
  filter(!is.na(LATpre_LL_KneeAngle) & !is.na(change_PI_LL))

cat(paste("Analysis includes", nrow(df_clean), "patients (dropped", nrow(df) - nrow(df_clean), "patients with missing data)\n"))

# Perform linear regression
model <- lm(change_PI_LL ~ LATpre_LL_KneeAngle, data = df_clean)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

# Calculate Pearson correlation coefficient
pearson_r <- cor(df_clean$LATpre_LL_KneeAngle, df_clean$change_PI_LL, use = "complete.obs")

# Create scatter plot with regression line
p <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = change_PI_LL)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative Knee Flexion (LATpre_LL_KneeAngle)",
    y = "Change in PI-LL Mismatch (LAT1Y_PI_LL - LATpre_PI_LL)",
    title = "Preoperative Knee Flexion vs Change in PI-LL Mismatch"
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
if (!dir.exists("results")) {
  dir.create("results")
}
ggsave("results/analysis1_kf_correction_PILL.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Saved plot to results/analysis1_kf_correction_PILL.png\n")

