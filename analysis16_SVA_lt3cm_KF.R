#!/usr/bin/env Rscript
# Analysis 16: Among patients with |6-week SVA| < 3 cm, knee flexion relationships
# Mirrors the logic of analysis2_correctedLL_kf.R (satisfactory PI-LL) for SVA balance.
#
# SVA is stored in mm in CADS (consistent with change_SVA <- (...)/10 for cm elsewhere).
# Threshold: |LAT6W_SVA_C2_S1| < 30 mm = 3 cm.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
SVA_THRESHOLD_MM <- 30

db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"
df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)

df_filtered <- df %>%
  filter(!is.na(LAT6W_SVA_C2_S1), abs(LAT6W_SVA_C2_S1) < SVA_THRESHOLD_MM)

cat(sprintf(
  "Filtered to %d patients with |6-week SVA| < %d mm (%.1f cm)\n",
  nrow(df_filtered), SVA_THRESHOLD_MM, SVA_THRESHOLD_MM / 10
))

df_filtered$knee_angle_diff <- df_filtered$LAT6W_LL_KneeAngle - df_filtered$LATpre_LL_KneeAngle
df_filtered$sva_6w_cm <- df_filtered$LAT6W_SVA_C2_S1 / 10

safe_cor_lm <- function(x, y, label) {
  ok <- !is.na(x) & !is.na(y)
  n <- sum(ok)
  if (n < 3) {
    cat(sprintf("%s: insufficient data (n = %d)\n", label, n))
    return(invisible(NULL))
  }
  r <- cor(x[ok], y[ok], use = "complete.obs")
  fit <- lm(y[ok] ~ x[ok])
  sm <- summary(fit)
  p <- sm$coefficients[2, 4]
  r2 <- sm$r.squared
  cat(sprintf("%s: n = %d, r = %.4f, R² = %.4f, p = %.4e\n", label, n, r, r2, p))
  invisible(list(n = n, r = r, r2 = r2, p = p))
}

cat("\n=== Key associations (same estimands as PI-LL <10° writeup) ===\n")
cat("1) 6-week KF vs preoperative KF\n")
safe_cor_lm(
  df_filtered$LATpre_LL_KneeAngle,
  df_filtered$LAT6W_LL_KneeAngle,
  "LAT6W_LL_KneeAngle ~ LATpre_LL_KneeAngle"
)

cat("\n2) 6-week KF vs change in KF (6W − preop)\n")
safe_cor_lm(
  df_filtered$knee_angle_diff,
  df_filtered$LAT6W_LL_KneeAngle,
  "LAT6W_LL_KneeAngle ~ knee_angle_diff"
)

cat("\n2b) Preoperative KF vs change in KF (alternative estimand)\n")
safe_cor_lm(
  df_filtered$LATpre_LL_KneeAngle,
  df_filtered$knee_angle_diff,
  "knee_angle_diff ~ LATpre_LL_KneeAngle"
)

cat("\n3) ΔKF vs 6-week SVA (cm) within this |SVA|<3 cm cohort\n")
safe_cor_lm(
  df_filtered$sva_6w_cm,
  df_filtered$knee_angle_diff,
  "knee_angle_diff ~ SVA_6w_cm"
)

create_plot <- function(data, x_var, y_var, x_label, y_label, title, filename) {
  data_clean <- data %>%
    filter(!is.na(!!sym(x_var)) & !is.na(!!sym(y_var)))

  if (nrow(data_clean) == 0) {
    cat(paste("Warning: No valid data for", title, "\n"))
    return(NULL)
  }

  formula <- as.formula(paste(y_var, "~", x_var))
  model <- lm(formula, data = data_clean)
  summary_model <- summary(model)
  r_squared <- summary_model$r.squared
  p_value <- summary_model$coefficients[2, 4]
  pearson_r <- cor(data_clean[[x_var]], data_clean[[y_var]], use = "complete.obs")

  p <- ggplot(data_clean, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0(
        "r = ", round(pearson_r, 3),
        "\nR² = ", round(r_squared, 3),
        "\np = ", formatC(p_value, format = "e", digits = 2)
      ),
      hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold"
    )

  if (!dir.exists("results")) dir.create("results")
  filepath <- file.path("results", filename)
  ggsave(filepath, plot = p, width = 10, height = 8, dpi = 300)
  cat(paste("Saved plot to", filepath, "\n"))
  invisible(p)
}

subtitle <- sprintf("|6-week SVA| < %.1f cm", SVA_THRESHOLD_MM / 10)

create_plot(
  df_filtered,
  "sva_6w_cm",
  "LATpre_LL_KneeAngle",
  "6-week SVA (cm)",
  "Preoperative knee flexion",
  paste0("6-week SVA vs preoperative KF\n(", subtitle, ")"),
  "analysis16_SVA_vs_preop_KF.png"
)

create_plot(
  df_filtered,
  "sva_6w_cm",
  "LAT6W_LL_KneeAngle",
  "6-week SVA (cm)",
  "6-week knee flexion",
  paste0("6-week SVA vs 6-week KF\n(", subtitle, ")"),
  "analysis16_SVA_vs_6W_KF.png"
)

create_plot(
  df_filtered,
  "sva_6w_cm",
  "knee_angle_diff",
  "6-week SVA (cm)",
  "Change in knee flexion (6-week − preop)",
  paste0("6-week SVA vs ΔKF\n(", subtitle, ")"),
  "analysis16_SVA_vs_delta_KF.png"
)

cat("\n=== Analysis 16 complete ===\n")
