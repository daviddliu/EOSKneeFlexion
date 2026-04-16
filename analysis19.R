#!/usr/bin/env Rscript
#
# Analysis 19: Preoperative knee flexion vs Planned ΔLL
#
# Planned ΔLL = −(planned_PI_LL − LATpre_PI_LL) = LATpre_PI_LL − planned_PI_LL
#   (flip sign of planned − preop PI–LL), after recoding small nonnegative planned
#   goals in [0, 10]° to 0 for analysis.
# Exclusions: |LATpre_PI_LL| > PREOP_ABS_PI_LL_GT; preop SVA C2–C7 < PREOP_SVA_C2_C7_LT;
#   Planned ΔLL >= PLANNED_DLL_MIN_KEEP.
#
# Plus PCA adjustment (BIC-chosen PCs) on preop radiographs + age; KF and
# LATpre_PI_LL are excluded from the PC space.
#
# Mild vs cervical-deformity strata (no global SVA C2–C7 exclusion): analysis19_cervical_deformity.R

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

if ("LATpre_L1PA" %in% names(df) && "LATpre_T4PA" %in% names(df)) {
  df$LATpre_T4_L1_PA <- df$LATpre_L1PA - df$LATpre_T4PA
} else {
  df$LATpre_T4_L1_PA <- NA_real_
}

age_var <- NULL
if ("demo_Age" %in% names(df)) {
  age_var <- "demo_Age"
} else if ("demo_age" %in% names(df)) {
  age_var <- "demo_age"
} else if ("age" %in% names(df)) {
  age_var <- "age"
} else if ("Age" %in% names(df)) {
  age_var <- "Age"
}

get_preop_confounders_no_pill <- function() {
  v <- c(
    "LATpre_S1PI", "LATpre_L1_S1", "LATpre_L4_S1", "LATpre_T2_T12",
    "LATpre_S1PT", "LATpre_SVA_C2_S1", "LATpre_T4_L1_PA"
  )
  if (!is.null(age_var)) {
    v <- c(v, age_var)
  }
  v
}

pca_adjust_planned_DLL <- function(dat, confounder_vars) {
  x_vec <- as.numeric(dat[["LATpre_LL_KneeAngle"]])
  y_vec <- as.numeric(dat[["planned_DLL"]])
  avail <- confounder_vars[confounder_vars %in% names(dat)]

  cc <- !is.na(x_vec) & !is.na(y_vec)
  for (v in avail) {
    cc <- cc & !is.na(dat[[v]])
  }
  x_clean <- x_vec[cc]
  y_clean <- y_vec[cc]
  if (length(x_clean) < 5) {
    return(NULL)
  }

  X <- as.matrix(dat[cc, avail, drop = FALSE])
  p <- ncol(X)
  n <- nrow(X)
  if (n < p + 3) {
    return(NULL)
  }

  Xs <- scale(X)
  pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
  pcx <- pca$x
  max_k <- min(p, n - 2L, 15L)
  if (max_k < 1L) {
    return(NULL)
  }

  best_bic <- Inf
  best_k <- 1L
  best_fit <- NULL

  for (k in seq_len(max_k)) {
    dt <- data.frame(y = y_clean, x = x_clean)
    for (i in seq_len(k)) {
      dt[[paste0("PC", i)]] <- pcx[, i]
    }
    f <- as.formula(paste("y ~ x +", paste(paste0("PC", seq_len(k)), collapse = " + ")))
    fit <- lm(f, data = dt)
    bic <- BIC(fit)
    if (bic < best_bic) {
      best_bic <- bic
      best_k <- k
      best_fit <- fit
    }
  }

  if (is.null(best_fit)) {
    return(NULL)
  }

  dt_best <- data.frame(y = y_clean, x = x_clean)
  for (i in seq_len(best_k)) {
    dt_best[[paste0("PC", i)]] <- pcx[, i]
  }
  sm <- summary(best_fit)
  slope <- sm$coefficients[2L, 1L]
  p_x <- sm$coefficients[2L, 4L]
  t_x <- sm$coefficients[2L, 3L]
  df_r <- sm$df[2L]
  r_partial <- t_x / sqrt(t_x^2 + df_r)

  f_pc_only <- as.formula(paste("y ~", paste(paste0("PC", seq_len(best_k)), collapse = " + ")))
  r2_pc <- summary(lm(f_pc_only, data = dt_best))$r.squared
  r2_full <- sm$r.squared
  r2_incr <- r2_full - r2_pc

  y_hat_pc <- predict(lm(f_pc_only, data = dt_best))
  partial_resid <- y_clean - y_hat_pc

  list(
    n = n,
    n_components = best_k,
    r_partial = unname(r_partial),
    r2_increment_x = unname(r2_incr),
    p_x = unname(p_x),
    slope_x = unname(slope),
    plot_df = tibble::tibble(
      LATpre_LL_KneeAngle = x_clean,
      partial_residual_planned_DLL = partial_resid
    )
  )
}

MAX_DEG <- 100

demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
}
if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
  demo_id_col <- "demo_id"
}

n_raw <- sum(!is.na(df$planned_PI_LL_parsed))
n_recode <- sum(
  !is.na(df$planned_PI_LL_parsed) &
    df$planned_PI_LL_parsed >= 0 &
    df$planned_PI_LL_parsed <= PLANNED_PILL_RECODE_MAX,
  na.rm = TRUE
)
n_lt <- sum(grepl("^<", trimws(df$demo_AlignGoal_PILL_raw), ignore.case = TRUE), na.rm = TRUE)
cat("\n=== Analysis 19: KF vs Planned DLL (preop minus planned PI-LL) ===\n")
cat(sprintf(
  "Planned DLL = preop − planned; exclude Planned ΔLL < %g°.\n",
  PLANNED_DLL_MIN_KEEP
))
cat(sprintf(
  paste0(
    "Inclusion: |preop PI–LL| > %d° (non-missing); preop SVA C2–C7 < %d (non-missing); ",
    "Planned ΔLL ≥ %g°\n"
  ),
  PREOP_ABS_PI_LL_GT,
  PREOP_SVA_C2_C7_LT,
  PLANNED_DLL_MIN_KEEP
))
cat(sprintf("Patients with parsed planned PI–LL: %d\n", n_raw))
cat(sprintf(
  "Parsed planned in [0, %d]° recoded to 0 for analysis: %d patient-rows\n",
  PLANNED_PILL_RECODE_MAX,
  n_recode
))
cat(sprintf("Rows with leading '<' in raw goal (parsed as bound or default 10): ~%d\n", n_lt))

df_complete <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle),
    !is.na(LATpre_PI_LL),
    !is.na(planned_PI_LL_parsed),
    !is.na(planned_DLL),
    planned_DLL >= PLANNED_DLL_MIN_KEEP
  )

n_drop_dll <- sum(
  !is.na(df$LATpre_LL_KneeAngle) &
    !is.na(df$LATpre_PI_LL) &
    !is.na(df$planned_PI_LL_parsed) &
    !is.na(df$planned_DLL) &
    df$planned_DLL < PLANNED_DLL_MIN_KEEP,
  na.rm = TRUE
)
if (n_drop_dll > 0) {
  cat(sprintf(
    "Excluded %d row(s) with Planned ΔLL < %g° (before |PI–LL| > %d° restriction)\n",
    n_drop_dll,
    PLANNED_DLL_MIN_KEEP,
    PREOP_ABS_PI_LL_GT
  ))
}

cat(sprintf(
  "Complete cases (KF + preop + planned + Planned ΔLL ≥ %g°): n = %d\n",
  PLANNED_DLL_MIN_KEEP,
  nrow(df_complete)
))

df_restrict <- restrict_planned_dll_analysis_cohort(df_complete)

cat(sprintf(
  paste0(
    "Restricted to |preop PI–LL| > %d° and preop SVA C2–C7 < %d: n = %d ",
    "(excluded %d from complete-case pool)\n"
  ),
  PREOP_ABS_PI_LL_GT,
  PREOP_SVA_C2_C7_LT,
  nrow(df_restrict),
  nrow(df_complete) - nrow(df_restrict)
))

df_clean <- df_restrict %>%
  filter(
    abs(LATpre_LL_KneeAngle) <= MAX_DEG,
    abs(LATpre_PI_LL) <= MAX_DEG,
    abs(planned_PI_LL) <= MAX_DEG,
    abs(planned_DLL) <= MAX_DEG
  )

n_outlier <- nrow(df_restrict) - nrow(df_clean)
if (n_outlier > 0) {
  cat(sprintf(
    "Excluded %d patient(s) with any KF / PI–LL / planned / |Planned ΔLL| > %d°\n",
    n_outlier,
    MAX_DEG
  ))
}

cat(sprintf(
  "After |value| ≤ %d° outlier screen: n = %d\n",
  MAX_DEG,
  nrow(df_clean)
))
cat(sprintf("Analysis sample: n = %d\n", nrow(df_clean)))

DIST_REVIEW_LT <- 5
cat("\n=== Review: small |Planned ΔLL| ===\n")
cat(sprintf(
  "Patients in analysis sample with |Planned ΔLL| < %d°: ",
  DIST_REVIEW_LT
))
sub_lt <- df_clean %>% filter(abs(planned_DLL) < DIST_REVIEW_LT)
cat(sprintf(
  "n = %d of %d (%.1f%%)\n",
  nrow(sub_lt),
  nrow(df_clean),
  if (nrow(df_clean) > 0) 100 * nrow(sub_lt) / nrow(df_clean) else NA_real_
))
if (nrow(sub_lt) > 0) {
  review_tbl <- sub_lt %>%
    mutate(demo_id = .data[[demo_id_col]]) %>%
    transmute(
      demo_id,
      LATpre_PI_LL,
      planned_PI_LL_parsed,
      planned_PI_LL_analysis = planned_PI_LL,
      demo_AlignGoal_PILL_raw,
      planned_minus_preop_PI_LL,
      planned_DLL,
      LATpre_LL_KneeAngle
    ) %>%
    arrange(abs(planned_DLL), demo_id)
  print(as.data.frame(review_tbl), row.names = FALSE, digits = 4)
  readr::write_csv(
    review_tbl,
    file.path(ensure_planned_results_dir(), "analysis19_planned_DLL_abs_lt5deg_patients.csv")
  )
  cat("\nFull table: planned_results/analysis19_planned_DLL_abs_lt5deg_patients.csv\n")
}

model <- lm(planned_DLL ~ LATpre_LL_KneeAngle, data = df_clean)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]
pearson_r <- cor(
  df_clean$LATpre_LL_KneeAngle,
  df_clean$planned_DLL,
  use = "complete.obs"
)

cat(sprintf("Pearson r: %.4f\n", pearson_r))
cat(sprintf("R²: %.4f\n", r_squared))
cat(sprintf("p (slope): %.4e\n", p_value))

conf_preop <- get_preop_confounders_no_pill()
conf_avail <- conf_preop[conf_preop %in% names(df_clean)]
cat("\n=== PCA confounder adjustment (preop radiographs + age) ===\n")
cat("PC space excludes KF and LATpre_PI_LL (PI–LL defines part of the outcome).\n")
cat(sprintf("Confounders in PCA: %s\n", paste(conf_avail, collapse = ", ")))

pca_res <- pca_adjust_planned_DLL(df_clean, conf_preop)
if (is.null(pca_res)) {
  cat("PCA adjustment skipped (insufficient complete cases or rank).\n\n")
} else {
  cat(sprintf("Complete cases for PCA model: n = %d\n", pca_res$n))
  cat(sprintf("BIC-selected PCs: %d\n", pca_res$n_components))
  cat(sprintf("Partial r (KF vs outcome | PCs): %.4f\n", pca_res$r_partial))
  cat(sprintf("Incremental R² for KF after PCs: %.4f\n", pca_res$r2_increment_x))
  cat(sprintf("Slope (KF | PCs): %.6f\n", pca_res$slope_x))
  cat(sprintf("p (KF | PCs): %.4e\n\n", pca_res$p_x))
}

p <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = planned_DLL)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative knee flexion (LATpre_LL_KneeAngle, °)",
    y = "Planned \u0394LL (\u00b0)",
    title = "Preoperative knee flexion vs Planned \u0394LL",
    subtitle = sprintf(
      "|preop PI–LL| > %d\u00b0 | SVA C2–C7 < %d | Planned \u0394LL \u2265 %g\u00b0",
      PREOP_ABS_PI_LL_GT,
      PREOP_SVA_C2_C7_LT,
      PLANNED_DLL_MIN_KEEP
    )
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

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
  file.path(pr_dir, "analysis19_kf_vs_planned_DLL.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)
cat("Saved plot to planned_results/analysis19_kf_vs_planned_DLL.png\n")

if (!is.null(pca_res)) {
  ann_pca <- paste0(
    "Partial r = ", round(pca_res$r_partial, 3), "\n",
    "Partial R\u00b2 = ", round(pca_res$r2_increment_x, 3), "\n",
    "Coefficient = ", sprintf("%.4f", pca_res$slope_x), "\n",
    "p = ", formatC(pca_res$p_x, format = "e", digits = 2)
  )
  p_pca_plot <- ggplot(
    pca_res$plot_df,
    aes(x = LATpre_LL_KneeAngle, y = partial_residual_planned_DLL)
  ) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", linetype = "solid") +
    labs(
      x = "Preoperative knee flexion (\u00b0)",
      y = "Residuals in planned ∆LL (°) from PCA-based model",
      title = "Residual relationship in knee flexion with planned ∆LL: \n Variance subtracted from PCA-based model)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = ann_pca,
      hjust = 1.05,
      vjust = 1.15,
      size = 4,
      fontface = "bold"
    )
  ggsave(
    file.path(pr_dir, "analysis19_kf_vs_planned_DLL_PCA_partial_residual.png"),
    plot = p_pca_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  cat("Saved plot to planned_results/analysis19_kf_vs_planned_DLL_PCA_partial_residual.png\n")
}

# -----------------------------------------------------------------------------
# SVA C2–C7 strata: separate unadjusted and PCA-adjusted plots (< 40 vs > 40)
# Same complete-case pool as df_complete, PI–LL rule only (no SVA restriction),
# then split on non-missing LATpre_SVA_C2_C7; same |value| ≤ MAX_DEG screen.
# -----------------------------------------------------------------------------
cat("\n=== SVA C2–C7 stratum plots (unadjusted + PCA) ===\n")

df_pi_for_strata <- df_complete %>%
  filter(!is.na(LATpre_PI_LL), abs(LATpre_PI_LL) > PREOP_ABS_PI_LL_GT)

stratum_screen <- function(d) {
  filter(
    d,
    abs(LATpre_LL_KneeAngle) <= MAX_DEG,
    abs(LATpre_PI_LL) <= MAX_DEG,
    abs(planned_PI_LL) <= MAX_DEG,
    abs(planned_DLL) <= MAX_DEG
  )
}

strata_list <- list(
  list(
    id = "SVA_C2C7_lt40",
    lab = "preop SVA C2–C7 < 40",
    dat = stratum_screen(
      df_pi_for_strata %>%
        filter(!is.na(LATpre_SVA_C2_C7), LATpre_SVA_C2_C7 < 40)
    )
  ),
  list(
    id = "SVA_C2C7_gt40",
    lab = "preop SVA C2–C7 > 40",
    dat = stratum_screen(
      df_pi_for_strata %>%
        filter(!is.na(LATpre_SVA_C2_C7), LATpre_SVA_C2_C7 > 40)
    )
  )
)

for (st in strata_list) {
  d_str <- st$dat
  cat(sprintf("Stratum %s: n = %d\n", st$id, nrow(d_str)))

  if (nrow(d_str) < 3L) {
    cat(sprintf("  Skip plots (n < 3).\n"))
    next
  }

  sm_u <- summary(lm(planned_DLL ~ LATpre_LL_KneeAngle, data = d_str))
  r_u <- cor(d_str$LATpre_LL_KneeAngle, d_str$planned_DLL, use = "complete.obs")
  ann_u <- paste0(
    "r = ", round(r_u, 3), "\n",
    "R² = ", round(sm_u$r.squared, 3), "\n",
    "p = ", formatC(sm_u$coefficients[2, 4], format = "e", digits = 2)
  )

  sub_u <- sprintf(
    "|preop PI–LL| > %d\u00b0 | %s | Planned \u0394LL \u2265 %g\u00b0",
    PREOP_ABS_PI_LL_GT,
    st$lab,
    PLANNED_DLL_MIN_KEEP
  )

  p_u <- ggplot(d_str, aes(x = LATpre_LL_KneeAngle, y = planned_DLL)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = "Preoperative knee flexion (LATpre_LL_KneeAngle, \u00b0)",
      y = "Planned \u0394LL (\u00b0)",
      title = "KF vs Planned \u0394LL (unadjusted)",
      subtitle = sub_u
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotate("text", x = Inf, y = Inf, label = ann_u, hjust = 1.1, vjust = 1.5, size = 4, fontface = "bold")

  fn_u <- file.path(pr_dir, sprintf("analysis19_kf_vs_planned_DLL_unadjusted_%s.png", st$id))
  ggsave(fn_u, plot = p_u, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved unadjusted: %s\n", fn_u))

  pca_s <- pca_adjust_planned_DLL(d_str, conf_preop)
  if (is.null(pca_s)) {
    cat("  PCA stratum plot skipped (insufficient complete cases or rank).\n")
    next
  }

  pr_s <- cor(
    pca_s$plot_df$LATpre_LL_KneeAngle,
    pca_s$plot_df$partial_residual_planned_DLL,
    use = "complete.obs"
  )
  sub_pca <- sprintf(
    paste0(
      "%s | |preop PI\u2013LL| > %d\u00b0 | Planned \u0394LL \u2265 %g\u00b0 | ",
      "r (partial resid vs KF) = %.3f | p (KF | PCs) = %s"
    ),
    st$lab,
    PREOP_ABS_PI_LL_GT,
    PLANNED_DLL_MIN_KEEP,
    pr_s,
    formatC(pca_s$p_x, format = "e", digits = 2)
  )

  p_pca_s <- ggplot(
    pca_s$plot_df,
    aes(x = LATpre_LL_KneeAngle, y = partial_residual_planned_DLL)
  ) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", linetype = "solid") +
    labs(
      x = "Preoperative knee flexion (\u00b0)",
      y = "Partial residual: Planned \u0394LL | PCs",
      title = "KF vs Planned \u0394LL (PCA-adjusted confounders)",
      subtitle = sub_pca
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  fn_p <- file.path(pr_dir, sprintf("analysis19_kf_vs_planned_DLL_PCA_%s.png", st$id))
  ggsave(fn_p, plot = p_pca_s, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved PCA-adjusted: %s\n", fn_p))
}
