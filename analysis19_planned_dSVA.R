#!/usr/bin/env Rscript
#
# Analysis 19 (SVA analog): Preoperative knee flexion vs planned change in SVA
#
# Planned ΔSVA (cm) = (LATpre_SVA_C2_S1 / 10) − planned_SVA (same as analysis1_kf_correction_SVA.R).
# Cohort parallels analysis19.R: same PI–LL / SVA C2–C7 restriction and knee–angle caps;
#   replaces Planned ΔLL screens with |planned_dSVA| ≤ MAX_PLANNED_DSVA_ABS_CM.
#
# PCA: BIC-chosen PCs on preop radiographs + age; KF and LATpre_PI_LL excluded from PC space
#   (same confounder set as analysis19).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

MAX_PLANNED_DSVA_ABS_CM <- 100

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)

df$planned_dSVA <- (df$LATpre_SVA_C2_S1 / 10) - df$planned_SVA

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

pca_adjust_planned_dSVA <- function(dat, confounder_vars) {
  x_vec <- as.numeric(dat[["LATpre_LL_KneeAngle"]])
  y_vec <- as.numeric(dat[["planned_dSVA"]])
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
      partial_residual_planned_dSVA = partial_resid
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

cat("\n=== Analysis 19 (planned ΔSVA): KF vs planned change in SVA (cm) ===\n")
cat(sprintf(
  "Planned ΔSVA = preop C2–S1 SVA (cm) − planned SVA (cm); |Δ| ≤ %g cm.\n",
  MAX_PLANNED_DSVA_ABS_CM
))
cat(sprintf(
  paste0(
    "Inclusion (then restrict_planned_dll): |preop PI–LL| > %d°; preop SVA C2–C7 < %d; ",
    "non-missing KF, PILL parse, preop/planned SVA\n"
  ),
  PREOP_ABS_PI_LL_GT,
  PREOP_SVA_C2_C7_LT
))

df_complete <- df %>%
  filter(
    !is.na(LATpre_LL_KneeAngle),
    !is.na(LATpre_PI_LL),
    !is.na(planned_PI_LL_parsed),
    !is.na(LATpre_SVA_C2_S1),
    !is.na(planned_SVA),
    !is.na(planned_dSVA),
    abs(planned_dSVA) <= MAX_PLANNED_DSVA_ABS_CM
  )

cat(sprintf(
  "Complete cases (KF + PI–LL + PILL parse + SVA + |planned ΔSVA| ≤ %g cm): n = %d\n",
  MAX_PLANNED_DSVA_ABS_CM,
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
    abs(planned_dSVA) <= MAX_PLANNED_DSVA_ABS_CM
  )

n_outlier <- nrow(df_restrict) - nrow(df_clean)
if (n_outlier > 0) {
  cat(sprintf(
    "Excluded %d patient(s) with KF / PI–LL / planned PI–LL or |planned ΔSVA| out of range\n",
    n_outlier
  ))
}

cat(sprintf(
  "After |KF|, |PI–LL|, |planned PI–LL| ≤ %d° and |planned ΔSVA| ≤ %g cm: n = %d\n",
  MAX_DEG,
  MAX_PLANNED_DSVA_ABS_CM,
  nrow(df_clean)
))

DIST_REVIEW_LT <- 5
cat(sprintf(
  "\n=== Review: small |planned ΔSVA| (< %g cm) ===\n",
  DIST_REVIEW_LT
))
sub_lt <- df_clean %>% filter(abs(planned_dSVA) < DIST_REVIEW_LT)
cat(sprintf("n = %d of %d\n", nrow(sub_lt), nrow(df_clean)))
if (nrow(sub_lt) > 0) {
  review_tbl <- sub_lt %>%
    mutate(demo_id = .data[[demo_id_col]]) %>%
    transmute(
      demo_id,
      LATpre_PI_LL,
      planned_PI_LL_parsed,
      LATpre_SVA_C2_S1_mm = LATpre_SVA_C2_S1,
      planned_SVA_cm = planned_SVA,
      planned_dSVA,
      demo_AlignGoal_SVA_raw,
      LATpre_LL_KneeAngle
    ) %>%
    arrange(abs(planned_dSVA), demo_id)
  print(as.data.frame(review_tbl), row.names = FALSE, digits = 4)
  readr::write_csv(
    review_tbl,
    file.path(ensure_planned_results_dir(), "analysis19_planned_dSVA_abs_lt5cm_patients.csv")
  )
  cat("Full table: planned_results/analysis19_planned_dSVA_abs_lt5cm_patients.csv\n")
}

model <- lm(planned_dSVA ~ LATpre_LL_KneeAngle, data = df_clean)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]
pearson_r <- cor(
  df_clean$LATpre_LL_KneeAngle,
  df_clean$planned_dSVA,
  use = "complete.obs"
)

cat(sprintf("Pearson r: %.4f\n", pearson_r))
cat(sprintf("R²: %.4f\n", r_squared))
cat(sprintf("p (slope): %.4e\n", p_value))

conf_preop <- get_preop_confounders_no_pill()
conf_avail <- conf_preop[conf_preop %in% names(df_clean)]
cat("\n=== PCA confounder adjustment (preop radiographs + age) ===\n")
cat("PC space excludes KF and LATpre_PI_LL (same as analysis19).\n")
cat(sprintf("Confounders in PCA: %s\n", paste(conf_avail, collapse = ", ")))

pca_res <- pca_adjust_planned_dSVA(df_clean, conf_preop)
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

pr_dir <- ensure_planned_results_dir()

p <- ggplot(df_clean, aes(x = LATpre_LL_KneeAngle, y = planned_dSVA)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
  labs(
    x = "Preoperative knee flexion (LATpre_LL_KneeAngle, °)",
    y = "Planned ΔSVA (cm)",
    title = "Preoperative knee flexion vs planned ΔSVA",
    subtitle = sprintf(
      "|preop PI–LL| > %d° | SVA C2–C7 < %d | |planned ΔSVA| ≤ %g cm",
      PREOP_ABS_PI_LL_GT,
      PREOP_SVA_C2_C7_LT,
      MAX_PLANNED_DSVA_ABS_CM
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

ggsave(
  file.path(pr_dir, "analysis19_planned_dSVA_kf_vs_planned_dSVA.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)
cat("Saved plot to planned_results/analysis19_planned_dSVA_kf_vs_planned_dSVA.png\n")

if (!is.null(pca_res)) {
  ann_pca <- paste0(
    "Partial r = ", round(pca_res$r_partial, 3), "\n",
    "Partial R\u00b2 = ", round(pca_res$r2_increment_x, 3), "\n",
    "Coefficient = ", sprintf("%.4f", pca_res$slope_x), "\n",
    "p = ", formatC(pca_res$p_x, format = "e", digits = 2)
  )
  p_pca_plot <- ggplot(
    pca_res$plot_df,
    aes(x = LATpre_LL_KneeAngle, y = partial_residual_planned_dSVA)
  ) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", linetype = "solid") +
    labs(
      x = "Preoperative knee flexion (°)",
      y = "Residuals in planned ΔSVA (cm) from PCA-based model",
      title = "Residual relationship: knee flexion vs planned ΔSVA\n(variance subtracted via PCA confounders)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
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
    file.path(pr_dir, "analysis19_planned_dSVA_PCA_partial_residual.png"),
    plot = p_pca_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  cat("Saved plot to planned_results/analysis19_planned_dSVA_PCA_partial_residual.png\n")
}

# SVA C2–C7 strata (same logic as analysis19)
cat("\n=== SVA C2–C7 stratum plots (unadjusted + PCA) ===\n")

df_pi_for_strata <- df_complete %>%
  filter(!is.na(LATpre_PI_LL), abs(LATpre_PI_LL) > PREOP_ABS_PI_LL_GT)

stratum_screen <- function(d) {
  filter(
    d,
    abs(LATpre_LL_KneeAngle) <= MAX_DEG,
    abs(LATpre_PI_LL) <= MAX_DEG,
    abs(planned_PI_LL) <= MAX_DEG,
    abs(planned_dSVA) <= MAX_PLANNED_DSVA_ABS_CM
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

  sm_u <- summary(lm(planned_dSVA ~ LATpre_LL_KneeAngle, data = d_str))
  r_u <- cor(d_str$LATpre_LL_KneeAngle, d_str$planned_dSVA, use = "complete.obs")
  ann_u <- paste0(
    "r = ", round(r_u, 3), "\n",
    "R² = ", round(sm_u$r.squared, 3), "\n",
    "p = ", formatC(sm_u$coefficients[2, 4], format = "e", digits = 2)
  )

  sub_u <- sprintf(
    "|preop PI–LL| > %d° | %s | |planned ΔSVA| ≤ %g cm",
    PREOP_ABS_PI_LL_GT,
    st$lab,
    MAX_PLANNED_DSVA_ABS_CM
  )

  p_u <- ggplot(d_str, aes(x = LATpre_LL_KneeAngle, y = planned_dSVA)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") +
    labs(
      x = "Preoperative knee flexion (°)",
      y = "Planned ΔSVA (cm)",
      title = "KF vs planned ΔSVA (unadjusted)",
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

  fn_u <- file.path(pr_dir, sprintf("analysis19_planned_dSVA_unadjusted_%s.png", st$id))
  ggsave(fn_u, plot = p_u, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved unadjusted: %s\n", fn_u))

  pca_s <- pca_adjust_planned_dSVA(d_str, conf_preop)
  if (is.null(pca_s)) {
    cat("  PCA stratum plot skipped (insufficient complete cases or rank).\n")
    next
  }

  ann_s <- paste0(
    "Partial r = ", round(pca_s$r_partial, 3), "\n",
    "Partial R\u00b2 = ", round(pca_s$r2_increment_x, 3), "\n",
    "Coefficient = ", sprintf("%.4f", pca_s$slope_x), "\n",
    "p = ", formatC(pca_s$p_x, format = "e", digits = 2)
  )

  p_pca_s <- ggplot(
    pca_s$plot_df,
    aes(x = LATpre_LL_KneeAngle, y = partial_residual_planned_dSVA)
  ) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", linetype = "solid") +
    labs(
      x = "Preoperative knee flexion (°)",
      y = "Partial residual: planned ΔSVA | PCs (cm)",
      title = "KF vs planned ΔSVA (PCA-adjusted)",
      subtitle = sprintf("%s | n = %d", st$lab, nrow(d_str))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = ann_s,
      hjust = 1.05, vjust = 1.15, size = 4, fontface = "bold"
    )

  fn_p <- file.path(pr_dir, sprintf("analysis19_planned_dSVA_PCA_%s.png", st$id))
  ggsave(fn_p, plot = p_pca_s, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved PCA-adjusted: %s\n", fn_p))
}

cat("\n=== analysis19_planned_dSVA completed ===\n")
