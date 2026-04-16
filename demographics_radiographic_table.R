#!/usr/bin/env Rscript
#
# Table: demographics (age, sex, BMI) and radiographic / planned / 6-week measures.
# Cohort: PJK-excluded COMBINE (same load as other scripts). Each row reports n
# for non-missing values for that variable (listwise n differs by field).
#
# Sex: COMBINE demo_gender is0 / 1 (see console note; typically 0 = female,
# 1 = male — verify against your codebook).
#
# Units: angles in degrees; SVA (C2–S1) in cm (database mm / 10).
# Planned ΔSVA (cm) = LATpre_SVA_C2_S1/10 − planned_SVA (same as analysis1_kf_correction_SVA.R).
# Knee flexion: LATpre_LL_KneeAngle (baseline), LAT6W_LL_KneeAngle (6 weeks).

suppressPackageStartupMessages({
  library(tidyverse)
})

options(warn = -1)

source("utils/utils.R")

EXCLUDE_PJK <- TRUE
db_path <- "/Users/ddliu/Desktop/ISSG/Retrospective_projects/Databases/CADS database - 2025.10.10.xlsx"

df <- load_combine_data(db_path, exclude_pjk = EXCLUDE_PJK)
df <- attach_planned_dll(df, db_path)

df$sv_pre_cm <- df$LATpre_SVA_C2_S1 / 10
df$sv_6w_cm <- df$LAT6W_SVA_C2_S1 / 10
if ("planned_SVA" %in% names(df) && "LATpre_SVA_C2_S1" %in% names(df)) {
  df$planned_dSVA_cm <- (df$LATpre_SVA_C2_S1 / 10) - df$planned_SVA
} else {
  df$planned_dSVA_cm <- NA_real_
}

# Planned ΔSVA summary: |Δ| cap as in analysis19_planned_dSVA.R (avoids extreme parsed goals)
MAX_PLANNED_DSVA_ABS_CM <- 100
ok_plan_sva <- !is.na(df$planned_dSVA_cm) & abs(df$planned_dSVA_cm) <= MAX_PLANNED_DSVA_ABS_CM

# Planned ΔLL summary: same outlier / floor screen as analysis19 (avoids extreme parsed goals)
MAX_PLANNED_DLL_ABS <- 100
planned_dll_for_summary <- df$planned_DLL
ok_plan <- !is.na(planned_dll_for_summary) &
  planned_dll_for_summary >= PLANNED_DLL_MIN_KEEP &
  abs(planned_dll_for_summary) <= MAX_PLANNED_DLL_ABS

age_col <- NULL
for (nm in c("demo_Age", "demo_age", "age", "Age")) {
  if (nm %in% names(df)) {
    age_col <- nm
    break
  }
}

summarize_continuous <- function(x, label, section) {
  x <- as.numeric(x)
  ok <- !is.na(x) & is.finite(x)
  n <- sum(ok)
  if (n < 1L) {
    return(tibble::tibble(
      section = section,
      parameter = label,
      n = n,
      mean = NA_real_,
      sd = NA_real_,
      value = "—"
    ))
  }
  m <- mean(x[ok])
  s <- stats::sd(x[ok])
  tibble::tibble(
    section = section,
    parameter = label,
    n = n,
    mean = m,
    sd = s,
    value = sprintf("%.2f (%.2f)", m, s)
  )
}

rows_demo <- list()

if (!is.null(age_col)) {
  rows_demo[[length(rows_demo) + 1L]] <- summarize_continuous(
    df[[age_col]], "Age (y)", "Demographics"
  )
} else {
  rows_demo[[length(rows_demo) + 1L]] <- tibble::tibble(
    section = "Demographics", parameter = "Age (y)", n = NA_integer_,
    mean = NA_real_, sd = NA_real_, value = "column missing"
  )
}

if ("BL_PE_BMI" %in% names(df)) {
  rows_demo[[length(rows_demo) + 1L]] <- summarize_continuous(
    df$BL_PE_BMI, "BMI (kg/m²)", "Demographics"
  )
} else {
  rows_demo[[length(rows_demo) + 1L]] <- tibble::tibble(
    section = "Demographics", parameter = "BMI (kg/m²)", n = NA_integer_,
    mean = NA_real_, sd = NA_real_, value = "column missing"
  )
}

if ("demo_gender" %in% names(df)) {
  g <- df$demo_gender
  ok <- !is.na(g)
  n_g <- sum(ok)
  if (n_g > 0L) {
    # COMBINE convention (verify in your codebook): 0 = female, 1 = male
    sex_label <- function(lev) {
      if (identical(as.numeric(lev), 0)) {
        return("Female, n (% of known sex)")
      }
      if (identical(as.numeric(lev), 1)) {
        return("Male, n (% of known sex)")
      }
      sprintf("Sex, code %s, n (%% of known sex)", lev)
    }
    for (lev in sort(unique(as.numeric(g[ok])))) {
      cnt <- sum(as.numeric(g[ok]) == lev, na.rm = TRUE)
      pct <- 100 * cnt / n_g
      lab <- sex_label(lev)
      rows_demo[[length(rows_demo) + 1L]] <- tibble::tibble(
        section = "Demographics",
        parameter = lab,
        n = cnt,
        mean = NA_real_,
        sd = NA_real_,
        value = sprintf("%d (%.1f%%)", cnt, pct)
      )
    }
  } else {
    rows_demo[[length(rows_demo) + 1L]] <- tibble::tibble(
      section = "Demographics", parameter = "Sex", n = 0L,
      mean = NA_real_, sd = NA_real_, value = "—"
    )
  }
} else {
  rows_demo[[length(rows_demo) + 1L]] <- tibble::tibble(
    section = "Demographics", parameter = "Sex", n = NA_integer_,
    mean = NA_real_, sd = NA_real_, value = "column missing"
  )
}

rows_radio <- list(
  summarize_continuous(df$LATpre_L1_S1, "Baseline lumbar lordosis L1–S1 (°)", "Radiographic"),
  summarize_continuous(df$LATpre_PI_LL, "Baseline PI–LL mismatch (°)", "Radiographic"),
  summarize_continuous(df$LATpre_S1PT, "Baseline pelvic tilt (°)", "Radiographic"),
  summarize_continuous(df$sv_pre_cm, "Baseline SVA C2–S1 (cm)", "Radiographic"),
  if ("LATpre_LL_KneeAngle" %in% names(df)) {
    summarize_continuous(
      df$LATpre_LL_KneeAngle,
      "Baseline knee flexion (LATpre_LL_KneeAngle, °)",
      "Radiographic"
    )
  } else {
    tibble::tibble(
      section = "Radiographic",
      parameter = "Baseline knee flexion (LATpre_LL_KneeAngle, °)",
      n = NA_integer_,
      mean = NA_real_,
      sd = NA_real_,
      value = "column missing"
    )
  },
  summarize_continuous(
    ifelse(ok_plan, planned_dll_for_summary, NA_real_),
    sprintf(
      "Planned change in LL — planned ΔLL (°); |Δ|≤%d°, ≥%.0f° (cf. analysis19)",
      MAX_PLANNED_DLL_ABS,
      PLANNED_DLL_MIN_KEEP
    ),
    "Radiographic"
  ),
  if ("planned_SVA" %in% names(df)) {
    summarize_continuous(
      ifelse(ok_plan_sva, df$planned_dSVA_cm, NA_real_),
      sprintf(
        "Planned change in SVA (cm); preop C2–S1/10 − planned_SVA; |Δ|≤%d cm (cf. analysis19_planned_dSVA)",
        MAX_PLANNED_DSVA_ABS_CM
      ),
      "Radiographic"
    )
  } else {
    tibble::tibble(
      section = "Radiographic",
      parameter = "Planned change in SVA (cm)",
      n = NA_integer_,
      mean = NA_real_,
      sd = NA_real_,
      value = "planned_SVA missing (attach_planned_dll)"
    )
  },
  summarize_continuous(df$LAT6W_L1_S1, "6-week lumbar lordosis L1–S1 (°)", "Radiographic"),
  summarize_continuous(df$LAT6W_PI_LL, "6-week PI–LL mismatch (°)", "Radiographic"),
  summarize_continuous(df$LAT6W_S1PT, "6-week pelvic tilt (°)", "Radiographic"),
  if ("LAT6W_LL_KneeAngle" %in% names(df)) {
    summarize_continuous(
      df$LAT6W_LL_KneeAngle,
      "6-week knee flexion (LAT6W_LL_KneeAngle, °)",
      "Radiographic"
    )
  } else {
    tibble::tibble(
      section = "Radiographic",
      parameter = "6-week knee flexion (LAT6W_LL_KneeAngle, °)",
      n = NA_integer_,
      mean = NA_real_,
      sd = NA_real_,
      value = "column missing"
    )
  },
  summarize_continuous(df$sv_6w_cm, "6-week SVA C2–S1 (cm)", "Radiographic")
)

out <- dplyr::bind_rows(!!!rows_demo, !!!rows_radio)

pr_dir <- ensure_planned_results_dir()
csv_path <- file.path(pr_dir, "demographics_radiographic_table.csv")
readr::write_csv(out, csv_path)

cat("\n=== Demographic & radiographic table ===\n")
cat(sprintf("Cohort: PJK excluded, N rows = %d\n", nrow(df)))
if ("demo_gender" %in% names(df)) {
  cat("demo_gender table (raw COMBINE codes):\n")
  print(table(df$demo_gender, useNA = "ifany"))
  cat("(Verify codebook: which code is female / male.)\n")
}
cat("\n")
print(as.data.frame(out[, c("section", "parameter", "n", "value")]), row.names = FALSE, right = FALSE)
cat(sprintf("\nWrote: %s\n", csv_path))
cat("\nFootnotes for manuscript:\n")
cat("  • Value column: mean (SD) for continuous variables; n (%) for sex categories.\n")
cat("  • n = non-missing count for that row’s variable.\n")
cat("  • SVA in cm (C2–S1; database mm scaled /10); planned ΔSVA = preop/10 − planned_SVA (cm);\n")
cat("    mean (SD) uses |planned ΔSVA| ≤ ", MAX_PLANNED_DSVA_ABS_CM, " cm.\n", sep = "")
cat("  • Knee flexion: LATpre_LL_KneeAngle (baseline), LAT6W_LL_KneeAngle (6 weeks), °.\n")
cat("  • planned ΔLL = LATpre_PI_LL − planned_PI_LL after attach_planned_dll();\n")
cat("    mean (SD) uses rows with planned ΔLL ≥ ", PLANNED_DLL_MIN_KEEP, "° and |planned ΔLL| ≤ ",
 MAX_PLANNED_DLL_ABS, "°.\n", sep = "")
