# Planned ΔLL (planned_DLL): surgeon goal vs preoperative PI–LL mismatch
#   planned_DLL = LATpre_PI_LL - planned_PI_LL
# where planned_PI_LL is parsed from Patients_Attributes$demo_AlignGoal_PILL and
# nonnegative values in [0, PLANNED_PILL_RECODE_MAX] are recoded to 0.
#
# Planned SVA: parsed from Patients_Attributes$demo_AlignGoal_SVA — strip "cm",
# strip "<" and ">", read numeric; if value > 10, divide by 10 (e.g. mm → cm).

PLANNED_PILL_RECODE_MAX <- 0

# Planned ΔLL must be >= this value (°) in planned_DLL analyses.
PLANNED_DLL_MIN_KEEP <- -30

# Planned-ΔLL regression / PCA cohort: |preop PI–LL| above this threshold (°)
PREOP_ABS_PI_LL_GT <- 0

# Preop sagittal vertical axis C2–C7 (same units as COMBINE); must be strictly below this.
PREOP_SVA_C2_C7_LT <- 40

#' Keep rows with |LATpre_PI_LL| > PREOP_ABS_PI_LL_GT and
#' LATpre_SVA_C2_C7 < PREOP_SVA_C2_C7_LT (non-missing for both).
restrict_planned_dll_analysis_cohort <- function(df) {
  if (!"LATpre_PI_LL" %in% names(df)) {
    stop("restrict_planned_dll_analysis_cohort: LATpre_PI_LL required")
  }
  if (!"LATpre_SVA_C2_C7" %in% names(df)) {
    stop("restrict_planned_dll_analysis_cohort: LATpre_SVA_C2_C7 required")
  }
  dplyr::filter(
    df,
    !is.na(.data$LATpre_PI_LL),
    abs(.data$LATpre_PI_LL) > PREOP_ABS_PI_LL_GT,
    !is.na(.data$LATpre_SVA_C2_C7),
    .data$LATpre_SVA_C2_C7 < PREOP_SVA_C2_C7_LT
  )
}

parse_planned_pill <- function(x) {
  vapply(x, USE.NAMES = FALSE, FUN.VALUE = numeric(1), FUN = function(raw) {
    if (is.null(raw) || (length(raw) == 1L && is.na(raw))) {
      return(NA_real_)
    }
    if (is.numeric(raw)) {
      return(as.numeric(raw))
    }
    s <- trimws(as.character(raw))
    if (!nzchar(s)) {
      return(NA_real_)
    }
    s <- tolower(s)
    s <- gsub("\u00b0", "", s, fixed = TRUE)
    s <- trimws(gsub("degrees?|deg", "", s, ignore.case = TRUE))
    if (grepl("^<", s)) {
      rest <- sub("^<\\s*", "", s)
      val <- suppressWarnings(as.numeric(rest))
      if (!is.na(val)) {
        return(val)
      }
      return(10)
    }
    if (grepl("^>", s)) {
      rest <- sub("^>\\s*", "", s)
      return(suppressWarnings(as.numeric(rest)))
    }
    val <- suppressWarnings(as.numeric(s))
    if (!is.na(val)) {
      return(val)
    }
    m <- regexpr("-?[0-9]+[.]?[0-9]*", s, perl = TRUE)
    if (m > 0) {
      return(suppressWarnings(as.numeric(regmatches(s, m))))
    }
    NA_real_
  })
}

parse_planned_sva <- function(x) {
  scale_if_large <- function(v) {
    if (is.na(v)) {
      return(NA_real_)
    }
    if (v > 10) {
      return(v / 10)
    }
    v
  }
  vapply(x, USE.NAMES = FALSE, FUN.VALUE = numeric(1), FUN = function(raw) {
    if (is.null(raw) || (length(raw) == 1L && is.na(raw))) {
      return(NA_real_)
    }
    if (is.numeric(raw)) {
      return(scale_if_large(as.numeric(raw)))
    }
    s <- trimws(as.character(raw))
    if (!nzchar(s)) {
      return(NA_real_)
    }
    s <- tolower(s)
    s <- gsub("cm", "", s, fixed = TRUE)
    s <- gsub("[<>]", "", s)
    s <- trimws(s)
    val <- suppressWarnings(as.numeric(s))
    if (is.na(val)) {
      m <- regexpr("-?[0-9]+[.]?[0-9]*", s, perl = TRUE)
      if (m > 0) {
        val <- suppressWarnings(as.numeric(regmatches(s, m)))
      }
    }
    scale_if_large(val)
  })
}

#' Join parsed planned PI–LL and planned_DLL onto COMBINE rows.
#'
#' @param df Data frame from load_combine_data (project root working directory).
#' @param db_path Path to CADS Excel workbook.
#' @return df with demo_AlignGoal_PILL_raw, planned_PI_LL_parsed, planned_PI_LL,
#'   planned_minus_preop_PI_LL, planned_DLL added or updated; and when present in
#'   Patients_Attributes, demo_AlignGoal_SVA_raw and planned_SVA.
attach_planned_dll <- function(df, db_path) {
  pat_attr <- readxl::read_excel(
    db_path,
    sheet = "Patients_Attributes",
    col_types = "text",
    .name_repair = ~ make.unique(.x, sep = "__")
  )
  if (!"demo_AlignGoal_PILL" %in% names(pat_attr)) {
    stop("Patients_Attributes sheet must contain demo_AlignGoal_PILL")
  }
  if (!"demo_id" %in% names(pat_attr)) {
    stop("Patients_Attributes sheet must contain demo_id")
  }

  align <- pat_attr %>%
    dplyr::select(demo_id, demo_AlignGoal_PILL_raw = demo_AlignGoal_PILL)
  if ("demo_AlignGoal_SVA" %in% names(pat_attr)) {
    align <- align %>%
      dplyr::left_join(
        pat_attr %>%
          dplyr::select(demo_id, demo_AlignGoal_SVA_raw = demo_AlignGoal_SVA),
        by = "demo_id"
      )
  } else {
    align$demo_AlignGoal_SVA_raw <- NA_character_
  }
  align <- align %>%
    dplyr::mutate(
      planned_PI_LL = parse_planned_pill(.data$demo_AlignGoal_PILL_raw),
      planned_SVA = parse_planned_sva(.data$demo_AlignGoal_SVA_raw)
    )

  demo_id_col <- grep("^demo_id$", names(df), value = TRUE)[1]
  if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
    demo_id_col <- grep("^demo_id", names(df), value = TRUE)[1]
  }
  if (is.na(demo_id_col) || !nzchar(demo_id_col)) {
    demo_id_col <- "demo_id"
  }

  if (!"LATpre_PI_LL" %in% names(df)) {
    stop("COMBINE data must contain LATpre_PI_LL for planned_DLL")
  }

  df %>%
    dplyr::left_join(align, by = setNames("demo_id", demo_id_col)) %>%
    dplyr::mutate(
      planned_PI_LL_parsed = .data$planned_PI_LL,
      planned_PI_LL = dplyr::if_else(
        is.na(.data$planned_PI_LL_parsed),
        NA_real_,
        dplyr::if_else(
          .data$planned_PI_LL_parsed >= 0 &
            .data$planned_PI_LL_parsed <= PLANNED_PILL_RECODE_MAX,
          0,
          .data$planned_PI_LL_parsed
        )
      ),
      planned_minus_preop_PI_LL = .data$planned_PI_LL - .data$LATpre_PI_LL,
      planned_DLL = -.data$planned_minus_preop_PI_LL
    )
}

ensure_planned_results_dir <- function(subdir = NULL) {
  base <- "planned_results"
  if (!dir.exists(base)) {
    dir.create(base, recursive = TRUE)
  }
  if (!is.null(subdir)) {
    path <- file.path(base, subdir)
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    return(path)
  }
  base
}
