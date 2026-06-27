# =============================================================================
# 03_prep_cohort_counts.R — cohort-level counts for the manuscript stats manifest
# -----------------------------------------------------------------------------
# Read the harmonised analysis-ready data and compute per-site and per-study
# counts. Writes the `cohort` section of the manuscript stats manifest.
#
# Phase 1 seed: minimal keys (total n + per-site + per-published-cohort).
# Extend with more cohort-level stats in later phases as the manuscript
# references them.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

# ---- load data -------------------------------------------------------------

data_path <- here::here("data", "processed", "bAVM_analysis_ready.rds")
if (!file.exists(data_path)) {
  stop("cohort_counts: missing ", data_path, call. = FALSE)
}

df <- readRDS(data_path)

# ---- resolve study / site columns (defensive) ------------------------------
# Prefer `study_clean` (canonical post-cleaning column) over the raw
# `study`. Audit 01 found these two columns disagree on 4 rows and the
# downstream fragments (cohort vs fig1) used different columns,
# producing per-series count drift (e.g. BCH 109 vs 110, UAB 71 vs 73).

study_col <- intersect(c("study_clean", "study", "cohort", "source"),
                       names(df))[1]
if (is.na(study_col)) {
  stop("cohort_counts: could not find a 'study' column in ",
       basename(data_path), call. = FALSE)
}

study_vals <- as.character(df[[study_col]])
study_vals[is.na(study_vals)] <- "UNKNOWN"

# Case-insensitive label matching so future renames of study factors don't
# silently break these counts.
.count_like <- function(pattern) {
  sum(grepl(pattern, study_vals, ignore.case = TRUE))
}

# ---- tissue-tested subset --------------------------------------------------

tissue_col <- intersect(
  c("tissue_tested", "geno_variant", "genotype_tested"),
  names(df)
)[1]

n_tissue_tested <- if (!is.na(tissue_col) && tissue_col == "geno_variant") {
  sum(!is.na(df[["geno_variant"]]))
} else if (!is.na(tissue_col)) {
  sum(as.logical(df[[tissue_col]]), na.rm = TRUE)
} else {
  NA_integer_
}

# ---- exclusion bookkeeping (from 01_clean_master.R) -----------------------
# Pre-exclusion total / excluded count / per-reason breakdown live in
# data/processed/cohort_exclusions.rds, written by the cleaning script.
# The Fig 1A CONSORT producer reads them via the manifest to render the
# pre-exclusion → exclusion arm → post-exclusion tier without hardcoded
# numbers. Falls back to NA if the cleaning script hasn't run since this
# file was added.
.exclusion_path <- here::here("data", "processed", "cohort_exclusions.rds")
if (file.exists(.exclusion_path)) {
  excl <- readRDS(.exclusion_path)
  n_pre_exclusion <- as.integer(excl$n_pre_exclusion)
  n_excluded      <- as.integer(excl$n_excluded)
  excluded_reasons <- excl$reasons  # named list: reason -> count
} else {
  n_pre_exclusion  <- NA_integer_
  n_excluded       <- NA_integer_
  excluded_reasons <- list()
}

# ---- assemble stats list ---------------------------------------------------

cohort_stats <- list(
  n_pre_exclusion    = n_pre_exclusion,
  n_excluded         = n_excluded,
  excluded_reasons   = excluded_reasons,
  n_total_harmonised = nrow(df),
  n_tissue_tested    = as.integer(n_tissue_tested),
  n_bch              = .count_like("^bch$"),
  n_chop             = .count_like("^chop$"),
  n_uab              = .count_like("^uab$"),
  n_nikolaev         = .count_like("nikolaev"),
  n_priemer          = .count_like("priemer"),
  n_hong             = .count_like("hong"),
  n_goss             = .count_like("goss"),
  n_gao              = .count_like("gao")
)

# ---- write fragment --------------------------------------------------------

source(here::here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(section = "cohort", stats = cohort_stats)
