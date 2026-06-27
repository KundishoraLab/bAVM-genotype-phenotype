#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — regenerate every figure, supplementary table, and manuscript
# statistic from the cohort data.
#
#   Rscript analysis/pipeline/builders/run_all.R
#
# Steps:
#   (1) Reprocess the cohort data if the cleaner is newer than the cached
#       analysis-ready dataset.
#   (2) Run every producer in analysis/01_main_analysis/ in numeric order
#       (02 -> 37). Each writes its panels / tables and a stats fragment.
#   (3) Aggregate the fragments into results/stats/manuscript_stats.rds.
#   (4) Validate the stats manifest and panel-token uniqueness.
#   (5) Organize the manuscript-facing deliverable into results/figures/
#       (Figure N / Fig NA … + Extended Data Fig N + Tables).
#
# Figure numbers and panel letters are read from the frozen assignment cache
# results/stats/panel_assignments.rds (Fig 1-4, Extended Data 1-10).
#
# Flags:
#   --skip-data       skip step 1 (assume the processed data is current)
#   --skip-analysis   skip step 2 (assume the producer outputs are current)
#   --strict          treat stats-manifest schema warnings as errors
# =============================================================================

suppressPackageStartupMessages(library(here))

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(f) f %in% args
if (has_flag("--strict")) Sys.setenv(BAVM_STATS_STRICT = "error")

ROOT <- here::here()
setwd(ROOT)
message("=== bAVM genotype-phenotype: figure / table / stats regeneration ===")
message("root: ", ROOT)
message("")

# ---- (1) reprocess data -----------------------------------------------------
if (!has_flag("--skip-data")) {
  clean_script <- file.path(ROOT, "analysis", "00_data_prep", "01_clean_master.R")
  processed    <- file.path(ROOT, "data", "processed", "bAVM_analysis_ready.rds")
  if (!file.exists(clean_script)) {
    message("[1/5] no cleaner found; skipping data reprocess")
  } else if (!file.exists(processed) ||
             file.info(clean_script)$mtime > file.info(processed)$mtime) {
    message("[1/5] reprocessing raw data → ", basename(processed))
    source(clean_script, echo = FALSE)
  } else {
    message("[1/5] processed data up-to-date; skipping reprocess")
  }
} else {
  message("[1/5] --skip-data: skipping reprocess")
}

# ---- (2) run producers in numeric order ------------------------------------
if (!has_flag("--skip-analysis")) {
  scripts <- list.files(
    file.path(ROOT, "analysis", "01_main_analysis"),
    pattern = "^[0-9].*\\.R$", full.names = TRUE
  )
  scripts <- scripts[order(basename(scripts))]
  message(sprintf("[2/5] running %d producers (02 → 37) …", length(scripts)))
  for (s in scripts) {
    rel <- sub(paste0("^", ROOT, "/?"), "", s)
    message("  → ", rel)
    tryCatch(
      source(s, echo = FALSE),
      error = function(e)
        warning(sprintf("  [error in %s] %s", rel, conditionMessage(e)),
                call. = FALSE)
    )
  }
} else {
  message("[2/5] --skip-analysis: skipping producers")
}

# ---- (3) aggregate the stats manifest --------------------------------------
message("[3/5] building stats manifest …")
source(file.path(ROOT, "analysis", "pipeline", "builders", "build_stats_manifest.R"),
       echo = FALSE)

# ---- (4) validate ----------------------------------------------------------
message("[4/5] validating manifest + panel-token uniqueness …")
source(file.path(ROOT, "analysis", "pipeline", "validators", "check_stats_manifest.R"),
       echo = FALSE)
source(file.path(ROOT, "analysis", "pipeline", "validators", "check_panel_uniqueness.R"),
       echo = FALSE)
if (exists("check_panel_uniqueness")) check_panel_uniqueness()

# ---- (5) organize the manuscript-facing figure/table deliverable -----------
message("[5/5] organizing figures → results/figures/ …")
source(file.path(ROOT, "analysis", "pipeline", "builders", "organize_figures.R"),
       echo = FALSE)

message("")
message("=== done — figures, tables, and statistics under results/; clean deliverable in results/figures/ ===")
