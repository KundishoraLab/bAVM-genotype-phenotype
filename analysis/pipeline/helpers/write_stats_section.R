# =============================================================================
# write_stats_section.R — producer-side helper for the stats manifest
# -----------------------------------------------------------------------------
# Called at the end of every analysis/01_main_analysis/*.R script. Writes that
# script's contribution to results/stats/_manifest_fragments/<section>.rds
# along with provenance metadata (script path, SHA-256, timestamp).
#
# Usage at bottom of analysis/01_main_analysis/09_F1_km_age.R:
#
#   source(here::here("analysis", "pipeline", "helpers", "write_stats_section.R"))
#   write_stats_section(
#     section = "fig4a",
#     stats = list(
#       n                 = nrow(km1_df),
#       median_g12d       = unname(summary(fit1)$table["variant_group=KRAS G12D", "median"]),
#       median_g12v       = unname(summary(fit1)$table["variant_group=KRAS G12V", "median"]),
#       median_neg        = unname(summary(fit1)$table["variant_group=Negative", "median"]),
#       logrank_p         = sd_fit$chisq_log$pvalue,
#       pairwise_g12d_neg = pw_p["KRAS G12D", "Negative"]
#     )
#   )
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' required. install.packages('digest')", call. = FALSE)
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' required. install.packages('here')", call. = FALSE)
  }
})

#' Write one section's stats to a manifest fragment.
#'
#' @param section character. Section key, matching stats_schema.R entries
#'   (e.g., "cohort", "fig1", "fig4a", "fig5").
#' @param stats named list. Every key required by stats_schema$<section>
#'   should be present; missing keys will be flagged at render time by
#'   check_stats_manifest.R.
#' @param script_path character or NULL. Path to the producing R script.
#'   Auto-detected from the call stack if NULL (works from Rscript, source(),
#'   and most interactive sessions).
#'
#' @return invisibly, the path to the written fragment.
write_stats_section <- function(section, stats, script_path = NULL) {

  stopifnot(
    is.character(section), length(section) == 1L, nzchar(section),
    is.list(stats), length(stats) > 0L,
    !is.null(names(stats)), all(nzchar(names(stats)))
  )

  # ---- resolve script path -------------------------------------------------
  if (is.null(script_path)) {
    script_path <- .detect_caller_script()
  }

  # ---- compute script hash (SHA-256) --------------------------------------
  if (!is.null(script_path) && file.exists(script_path)) {
    script_sha <- digest::digest(file = script_path, algo = "sha256")
    rel_path <- .to_relative(script_path)
  } else {
    script_sha <- NA_character_
    rel_path   <- "<unknown>"
    warning("write_stats_section(): could not resolve producing script path; ",
            "hash check will be disabled for section '", section, "'.",
            call. = FALSE)
  }

  # ---- assemble fragment ---------------------------------------------------
  fragment <- list(
    section        = section,
    stats          = stats,
    producer_path  = rel_path,
    producer_sha   = script_sha,
    written_at     = Sys.time(),
    r_version      = paste(R.version$major, R.version$minor, sep = "."),
    hostname       = Sys.info()[["nodename"]]
  )

  # ---- write to disk -------------------------------------------------------
  frag_dir <- here::here("results", "stats", "_manifest_fragments")
  dir.create(frag_dir, showWarnings = FALSE, recursive = TRUE)
  frag_path <- file.path(frag_dir, paste0(section, ".rds"))
  saveRDS(fragment, frag_path)

  message(sprintf("  [stats] %-10s \u2190 %s (%d keys)",
                  section, basename(rel_path), length(stats)))
  invisible(frag_path)
}

# -----------------------------------------------------------------------------
# Internals
# -----------------------------------------------------------------------------

# Try, in order: Rscript commandArgs, sys.frames(), rstudioapi.
.detect_caller_script <- function() {

  # (1) Rscript
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("^--file=.+$", args))
  if (length(m) > 0L) {
    return(normalizePath(sub("^--file=", "", m[1]), mustWork = FALSE))
  }

  # (2) Call stack from source()
  frames <- sys.frames()
  for (fr in rev(frames)) {
    ofile <- tryCatch(get("ofile", envir = fr, inherits = FALSE),
                      error = function(e) NULL)
    if (!is.null(ofile) && is.character(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, mustWork = FALSE))
    }
  }

  # (3) RStudio active document
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      rstudioapi::hasFun("getSourceEditorContext")) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(normalizePath(ctx$path, mustWork = FALSE))
    }
  }

  NULL
}

.to_relative <- function(path) {
  root <- tryCatch(here::here(), error = function(e) getwd())
  path <- normalizePath(path, mustWork = FALSE)
  if (startsWith(path, root)) {
    rel <- substr(path, nchar(root) + 1L, nchar(path))
    sub("^/", "", rel)
  } else {
    path
  }
}
