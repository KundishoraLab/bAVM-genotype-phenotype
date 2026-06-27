# =============================================================================
# check_stats_manifest.R — hard freshness + integrity validator
# -----------------------------------------------------------------------------
# Runs three independent checks and halts with an actionable error if any
# fails. Intended to be sourced at the top of every parent .qmd and by
# run_all.R immediately after build_stats_manifest.R.
#
#   (1) Timestamp check — each section's fragment mtime must be newer than
#       its producer script's mtime and newer than the upstream processed
#       data file.
#   (2) Schema check — every required key in stats_schema.R must be present
#       in the manifest.
#   (3) Content-hash check — each section's recorded producer SHA must match
#       the current SHA of the producer R script.
#
# Strictness is controlled by the env var BAVM_STATS_STRICT. Default is
# "warn" so day-to-day producer edits don't immediately fail the build
# (touching an R script bumps its mtime and would invalidate every
# stats fragment whose timestamp is older). Promote to "error" for CI
# and release builds via:
#   - the `--strict` flag on run_all.R (sets BAVM_STATS_STRICT=error
#     for the duration of that invocation), or
#   - BAVM_STATS_STRICT=error Rscript ... directly.
# BAVM_STATS_STRICT=off skips all checks entirely (drafting only).
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' required.", call. = FALSE)
  }
})

source(here::here("analysis", "pipeline", "helpers", "stats_schema.R"))

# null-coalesce helper (base R lacks one) — defined early so schema check uses it
`%||%` <- function(a, b) if (is.null(a)) b else a

STRICT <- tolower(Sys.getenv("BAVM_STATS_STRICT", unset = "warn"))
stopifnot(STRICT %in% c("error", "warn", "off"))

.flag <- function(msg) {
  if (STRICT == "error") stop(msg, call. = FALSE)
  if (STRICT == "warn")  warning(msg, call. = FALSE)
  # STRICT == "off": silent
}

# Early-exit wrapper: the rest of this file only runs if STRICT != "off".
# (Using if/else rather than return() because this file is sourced, not
# inside a function.)
.run_checks <- STRICT != "off"
if (!.run_checks) {
  message("check_stats_manifest: BAVM_STATS_STRICT=off \u2014 skipping all checks.")
}

if (.run_checks) {

# ---- load manifest ---------------------------------------------------------

manifest_path <- here::here("results", "stats", "manuscript_stats.rds")
meta_path     <- here::here("results", "stats", "manuscript_stats_meta.rds")

if (!file.exists(manifest_path)) {
  .flag(paste0(
    "STATS MANIFEST MISSING: ", manifest_path, "\n",
    "  Run: Rscript analysis/pipeline/builders/build_stats_manifest.R"
  ))
}

stats <- if (file.exists(manifest_path)) readRDS(manifest_path) else list()
meta  <- if (file.exists(meta_path))     readRDS(meta_path)     else list()

# ---- paths we depend on ----------------------------------------------------

upstream_data <- here::here("data", "processed", "bAVM_analysis_ready.rds")
upstream_clean <- here::here("analysis", "00_data_prep", "01_clean_master.R")

upstream_mtime <- max(
  if (file.exists(upstream_data))  file.info(upstream_data)$mtime  else as.POSIXct(NA),
  if (file.exists(upstream_clean)) file.info(upstream_clean)$mtime else as.POSIXct(NA),
  na.rm = TRUE
)

# ---- (1) timestamp check ---------------------------------------------------

ts_errors <- character(0L)
for (section in names(meta)) {
  m <- meta[[section]]
  frag_path <- here::here("results", "stats", "_manifest_fragments",
                          paste0(section, ".rds"))
  if (!file.exists(frag_path)) next

  frag_mtime <- file.info(frag_path)$mtime

  producer_full <- here::here(m$producer_path)
  if (file.exists(producer_full)) {
    prod_mtime <- file.info(producer_full)$mtime
    if (!is.na(prod_mtime) && frag_mtime < prod_mtime) {
      ts_errors <- c(ts_errors, sprintf(
        "  \u2716 [%s] fragment (%s) older than producer %s (%s).\n      Re-run: source('%s')",
        section, format(frag_mtime, "%Y-%m-%d %H:%M"),
        m$producer_path,
        format(prod_mtime, "%Y-%m-%d %H:%M"),
        m$producer_path
      ))
    }
  }

  if (!is.na(upstream_mtime) && frag_mtime < upstream_mtime) {
    ts_errors <- c(ts_errors, sprintf(
      "  \u2716 [%s] fragment older than upstream data (%s).\n      Re-run producer: source('%s')",
      section,
      format(upstream_mtime, "%Y-%m-%d %H:%M"),
      m$producer_path
    ))
  }
}

if (length(ts_errors) > 0L) {
  .flag(paste0("STATS STALE (timestamp):\n", paste(ts_errors, collapse = "\n")))
}

# ---- (2) schema check ------------------------------------------------------

schema_errors <- character(0L)
for (section in names(stats_schema)) {
  required <- stats_schema[[section]]
  if (!section %in% names(stats)) {
    schema_errors <- c(schema_errors, sprintf(
      "  \u2716 [%s] entire section missing. Producer: %s",
      section,
      stats_producers[[section]] %||% "<unknown>"
    ))
    next
  }
  present <- names(stats[[section]])
  missing <- setdiff(required, present)
  if (length(missing) > 0L) {
    schema_errors <- c(schema_errors, sprintf(
      "  \u2716 [%s] missing keys: %s",
      section, paste(missing, collapse = ", ")
    ))
  }
}

if (length(schema_errors) > 0L) {
  .flag(paste0("STATS SCHEMA MISMATCH:\n", paste(schema_errors, collapse = "\n")))
}

# ---- (3) content-hash check ------------------------------------------------

hash_errors <- character(0L)
for (section in names(meta)) {
  m <- meta[[section]]
  if (is.na(m$producer_sha)) next
  producer_full <- here::here(m$producer_path)
  if (!file.exists(producer_full)) {
    hash_errors <- c(hash_errors, sprintf(
      "  \u2716 [%s] producer script no longer exists: %s",
      section, m$producer_path
    ))
    next
  }
  current_sha <- digest::digest(file = producer_full, algo = "sha256")
  if (current_sha != m$producer_sha) {
    hash_errors <- c(hash_errors, sprintf(
      "  \u2716 [%s] producer SHA mismatch (script edited since last run).\n      Re-run: source('%s')",
      section, m$producer_path
    ))
  }
}

if (length(hash_errors) > 0L) {
  .flag(paste0("STATS HASH MISMATCH:\n", paste(hash_errors, collapse = "\n")))
}

# ---- success ---------------------------------------------------------------

if (length(ts_errors) == 0L &&
    length(schema_errors) == 0L &&
    length(hash_errors) == 0L) {
  message(sprintf(
    "\u2713 manifest checks passed (sections: %d, strictness: %s)",
    length(stats) - ("meta" %in% names(stats)),
    STRICT
  ))
}

}  # end if (.run_checks)
