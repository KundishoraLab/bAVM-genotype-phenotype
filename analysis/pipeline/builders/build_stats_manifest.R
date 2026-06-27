# =============================================================================
# build_stats_manifest.R — aggregate fragments → manuscript_stats.rds
# -----------------------------------------------------------------------------
# Reads every fragment in results/stats/_manifest_fragments/, assembles them
# into a single nested named list, and writes:
#   - results/stats/manuscript_stats.rds       (the manifest itself)
#   - results/stats/manuscript_stats_meta.rds  (per-section provenance)
#
# The Quarto .qmd files load the first; check_stats_manifest.R validates both.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
})

# ---- locate + read fragments -----------------------------------------------

frag_dir <- here::here("results", "stats", "_manifest_fragments")
if (!dir.exists(frag_dir)) {
  dir.create(frag_dir, showWarnings = FALSE, recursive = TRUE)
}

frag_files <- list.files(frag_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(frag_files) == 0L) {
  warning("build_stats_manifest(): no fragments found in ", frag_dir,
          ". Manifest will contain only meta.", call. = FALSE)
}

# ---- assemble manifest + metadata ------------------------------------------

stats <- list()
meta_per_section <- list()

for (f in frag_files) {
  frag <- readRDS(f)
  stopifnot(
    is.list(frag),
    all(c("section", "stats", "producer_path",
          "producer_sha", "written_at") %in% names(frag))
  )
  stats[[frag$section]] <- frag$stats
  meta_per_section[[frag$section]] <- list(
    producer_path = frag$producer_path,
    producer_sha  = frag$producer_sha,
    written_at    = frag$written_at,
    r_version     = frag$r_version,
    hostname      = frag$hostname
  )
}

# ---- git SHA (best effort) -------------------------------------------------

git_sha <- tryCatch(
  {
    sha <- system2("git", c("rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE)
    if (length(sha) == 0L || !nzchar(sha[1])) NA_character_ else sha[1]
  },
  error = function(e) NA_character_,
  warning = function(w) NA_character_
)

# ---- stamp meta block ------------------------------------------------------

stats$meta <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  git_sha      = git_sha,
  n_sections   = length(frag_files)
)

# ---- write outputs ---------------------------------------------------------

out_manifest <- here::here("results", "stats", "manuscript_stats.rds")
out_meta     <- here::here("results", "stats", "manuscript_stats_meta.rds")

saveRDS(stats, out_manifest)
saveRDS(meta_per_section, out_meta)

message(sprintf(
  "\u2713 manifest: %d sections \u2192 %s",
  length(frag_files),
  .to_relative_maybe <- {
    root <- here::here()
    if (startsWith(out_manifest, root)) {
      sub(paste0("^", root, "/?"), "", out_manifest)
    } else out_manifest
  }
))
