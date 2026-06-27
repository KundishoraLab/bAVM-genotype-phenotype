# =============================================================================
# save_panel.R — canonical panel-save helper.
# -----------------------------------------------------------------------------
# Writes a ggplot/patchwork object as PDF + PNG (+ optionally an RDS cache
# for the composer pipeline) at exactly (w × h) inches. One source of truth
# replaces ~18 near-identical save_panel / save_trio definitions that used
# to live inline in every producer under analysis/01_main_analysis/.
#
# Why centralize:
#   * Deterministic-seed wrap around each ggsave() so geom_jitter() and any
#     other random-position layer lays dots in the SAME positions in the
#     PDF and PNG. Without this, the two outputs diverge — same data,
#     different dot offsets — and collaborators flipping between formats
#     see two visually different plots.
#   * Single edit point for any future change (resolution, cairo flag, RDS
#     persistence, etc.) instead of 18 copy-pasted bodies.
#
# Producer integration:
#   1. source() this file at the top of each producer:
#        source(here("analysis", "pipeline", "helpers", "save_panel.R"))
#   2. Use save_panel(dir, name, plot, w, h) directly when the producer
#      already passes a full directory path.
#   3. Wrap when the producer uses a subdir or out_root pattern, e.g.
#        save_panel <- function(subdir, name, plot, w, h)
#          save_panel_impl(file.path(out_root, subdir), name, plot, w, h,
#                          device = "cairo")
#      (See analysis/01_main_analysis/17_ed_scrna_phred_qc.R for an example.)
#
# Defaults match the most common producer (default ggsave device, write RDS).
# Cairo callers (most ED producers, the scRNA panels, the risk-score panels)
# pass device = "cairo" so the cairo_pdf + type="cairo" engine kicks in for
# correct unicode glyph rendering of β, ρ, ≤, ≥, etc.
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("save_panel: package 'ggplot2' is required. install.packages('ggplot2')",
         call. = FALSE)
  }
})

#' Save a ggplot/patchwork object as PDF + PNG (+ optionally RDS).
#'
#' @param dir       character. Directory to write into; created if missing.
#' @param name      character. File stem (no extension); resulting files are
#'                  `<dir>/<name>.pdf`, `<dir>/<name>.png`, and (if
#'                  `save_rds`) `<dir>/<name>.rds`. Used to derive the
#'                  per-panel RNG seed so jitter is stable across saves.
#' @param plot      ggplot/patchwork object.
#' @param w, h      numeric. Width and height in inches.
#' @param device    "default" (base ggsave defaults) or "cairo" (uses
#'                  cairo_pdf for the PDF and type="cairo" for the PNG so
#'                  unicode glyphs render correctly).
#' @param save_rds  logical. If TRUE, also write `<dir>/<name>.rds` with the
#'                  unrendered ggplot/patchwork object so the composer
#'                  pipeline (analysis/pipeline/helpers/compose_figure.R) can
#'                  re-theme the panel at native composite size without
#'                  re-fitting the analysis.
save_panel_impl <- function(dir, name, plot, w, h,
                            device   = c("default", "cairo"),
                            save_rds = TRUE) {
  device <- match.arg(device)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  png_path <- file.path(dir, paste0(name, ".png"))
  svg_path <- file.path(dir, paste0(name, ".svg"))

  # Deterministic seed: hash of `name` into the RNG range. Same panel name
  # always yields the same seed; different panels get different seeds so
  # their jitter patterns don't all visually align.
  .seed <- abs(sum(utf8ToInt(name))) %% .Machine$integer.max

  # RDS first: the composer pipeline reads from .rds, so we need it on
  # disk before any standalone PDF/PNG ggsave can crash on a font/device
  # issue and abort the producer. PDF/PNG failures stay non-fatal
  # (warning) instead of aborting the run.
  if (save_rds) {
    saveRDS(plot, file.path(dir, paste0(name, ".rds")))
  }

  .save_one <- function(path, ...) {
    tryCatch(
      ggplot2::ggsave(path, plot, width = w, height = h, ...),
      error = function(e)
        warning(sprintf("[save_panel] '%s' failed: %s - continuing (RDS saved).",
                        path, conditionMessage(e)),
                call. = FALSE)
    )
  }

  if (device == "cairo") {
    set.seed(.seed); .save_one(pdf_path, device = grDevices::cairo_pdf)
    set.seed(.seed); .save_one(png_path, dpi = 300, type = "cairo")
  } else {
    set.seed(.seed); .save_one(pdf_path)
    set.seed(.seed); .save_one(png_path, dpi = 300)
  }

  # SVG with fix_text_size=FALSE: text elements are <text> nodes in the SVG,
  # not width-locked via textLength/lengthAdjust, so Illustrator can select,
  # resize, and retype them directly.
  if (requireNamespace("svglite", quietly = TRUE)) {
    set.seed(.seed)
    .save_one(svg_path, device = svglite::svglite, fix_text_size = FALSE)
  }

  invisible(NULL)
}

# Default (non-cairo, with RDS) — matches the most common producer.
# Producers that need cairo or no-RDS pass the right flags; producers
# that use a subdir / out_root pattern wrap save_panel_impl in a thin
# local function (see header).
save_panel <- function(dir, name, plot, w, h, ...) {
  save_panel_impl(dir, name, plot, w, h, ...)
}

# Legacy alias used by 09_F1_km_age.R + 12_ED6_km_sex_stratified.R. Same
# semantics as save_panel; kept so those producers can adopt the
# centralized helper without changing their existing `save_trio()`
# call sites.
save_trio <- function(dir, name, plot, w, h, ...) {
  save_panel_impl(dir, name, plot, w, h, ...)
}
