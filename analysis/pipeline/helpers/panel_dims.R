# =============================================================================
# panel_dims.R — per-panel native-scale dimension resolver.
# -----------------------------------------------------------------------------
# Each panel saved by a producer should arrive at its composite cell at the
# EXACT footprint and font size it will appear in the final Nature-formatted
# composite — no post-hoc theme override, no composer-side scaling. This file
# defines the per-figure layout descriptors and the lookup helpers producers
# use to render at native dimensions.
#
# Anchored on the Nature display-item spec (codified in
# canonical figure formatting):
#
#   single column width = 89 mm  = NATIVE$widths_in$single   = 3.50 in
#   double column width = 183 mm = NATIVE$widths_in$double   = 7.20 in
#   maximum height      = 170 mm = NATIVE$max_height_in      = 6.69 in
#   axis-text font size = 5-7 pt at print
#
# Workflow:
#   1. Composer's layout (LAYOUT string + AREA_TO_TOKEN + relative widths
#      / heights + composite width_in / height_in) is registered here in
#      FIGURE_LAYOUTS.
#   2. Producer calls save_panel_native(token, plot) which:
#        a. Resolves the token's group via the panel-assignments cache.
#        b. Looks up FIGURE_LAYOUTS[[group]] to find the token's row/col.
#        c. Computes cell dimensions in inches.
#        d. Applies theme_avm_native(base_size = NATIVE$base_size).
#        e. Calls save_panel_impl() at the resolved (w, h).
#   3. Composer reads the panel RDS and assembles via compose_figure() with
#      native_theme = FALSE — the panels arrive pre-themed at their target
#      cell dimensions.
#
# Adding a new figure: append a list entry to FIGURE_LAYOUTS keyed by group
# name. Required fields: width_in, height_in, layout (patchwork design
# string), area_to_token (named char vector), widths (relative col weights),
# heights (relative row weights).
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here::here("analysis", "helper_scripts", "utils.R"))
source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here::here("analysis", "pipeline", "helpers", "save_panel.R"))

# Convenience: Nature's max height limit, exposed via NATIVE for symmetry
# with widths_in. Producers / composers MUST keep composite total height
# <= this value. Hard cap: 170 mm = 6.69 in.
if (is.null(NATIVE$max_height_in)) {
  NATIVE$max_height_in <- 6.69
}

# ── Per-figure layout descriptors ──────────────────────────────────────────
# Each entry mirrors the composer's existing LAYOUT / AREA_TO_TOKEN /
# widths / heights / composite_w_in / composite_h_in but at the journal
# native scale (width_in = 7.20 in for double-column).
FIGURE_LAYOUTS <- list(
  # Fig 2 — gxp_associations
  # v6.6 (2026-05-23): sequential a→e layout matching visual reading order.
  #   Row 1: A (cols 1-3)  SM grade stacked bar         [panel a]
  #          B (cols 4-5)  SM-components dumbbell       [panel b]
  #   Row 2: C (cols 1-2)  clinical-history dumbbell    [panel c]
  #          D (cols 3-4)  parietal-KRAS lollipop       [panel d]
  #          L (col 5)     legend stack (not a data panel — excluded)
  #   Row 3: E (cols 1-5)  high-risk OR forest          [panel e]
  #
  # 5-col grid; C and D are both 2/5-wide and share the same row-height
  # weight → native cell 2.88 × 1.657 in each. Panels c and d therefore
  # pre-render at identical dimensions.
  # Composite: 7.20 × 5.80 in (147 mm height, under Nature 170 mm cap).
  gxp_associations = list(
    width_in      = NATIVE$widths_in$double,
    height_in     = 5.80,
    layout        = "AAABB\nCCDDL\nEEEEE",
    area_to_token = c(
      A = "sm_grade_dist",
      B = "sm_components_dumbbell",
      C = "clinical_history_dumbbell",
      D = "parietal_kras_rupture",
      E = "hr_features_OR"
      # L = legend; not a data panel, excluded from area_to_token
    ),
    widths  = c(1, 1, 1, 1, 1),
    heights = c(1.0, 1.0, 1.5)
  ),

  # Fig 1 — cohort_natural_history (7 panels, native scale)
  # 2026-05-23 v6.7: CONSORT flow (A) spans rows 1-2 of col 1. Row 3
  # col 1 is a prose-legend gutter (Z = dummy, not a data panel).
  # F (km_rupture) and G (vaf_age_scatter) each occupy a single cell in
  # row 3 cols 2-3. Total height 6.50 in (165 mm).
  #
  #   Row 1:  A | B | C    consort | variant_landscape | vaf_by_variant
  #   Row 2:  A | D | E    (A continues) | km_presentation | age_density
  #   Row 3:  Z | F | G    prose gutter (not a panel) | km_rupture | vaf_age_scatter
  cohort_natural_history = list(
    width_in      = NATIVE$widths_in$double,
    height_in     = 6.50,
    layout        = "ABZ\nACD\nEFG",
    area_to_token = c(
      A = "consort_flow",
      B = "variant_landscape",
      C = "vaf_by_variant",
      D = "km_presentation",
      E = "age_density_variant",
      F = "km_rupture",
      G = "vaf_age_scatter"
    ),
    widths  = c(1, 1, 1),
    heights = c(0.9, 1.0, 1.0)
  ),

  # Fig 3 — rupture_score (2 panels, single-column half-width)
  # 2026-05-23: half-column (89 mm) portrait format. km_by_score (A)
  # stacked above rupture_lookup_heatmap (B). KM curve gets a natural
  # aspect ratio at single-column width; heatmap sits below.
  rupture_score = list(
    width_in      = NATIVE$widths_in$single,
    height_in     = 5.50,
    layout        = "A\nB",
    area_to_token = c(
      A = "km_by_score",
      B = "rupture_lookup_heatmap"
    ),
    widths  = c(1),
    heights = c(1.2, 1.0)
  )
  # Other groups (ed_*) to be added in subsequent
  # migration passes. Until a group is registered here, its producers may
  # continue using save_panel(slot_dir(token), ...) and the composer's
  # native_theme = TRUE override path.
)

# Find which figure group contains a panel token.
.find_group_for_token <- function(token) {
  for (g in names(FIGURE_LAYOUTS)) {
    cfg <- FIGURE_LAYOUTS[[g]]
    if (token %in% unname(cfg$area_to_token)) return(g)
  }
  NA_character_
}

# Walk a layout string to find the (row, col) span occupied by `area`.
# Returns a list with row_idx (integer vector), col_idx (integer vector).
.layout_span <- function(layout_str, area) {
  rows <- strsplit(layout_str, "\n", fixed = TRUE)[[1]]
  span_rows <- c(); span_cols <- c()
  for (i in seq_along(rows)) {
    chars <- strsplit(rows[i], "", fixed = TRUE)[[1]]
    hits  <- which(chars == area)
    if (length(hits) > 0L) {
      span_rows <- c(span_rows, i)
      span_cols <- c(span_cols, hits)
    }
  }
  list(row_idx = unique(span_rows), col_idx = unique(span_cols))
}

#' Resolve the native-scale cell dimensions for a panel token.
#'
#' Returns NULL if the token's group has no FIGURE_LAYOUTS entry yet
#' (caller falls back to the legacy path).
#'
#' @param token character; the panel token (e.g., "stenosis_waffle").
#' @return list with `w_in` and `h_in`, or NULL if unregistered.
panel_native_dims <- function(token) {
  g <- .find_group_for_token(token)
  if (is.na(g)) return(NULL)
  cfg <- FIGURE_LAYOUTS[[g]]
  # Find which area letter maps to this token.
  area_idx <- match(token, unname(cfg$area_to_token))
  if (is.na(area_idx)) return(NULL)
  area <- names(cfg$area_to_token)[area_idx]
  span <- .layout_span(cfg$layout, area)
  # Cell width = sum of relative col weights for the cols this area
  # spans, divided by the total weight, times composite width_in.
  w_in <- cfg$width_in *
          sum(cfg$widths[span$col_idx]) / sum(cfg$widths)
  h_in <- cfg$height_in *
          sum(cfg$heights[span$row_idx]) / sum(cfg$heights)
  list(w_in = w_in, h_in = h_in, group = g)
}

#' Save a panel at its native-scale composite cell dimensions with
#' theme_avm_native() pre-applied. Producers use this instead of
#' save_panel() so the standalone PNG/PDF/RDS triple matches what
#' the composer will embed at print.
#'
#' If the token has no FIGURE_LAYOUTS entry, falls back to the legacy
#' save_panel() with a console warning so unmigrated panels still write.
#'
#' @param token character; panel token (must be registered in
#'   panel_registry.R and in FIGURE_LAYOUTS for native sizing).
#' @param plot ggplot or patchwork object.
#' @param ... additional args forwarded to save_panel_impl() (e.g.,
#'   device = "cairo").
save_panel_native <- function(token, plot, ...) {
  d <- panel_native_dims(token)
  if (is.null(d)) {
    warning(sprintf(
      "[save_panel_native] token '%s' has no FIGURE_LAYOUTS entry yet; ",
      "falling back to legacy save_panel(). Migrate the panel's figure ",
      "to FIGURE_LAYOUTS to enable native sizing.", token
    ), call. = FALSE)
    return(invisible(NULL))
  }
  dir <- slot_dir(token)
  if (is.null(dir)) {
    stop(sprintf(
      "[save_panel_native] slot_dir('%s') is NULL; is the token ",
      "registered in panel_registry.R?", token
    ), call. = FALSE)
  }
  # NB: the producer is responsible for applying theme_avm_native()
  # before any custom overrides (e.g. axis.text = element_blank() on
  # a waffle, theme_void() on a forest). save_panel_native() does NOT
  # re-theme — adding theme_avm_native() here at the end would
  # clobber the producer's plot-specific theme overrides (ggplot
  # themes layer with the LAST addition winning). The canvas size
  # alone (w/h from panel_native_dims) is what makes the saved file
  # native scale; the producer must pair that with theme_avm_native()
  # so the in-canvas text is the print text size.
  save_panel_impl(dir, token, plot,
                  w = d$w_in, h = d$h_in,
                  device = "cairo", ...)
}

# slot_dir alias — provider scripts source utils.R which exports
# panel_slot_dir; some legacy producers define a local `slot_dir` wrapper.
# Provide one here for consistency so save_panel_native() works in any
# producer scope.
if (!exists("slot_dir", mode = "function")) {
  slot_dir <- function(token) {
    pa <- tryCatch(load_panel_assignments(), error = function(e) NULL)
    panel_slot_dir(token, pa)
  }
}
