# =============================================================================
# 32_ED8_assemble.R — composer for ed_anatomic_localization.
# -----------------------------------------------------------------------------
# v6.36 (2026-05-21): interaction_cleveland dropped from the figure —
# 27-cell facet grid couldn't render at Nature spec; the same data
# lives in SuppTable8 (location_interactions). Five-panel composite,
# with the three bottom panels (mosaic, parietal forest, meta forest)
# each getting their own full-width row.
#
# Five-panel composite, 4 rows × 2 cols (with D, E, F spanning full width):
#
#   Row 1: A B    per_lobe_prevalence  | per_variant_anatomy
#   Row 2: C C    rupture_categories    (full width)
#   Row 3: D D    parietal_forest       (full width)
#   Row 4: E E    rupture_meta_forest   (full width)
#
# Composite footprint: 7.20 × 6.69 in (Nature double-col × max height).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
.pa <- load_panel_assignments()

panel_root <- here("results", "ExtendedData", "ed_anatomic_localization")
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)

load_panel <- function(token) {
  d <- panel_slot_dir(token, .pa)
  readRDS(file.path(d, paste0(token, ".rds")))
}

# v6.38 (2026-05-21): theme_cell removed entirely. Every ED08 producer
# now bakes theme_nature_panel() (or theme_void via table_forest) in
# directly, so composer-side override is unnecessary AND harmful —
# theme_nature_panel uses %+replace% semantics that clobbered the
# producer's `theme(axis.text.x = element_text(angle = 45))` on Panel B
# and pulled axis ticks back onto the forest panels D + E.
p_A <- load_panel("per_lobe_prevalence")
p_B <- load_panel("per_variant_anatomy")
p_C <- load_panel("rupture_categories")
p_D <- load_panel("parietal_forest")

# ---- Compose ----------------------------------------------------------------
# v6.51 (2026-05-21): rupture_meta_forest dropped from the panel set
# (rupture prevalence no longer claimed in §3). 4 panels: row 1 = the
# two heatmaps side by side; rows 2 and 3 = mosaic + parietal forest
# each full-width. Reclaimed vertical room redistributed to the
# parietal forest (1.55 -> 2.20) which carries the headline §3
# parietal-by-rupture finding.
design <- "
AB
CC
DD
"

composite <- p_A + p_B + p_C + p_D +
  plot_layout(design = design,
              widths  = c(1, 1),
              heights = c(1.45, 1.45, 2.20)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d"))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family),
        # v6.40: tighten inter-panel whitespace. Default ggplot
        # plot.margin (~5.5 pt all sides) compounds across 4 rows and
        # eats vertical room that Panels C/D need for their content.
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))

# ---- Save -------------------------------------------------------------------
.ed_num <- .pa$group_number$EDFig[["ed_anatomic_localization"]]
out_stem <- file.path(panel_root, sprintf("ED%02d_composite", .ed_num))

.W <- 7.20
.H <- 6.69

ggsave(paste0(out_stem, ".pdf"), composite,
       width = .W, height = .H, device = cairo_pdf, family = NM$font_family)
ggsave(paste0(out_stem, ".png"), composite,
       width = .W, height = .H, dpi = 300, type = "cairo")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_stem, ".svg"), composite,
         width = .W, height = .H, device = svglite::svglite, fix_text_size = FALSE)
}

cat(sprintf("  ✓ ED%02d composite: %s.{pdf,png,svg}  (%.2f × %.2f in)\n",
            .ed_num, out_stem, .W, .H))
