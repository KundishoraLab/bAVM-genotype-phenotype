# =============================================================================
# 31_ED7_assemble.R — composer for the merged
# ed_vaf_deep_dive Extended Data figure.
# -----------------------------------------------------------------------------
# v6.29 (2026-05-20): rewritten in Fig 1 v6 style — direct patchwork,
# shared variant-shape legend in a right-margin gutter. Drops the
# vaf_outlier_combined 2x2 composite (its info now lives in
# SuppTable10_vaf_age_sensitivity).
#
# Seven-panel composite, 2 cols × 4 rows (legend in the row-1 right cell).
# Layout per author 2026-05-27 — a clockwise rotation of the legend,
# vaf_sm_total (b) and vaf_highrisk (g) lands the legend top-right:
#
#   Row 1: A L   a vaf_age_rupture_scatter | legend
#   Row 2: B G   b vaf_sm_total            | g vaf_highrisk
#   Row 3: C D   c vaf_sm_size             | d vaf_drainage
#   Row 4: E F   e vaf_eloquence           | f vaf_rupture
#
# NB: patchwork fills design AREAS in alphabetical letter order from the
# plot ADD order, so plots are added p_A..p_G then L, and each area letter
# is placed to hold its matching panel (area X = p_X; area L = legend).
# Panel tags follow the add order. Panel A is the §3 mirror of Fig 1
# Panel G (VAF -> age dose-response). L stacks the variant-shape legend
# (KRAS G12D / G12V / Other KRAS / BRAF) + the heatmap colorbar.
#
# Composite footprint: 7.20 × 6.69 in (Nature double-col × max height).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork); library(cowplot)
  library(grid)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
.pa <- load_panel_assignments()

panel_root <- here("results", "ExtendedData", "ed_vaf_deep_dive")
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)

# ---- Load panels -----------------------------------------------------------
load_panel <- function(token) {
  d <- panel_slot_dir(token, .pa)
  readRDS(file.path(d, paste0(token, ".rds")))
}

.strip <- theme(legend.position = "none")
.no_grid <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

p_A <- load_panel("vaf_age_rupture_scatter") + .strip + .no_grid
p_B <- load_panel("vaf_sm_total")            + .strip + .no_grid
p_C <- load_panel("vaf_sm_size")             + .strip + .no_grid
p_D <- load_panel("vaf_drainage")            + .strip + .no_grid
p_E <- load_panel("vaf_eloquence")           + .strip + .no_grid
p_F <- load_panel("vaf_rupture")             + .strip + .no_grid
p_G <- load_panel("vaf_highrisk")            + .strip + .no_grid

# ---- Shared variant legend (right-margin gutter) ---------------------------
# Build a synthetic ggplot mapping the 4 variant categories to their
# SHAPE_VARIANT codes (KRAS G12D filled circle / G12V open circle /
# Other KRAS open triangle / BRAF X). cowplot::get_legend() lifts the
# shape key out for stacking in column L.
.variant_levels <- c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF")
.variant_shapes <- c(
  "KRAS G12D"  = 16,
  "KRAS G12V"  = 1,
  "Other KRAS" = 2,
  "BRAF"       = 4
)
.legend_df <- data.frame(
  Variant = factor(.variant_levels, levels = .variant_levels),
  x = seq_along(.variant_levels), y = 0
)
.legend_plot <- ggplot(.legend_df, aes(x = x, y = y, shape = Variant)) +
  geom_point(colour = "#1A1A1A", size = 1.6, stroke = 0.5) +
  scale_shape_manual(values = .variant_shapes, name = "Variant",
                     breaks = .variant_levels) +
  theme_nature_panel() +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.margin    = margin(0, 0, 0, 0)) +
  guides(shape = guide_legend(override.aes = list(size = 1.8), ncol = 1))
.shared_legend <- get_legend(.legend_plot)

# v6.31 (2026-05-20): heatmap colorbar (panels D, E, F median-VAF fill
# scale: low=white, high=PAL_HEAT_HIGH = W_BLUE) extracted from a
# synthetic ggplot so it can be stacked in the right-margin gutter
# alongside the variant legend. Per-panel colorbar suppressed inside
# each heatmap (theme(legend.position = "none")), so the gutter is the
# only place the gradient is keyed.
.cb_df <- data.frame(x = 1, y = 1, v = c(0, 5, 10))
.colorbar_plot <- ggplot(.cb_df, aes(x = x, y = y, fill = v)) +
  geom_tile() +
  # v6.32 (2026-05-20): annotate ONLY min/max (0 and 10) on the colorbar
  # — intermediate ticks were visually busy at the small gutter footprint.
  scale_fill_gradient(low = "white", high = PAL_HEAT_HIGH,
                      name = "Median\nVAF (%)",
                      limits = c(0, 10),
                      breaks = c(0, 10)) +
  theme_nature_panel() +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.margin    = margin(0, 0, 0, 0)) +
  guides(fill = guide_colorbar(barheight = unit(1.2, "cm"),
                               barwidth  = unit(0.30, "cm"),
                               title.position = "top"))
.colorbar_legend <- get_legend(.colorbar_plot)

# v6.32 (2026-05-20): legend stack runs horizontally (variant shapes |
# heatmap colorbar) so the row-1 col-2 gutter reads as a single legend
# row at the top of the figure, not a tall vertical column.
.legend_stack <- cowplot::plot_grid(
  .shared_legend, .colorbar_legend,
  nrow = 1, align = "h", axis = "b",
  rel_widths = c(1.4, 1.0)
)
L <- wrap_elements(full = .legend_stack)

# ---- Compose ----------------------------------------------------------------
# v6.30 (2026-05-20): legend moved to row-1 col-2 instead of a tall
# right-margin gutter. Panel A (vaf_age_rupture_scatter) no longer
# horizontally stretched across both cols; it occupies row-1 col-1 at
# the same footprint as panels B-G below, giving the grid a uniform
# 2-col × 4-row layout.
design <- "
AL
BG
CD
EF
"

composite <- p_A + p_B + p_C + p_D + p_E + p_F + p_G + L +
  plot_layout(design = design,
              widths  = c(1, 1),
              heights = c(1.5, 1.5, 1.5, 1.5)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", ""))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family))

# ---- Save -------------------------------------------------------------------
.ed_num <- .pa$group_number$EDFig[["ed_vaf_deep_dive"]]
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
