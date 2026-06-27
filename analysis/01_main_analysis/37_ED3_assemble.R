# =============================================================================
# 37_ED3_assemble.R — Extended Data figure composer (ED03).
# -----------------------------------------------------------------------------
# Six per-variant dPCR amplitude-cluster scatter panels assembled from the
# RDS files written by analysis/01_main_analysis/23_ED3_dpcr_validation.R.
#
# v6.20 (2026-05-20): rebuilt to follow the Fig 1 v6 architecture —
# patchwork directly, with one shared colour legend extracted in the
# composer and placed in the right-margin gutter. Producer no longer
# carries an inline colour legend (all six panels share the background
# and WT-clamp-suppressed greys; the mutant colour is variant-specific
# and shown in the gutter beside the panel grid).
#
# Layout (2 panel cols + 1 legend col, 3 rows):
#     A B L
#     C D L
#     E F L
#   where A..F = dpcr_scatter_{g12d, g12v, g12c, g12a, v600e, tert} in
#   registry order and L = shared legend stack.
#
# Composite footprint: 7.20 x 6.69 in (Nature double-col x max height).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork); library(cowplot)
  library(grid); library(ggtext)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "pipeline", "helpers", "panel_dims.R"))
source(here("analysis", "helper_scripts", "utils.R"))

.pa <- load_panel_assignments()
panel_root <- here("results", "ExtendedData", "ed_dpcr_validation")

load_panel <- function(key) {
  readRDS(file.path(panel_root, key, paste0(key, ".rds")))
}

# ---- Load the six themed/legend-stripped panels ------------------------------
p_A <- load_panel("dpcr_scatter_g12d")
p_B <- load_panel("dpcr_scatter_g12v")
p_C <- load_panel("dpcr_scatter_g12c")
p_D <- load_panel("dpcr_scatter_g12a")
p_E <- load_panel("dpcr_scatter_v600e")
p_F <- load_panel("dpcr_scatter_tert")

# ---- Build the shared colour legend ------------------------------------------
# Seven categories: 6 dye-channel positive-partition colours (one per panel,
# keyed by fluorophore emission spectrum — PAL_FLUOROPHORE) + a single grey for
# negative (below-threshold) partitions, shared across all panels. Labels lead
# with the dye because colour now encodes the emission channel, not the variant
# tier. Cowplot `get_legend()` is run on a tiny synthetic ggplot whose only
# purpose is to carry the scale; nothing about it is rendered as data.
.legend_categories <- c(
  "FAM · KRAS G12D"              = PAL_FLUOROPHORE[["FAM"]],
  "HEX · KRAS G12V"              = PAL_FLUOROPHORE[["HEX"]],
  "ROX · KRAS G12C"              = PAL_FLUOROPHORE[["ROX"]],
  "Cy5 · KRAS G12A"              = PAL_FLUOROPHORE[["Cy5"]],
  "TAMRA · BRAF V600E"           = PAL_FLUOROPHORE[["TAMRA"]],
  "Cy5.5 · TERT (pos. control)"  = PAL_FLUOROPHORE[["Cy5.5"]],
  "Negative (below threshold)"   = "#D0D0D0"
)
.legend_df <- data.frame(
  category = factor(names(.legend_categories), levels = names(.legend_categories)),
  x = 0, y = seq_along(.legend_categories)
)
.legend_plot <- ggplot(.legend_df, aes(x = x, y = y, colour = category)) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = .legend_categories, name = NULL,
                      breaks = names(.legend_categories)) +
  theme_nature_panel() +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.key.size  = unit(0.32, "cm"),
        legend.margin    = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1.8),
                               ncol = 1))
.shared_legend <- get_legend(.legend_plot)
L <- wrap_elements(full = .shared_legend)

# ---- Compose ----------------------------------------------------------------
design <- "
ABL
CDL
EFL
"

composite <- p_A + p_B + p_C + p_D + p_E + p_F + L +
  plot_layout(design = design,
              widths  = c(1, 1, 0.55),
              heights = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d", "e", "f", ""))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family))

# ---- Save -------------------------------------------------------------------
out_stem <- file.path(panel_root,
                      sprintf("ED%02d_composite",
                              .pa$group_number$EDFig[["ed_dpcr_validation"]]))

.W <- 7.20   # Nature double-column width
.H <- 6.69   # Nature max height

ggsave(paste0(out_stem, ".pdf"), composite,
       width = .W, height = .H, device = cairo_pdf, family = NM$font_family)
ggsave(paste0(out_stem, ".png"), composite,
       width = .W, height = .H, dpi = 300, type = "cairo")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_stem, ".svg"), composite,
         width = .W, height = .H, device = svglite::svglite, fix_text_size = FALSE)
}

cat(sprintf("  ✓ ED03 composite: %s.{pdf,png,svg}  (%.2f × %.2f in)\n",
            out_stem, .W, .H))
