# =============================================================================
# 33_ED9_assemble.R — Extended Data figure composer for
# ed_per_variant_pheno.
# -----------------------------------------------------------------------------
# v6.41 (2026-05-21): rewritten in Fig 1 v6 / ED04 / ED08 style —
# direct patchwork, native cell dims, no composer-side theme override
# (producers bake theme_nature_panel; forests ride on table_forest's
# theme_void). Composite shrunk from 16 × 9 in to Nature double-col
# (7.20 × 6.69 in).
#
# Four-panel composite, 2 rows × 2 cols (uniform 2x2 grid):
#
#   Row 1: A B   sm_ordinal       | sm_comp_variant
#   Row 2: C D   hr_variant       | clinical_variant
#
# v6.42 (2026-05-21): switched from a full-row-1 forest + 3 narrow
# heatmaps to a balanced 2×2 — heatmaps were horizontally squished at
# 2.20-in wide × 5 variant cols. 3.30-in wide × 5 cols gives each cell
# real breathing room while the forest still has enough horizontal room
# for its table columns thanks to the table-forest compact mode.
#
# Composite footprint: 7.20 × 6.69 in (Nature double-col × max height).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
.pa <- load_panel_assignments()

panel_root <- here("results", "ExtendedData", "ed_per_variant_pheno")
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)

load_panel <- function(token) {
  d <- panel_slot_dir(token, .pa)
  readRDS(file.path(d, paste0(token, ".rds")))
}

p_A <- load_panel("sm_ordinal")           # forest, theme_void from table_forest
p_B <- load_panel("sm_comp_variant")      # heatmap
p_C <- load_panel("hr_variant")           # heatmap
p_D <- load_panel("clinical_variant")     # heatmap

# ---- Compose ----------------------------------------------------------------
design <- "
AB
CD
"

composite <- p_A + p_B + p_C + p_D +
  plot_layout(design = design,
              widths  = c(1, 1),
              heights = c(3.04, 3.04)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d"))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family),
        # v6.41: tight inter-panel margins (ED08 v6.40 pattern) so the
        # composite holds within Nature's 6.69 in max height after
        # adding the 4 panel-tag rows.
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))

# ---- Save -------------------------------------------------------------------
.ed_num <- .pa$group_number$EDFig[["ed_per_variant_pheno"]]
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
