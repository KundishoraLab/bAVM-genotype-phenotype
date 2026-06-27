# =============================================================================
# 34_ED10_assemble.R — Extended Data figure composer for
# ed_null_audit.
# -----------------------------------------------------------------------------
# v6.44 (2026-05-21): rewritten in Fig 1 v6 / ED08 / ED09 style —
# direct patchwork, native cell dims, no composer-side theme_nature_panel
# override (producers bake it in; forests ride on table_forest's theme).
# Composite shrunk from 14 × 11 in to Nature double-col (7.20 × 6.69 in).
#
# Three-panel composite, 2 rows × 2 cols (C spans full bottom row):
#
#   Row 1: A B    mde_all (forest)        | age_adj_correlation (forest)
#   Row 2: C C    stenosis_ascertainment  (full-width stacked bar)
#
# Composite footprint: 7.20 × 6.69 in (Nature double-col × max height).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
.pa <- load_panel_assignments()

panel_root <- here("results", "ExtendedData", "ed_null_audit")
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)

load_panel <- function(token) {
  d <- panel_slot_dir(token, .pa)
  readRDS(file.path(d, paste0(token, ".rds")))
}

.strip <- theme(legend.position = "none")
.no_grid <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

p_A <- load_panel("mde_all")                # forest (theme_void via table_forest)
p_B <- load_panel("stenosis_ascertainment") # stacked bar (theme_nature_panel)
# v6.48 (2026-05-21): strip per-panel "Unadjusted / Age-adjusted" model
# legend on the two age-adj forests — the shared legend lives in the
# right-margin gutter (column L).
p_C <- load_panel("age_adj_binary")     + .strip
# patchwork::free() prevents D's longer outcome labels from forcing B
# (stacked bar, short study codes) to indent its plot area in column 2.
p_D <- patchwork::free(load_panel("age_adj_continuous") + .strip)

# v6.48: build the shared model legend from a fresh copy of panel C
# (PAL_SCORE teal + black for Unadjusted vs Age-adjusted).
.legend_model <- cowplot::get_legend(
  load_panel("age_adj_binary") +
    theme_nature_panel() +
    theme(legend.position  = "right",
          legend.direction = "vertical",
          legend.margin    = margin(0, 0, 0, 0))
)
L <- patchwork::wrap_elements(full = .legend_model)

# ---- Compose ----------------------------------------------------------------
# v6.48 (2026-05-21): 2x2 panel grid + right-margin shared-legend gutter
# (Fig 1 v6 pattern). L spans both rows. Widths c(1, 1, 0.35) keeps the
# panel cells at Fig 1 C/G size (~3.15 in) while the legend column is
# narrow but readable.
design <- "
ABL
CDL
"

composite <- p_A + p_B + p_C + p_D + L +
  plot_layout(design = design,
              widths  = c(1, 1, 0.35),
              heights = c(3.04, 3.04)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d", ""))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family),
        # v6.44: tight inter-panel margins (ED08 / ED09 pattern).
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))

# ---- Save -------------------------------------------------------------------
.ed_num <- .pa$group_number$EDFig[["ed_null_audit"]]
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
