# =============================================================================
# 29_ED6_assemble.R — composer for ed_km_diagnostics.
# -----------------------------------------------------------------------------
# v6.23 (2026-05-20): rewritten to use the Fig 1 v6 KM polish 1-to-1.
#
# Architecture mirrors 26_F1_assemble.R for Fig 1
# panels D / F / H:
#   * Load each KM panel's *__curve.rds* directly (skip the save_km_panel
#     pointer to drop the at-risk table baked into the assembly).
#   * Apply theme_nature_panel() + .km_collapse_margin + .strip + .no_grid
#     at compose time (the producer's theme_avm bake-in gets overridden).
#   * Extract one shared legend per palette via cowplot::get_legend on a
#     fresh RDS load (a stripped panel won't carry a guide grob).
#   * Stack both legends in a right-margin gutter (column L).
#
# Five-panel composite, 2 panel cols + 1 legend col × 3 rows:
#
#   Row 1: A B L   rare_variants  | km_pres_sex_F   | legend stack
#   Row 2: C D L   km_pres_sex_M  | km_rupt_sex_F   | (cont.)
#   Row 3: E . L   km_rupt_sex_M  | (empty)         | (cont.)
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

panel_root <- here("results", "ExtendedData", "ed_km_diagnostics")
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)

# ---- Load curves directly (no at-risk tables) -------------------------------
load_curve <- function(token) {
  d <- panel_slot_dir(token, .pa)
  readRDS(file.path(d, paste0(token, "__curve.rds")))
}

# ---- Compose-time theme stack (mirrors Fig 1 v6 KM panels D/F/H) -----------
.strip   <- theme(legend.position = "none")
.no_grid <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
.km_collapse_margin <- theme(
  axis.title.y      = element_text(margin = margin(r = 0, l = 0)),
  axis.title.y.left = element_text(margin = margin(r = 0), hjust = 0.5),
  plot.margin       = margin(2, 2, 2, 0, "pt")
)

theme_km_cell <- function(p) {
  p + theme_nature_panel() + .strip + .km_collapse_margin + .no_grid
}

p_A <- theme_km_cell(load_curve("rare_variants"))
p_B <- theme_km_cell(load_curve("km_pres_sex_F"))
p_C <- theme_km_cell(load_curve("km_pres_sex_M"))
p_D <- theme_km_cell(load_curve("km_rupt_sex_F"))
p_E <- theme_km_cell(load_curve("km_rupt_sex_M"))

# ---- Legend extraction (FRESH copies, not the stripped panels) -------------
# Two palettes share the gutter:
#   PAL_KM    — KRAS G12D / KRAS G12V / Panel-negative (4 sex KMs)
#   PAL_RARE  — Other KRAS / BRAF V600E / Panel-negative (rare_variants)
.legend_theme <- theme(
  legend.position  = "right",
  legend.direction = "vertical",
  legend.margin    = margin(0, 0, 0, 0)
)

.leg_genotype <- get_legend(
  load_curve("km_pres_sex_F") + theme_nature_panel() + .legend_theme
)
.leg_rare <- get_legend(
  load_curve("rare_variants") + theme_nature_panel() + .legend_theme
)

.legend_stack <- cowplot::plot_grid(
  .leg_genotype, .leg_rare,
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(1, 1)
)
L <- wrap_elements(full = .legend_stack)

# ---- Compose ----------------------------------------------------------------
# 2 panel cols + 1 legend col × 3 rows. Bottom-right panel cell ('#') is
# intentionally empty; the legend column L spans all three rows.
design <- "
ABL
CDL
E#L
"

composite <- p_A + p_B + p_C + p_D + p_E + L +
  plot_layout(design = design,
              widths  = c(1, 1, 0.45),
              heights = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d", "e", ""))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family))

# ---- Save -------------------------------------------------------------------
.ed_num <- .pa$group_number$EDFig[["ed_km_diagnostics"]]
out_stem <- file.path(panel_root, sprintf("ED%02d_composite", .ed_num))

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

cat(sprintf("  ✓ ED%02d composite: %s.{pdf,png,svg}  (%.2f × %.2f in)\n",
            .ed_num, out_stem, .W, .H))
