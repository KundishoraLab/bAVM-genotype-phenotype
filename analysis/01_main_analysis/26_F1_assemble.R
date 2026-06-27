# =============================================================================
# 26_F1_assemble.R — Fig 1 (registry: cohort_natural_history)
# -----------------------------------------------------------------------------
# v6.7 (2026-05-23) — CONSORT (A) spans rows 1-2 of col 1; the genotype
# + patients-n legend stack (L) occupies row 3 col 1; F and G each
# occupy a single cell in row 3 cols 2-3. No right-margin legend column.
#
# Layout (9 patches, 3 rows x 3 cols via patchwork `design`):
#
#     A B L
#     A C D
#     E F G
#
#   a CONSORT flow                external PNG, spans rows 1-2 of col 1
#   b variant_landscape           horizontal bar (row 1, col 2)
#   L legend stack                Genotype (top) + Patients-n (bottom),
#                                 row 1, col 3
#   c vaf_by_variant              violin + box (row 2, col 2)
#   d km_presentation             KM curve     (row 2, col 3)
#   e age_density_variant         density      (row 3, col 1)
#   f km_rupture                  KM curve     (row 3, col 2)
#   g vaf_age_scatter             scatter      (row 3, col 3)
#
# Per the v6 contract:
#   - Producers continue to save plots with theme_avm() for review-quality
#     previews; the composer's `theme_nature_panel() + theme(...)` chain
#     replaces it at compose time.
#   - All inline legends stripped via `theme(legend.position = "none")` so
#     the L cell is the only place legends render.
#   - Legends extracted from FRESH `load_panel()` copies (separate from
#     the stripped objects used in the grid) so cowplot::get_legend()
#     returns the rendered guide grob.
#   - Categorical-y axis (B) blanks panel.grid.major.y AFTER
#     theme_nature_panel() restores it.
#   - patchwork `design` string matches plots to areas alphabetically;
#     the `+`-chain order is A,B,C,D,E,F,G,L to match.
#   - tag_levels = list with "" for L so only a..g get panel tags.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork); library(cowplot)
  library(grid); library(png)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "panel_dims.R"))
.cfg <- FIGURE_LAYOUTS$cohort_natural_history

load_panel <- function(token, file = NULL) {
  dir <- panel_slot_dir(token)
  path <- if (is.null(file)) file.path(dir, paste0(token, ".rds"))
          else              file.path(dir, file)
  readRDS(path)
}

# ---- CONSORT flow → patchwork-compatible grob -------------------------------
.consort_png_path <- file.path(panel_slot_dir("consort_flow"),
                               "Fig1A_consort_flow.png")
if (!file.exists(.consort_png_path)) {
  .consort_png_path <- file.path(panel_slot_dir("consort_flow"),
                                 "consort_flow.png")
}
p_A <- wrap_elements(
  full = grid::rasterGrob(png::readPNG(.consort_png_path),
                          interpolate = TRUE)
)

# ---- ggplot panels (themed + legendless) ------------------------------------
# theme_nature_panel() applied first (complete theme via %+replace%),
# then per-panel overrides (grid blanking, legend strip). Order matters.
.strip <- theme(legend.position = "none")
.cat_y <- theme(panel.grid.major.y = element_blank())
# v6.6 (2026-05-20): blanket gridline removal applied to every panel.
# theme_nature_panel() restores grey93 major gridlines; KM panels (D, F,
# H) inherit additional grids from survminer. The user prefers a clean
# gridless look across the entire figure, so .no_grid lands after
# theme_nature_panel + .cat_y on every ggplot panel.
.no_grid <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

# v6.4 (2026-05-21): more aggressive collapse — r = 0 pulls the title
# flush against the y-axis line; axis.title.y.left specificity is
# needed because survminer bakes a non-default value into its $plot
# that axis.title.y alone may not reach. plot.margin l = 0 removes the
# left outer padding too. If the title ends up clipping at the cell
# edge after this lands, the fallback is plot.margin l = 4 plus a
# negative title r-margin (-2) to pull the title back toward the axis.
.km_collapse_margin <- theme(
  axis.title.y      = element_text(margin = margin(r = 0, l = 0)),
  axis.title.y.left = element_text(margin = margin(r = 0), hjust = 0.5),
  plot.margin       = margin(2, 2, 2, 0, "pt")
)

p_B <- patchwork::free(
  load_panel("variant_landscape") +
    theme_nature_panel() + .cat_y + .strip + .no_grid +
    # v6.2 (2026-05-21): B's y-axis tick labels are the mutation names
    # (longest = "Panel-negative", ~70 pt wide at 6 pt text). Give them
    # a small right-margin so they sit cleanly against the plot edge.
    # patchwork::free() prevents B's wide y-axis from forcing panels C
    # and F (same column) to indent their plot areas to match.
    theme(axis.text.y = element_text(margin = margin(r = 3)))
)
p_C <- load_panel("vaf_by_variant") +
  theme_nature_panel() + .strip + .no_grid
p_D <- load_panel("km_presentation", "km_presentation__curve.rds") +
  theme_nature_panel() + .strip + .km_collapse_margin + .no_grid
p_E <- load_panel("age_density_variant") +
  theme_nature_panel() + .strip + .no_grid
p_F <- load_panel("km_rupture", "km_rupture__curve.rds") +
  theme_nature_panel() + .strip + .km_collapse_margin + .no_grid
p_G <- load_panel("vaf_age_scatter") +
  theme_nature_panel() + .strip + .no_grid +
  # v6.4 (2026-05-21): removed the r = 8 margin from v6.2 that
  # artificially detached G's title from its y-axis. B/G column-2
  # alignment is handled by B's tick-label r-margin (above) instead.
  theme(axis.title.y = element_text(margin = margin(r = 0)))

# ---- Legend extraction from FRESH copies ------------------------------------
# Genotype (from D's clean copy — same 3-arm scale as F).
.leg_genotype <- get_legend(
  load_panel("km_presentation", "km_presentation__curve.rds") +
    theme_nature_panel() +
    theme(legend.position  = "right",
          legend.direction = "vertical",
          legend.margin    = margin(0, 0, 0, 0))
)

# v6.6 (2026-05-20): Patients-n size legend from B's clean copy.
# Producer exposes guide_legend on scale_size_area; strip color/fill
# guides so only the size key (open grey circles) appears.
.leg_bsize <- get_legend(
  load_panel("variant_landscape") +
    theme_nature_panel() +
    guides(colour = "none", fill = "none") +
    theme(legend.position  = "right",
          legend.direction = "vertical",
          legend.margin    = margin(0, 0, 0, 0))
)

# Legend stack L (Patients-n on top, Genotype below).
.legend_stack <- cowplot::plot_grid(
  .leg_bsize, .leg_genotype,
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(1.1, 1.3)
)
L <- wrap_elements(full = .legend_stack)

# ---- Compose ----------------------------------------------------------------
design <- "
ABL
ACD
EFG
"

composite <- p_A + p_B + p_C + p_D + p_E + p_F + p_G + L +
  plot_layout(design = design,
              widths  = c(1, 1, 1),
              heights = c(1.0, 1.0, 1.0)) +
  plot_annotation(
    tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", ""))
  ) &
  theme(plot.tag    = element_text(size = NM$label_pt, face = "bold",
                                   family = NM$font_family),
        plot.margin = margin(2, 8, 2, 8, unit = "pt"))

# ---- Save -------------------------------------------------------------------
panel_root <- here("results",
                   paste0("Figure", figure_number("cohort_natural_history")))
out_stem   <- file.path(panel_root,
                        sprintf("Fig%d_composite",
                                figure_number("cohort_natural_history")))

.W <- .cfg$width_in   # 7.20 in (183 mm)
.H <- .cfg$height_in  # 6.50 in (165 mm) — under Nature 170 mm cap

ggsave(paste0(out_stem, ".pdf"), composite,
       width = .W, height = .H, device = cairo_pdf, family = NM$font_family)
ggsave(paste0(out_stem, ".png"), composite,
       width = .W, height = .H, dpi = 300, type = "cairo")
# v6.49 (2026-05-21): emit SVG alongside PDF/PNG so Fig 1 is round-
# trippable into Adobe Illustrator with native text (same pattern
# applied across every ED composite).
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_stem, ".svg"), composite,
         width = .W, height = .H, device = svglite::svglite, fix_text_size = FALSE)
}
cat(sprintf(
  "  ✓ composite: %s.{pdf,png}  (%.2f x %.2f in)\n",
  out_stem, .W, .H
))
