# =============================================================================
# 35_F3_assemble.R — Fig 3 (registry: rupture_score)
# -----------------------------------------------------------------------------
# 2026-05-23: Half-column (single-column, 89 mm / 3.50 in wide) portrait
# layout. A stacked above B so the KM curve renders at a natural square-ish
# aspect without horizontal stretching.
#
# Layout (2 panels, 2 rows × 1 col):
#
#     A
#     B
#
#   a km_by_score              KM by integer score     (top)
#   b rupture_lookup_heatmap   cumulative-risk grid     (bottom)
#
# Legends: kept inline — score-level colour legend on A, fill colourbar
# on B, both positioned at bottom to avoid wasted right-margin space in
# a narrow single-column format.
#
# Theme handling: theme_nature_panel() applied to both panels; legends
# are NOT stripped.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "panel_dims.R"))
.cfg <- FIGURE_LAYOUTS$rupture_score

load_panel <- function(token, file = NULL) {
  dir <- panel_slot_dir(token)
  path <- if (is.null(file)) file.path(dir, paste0(token, ".rds"))
          else               file.path(dir, file)
  readRDS(path)
}

.no_grid <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

.km_collapse_margin <- theme(
  axis.title.y      = element_text(margin = margin(r = 0, l = 0)),
  axis.title.y.left = element_text(margin = margin(r = 0), hjust = 0.5),
  plot.margin       = margin(2, 2, 2, 0, "pt")
)

p_A <- load_panel("km_by_score", "km_by_score__curve.rds") +
  theme_nature_panel() + .km_collapse_margin + .no_grid +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.margin    = margin(0, 0, 0, 0))

p_B <- load_panel("rupture_lookup_heatmap") +
  theme_nature_panel() + .no_grid +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin   = margin(0, 0, 0, 0)) +
  guides(fill = guide_colorbar(
    direction  = "horizontal",
    barwidth   = unit(2.0, "cm"),
    barheight  = unit(0.2, "cm"),
    title.position = "top"
  ))

# ---- Compose ----------------------------------------------------------------
composite <- p_A + p_B +
  plot_layout(ncol = 1, heights = c(1.2, 1.0)) +
  plot_annotation(
    tag_levels = list(c("a", "b"))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family))

# ---- Save -------------------------------------------------------------------
panel_root <- here("results",
                   paste0("Figure", figure_number("rupture_score")))
out_stem   <- file.path(panel_root,
                        sprintf("Fig%d_composite",
                                figure_number("rupture_score")))

.W <- .cfg$width_in   # 3.50 in (89 mm, single column)
.H <- .cfg$height_in  # 5.50 in (140 mm) — under Nature 170 mm cap

ggsave(paste0(out_stem, ".pdf"), composite,
       width = .W, height = .H, device = cairo_pdf, family = NM$font_family)
ggsave(paste0(out_stem, ".png"), composite,
       width = .W, height = .H, dpi = 300, type = "cairo")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_stem, ".svg"), composite,
         width = .W, height = .H, device = svglite::svglite, fix_text_size = FALSE)
}
cat(sprintf(
  "  ✓ composite: %s.{pdf,png}  (%.2f x %.2f in)\n",
  out_stem, .W, .H
))
