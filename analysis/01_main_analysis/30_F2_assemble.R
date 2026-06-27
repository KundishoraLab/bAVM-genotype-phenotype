# =============================================================================
# 30_F2_assemble.R — Fig 2 (registry: gxp_associations)
# -----------------------------------------------------------------------------
# v6.6 (2026-05-23) — fully sequential a→e labeling achieved by:
#   (1) swapping prose citation order (parietal cited before hr_features_OR), and
#   (2) rearranging the patchwork cells to match visual reading order.
#
# Layout (3 rows × 5 cols, patchwork `design` string):
#
#     A A A B B
#     C C D D L
#     E E E E E
#
#   Row 1: A (cols 1-3)  SM grade stacked bar          [panel label a]
#          B (cols 4-5)  SM-components dumbbell        [panel label b]
#   Row 2: C (cols 1-2)  clinical-history dumbbell     [panel label c]
#          D (cols 3-4)  parietal-KRAS lollipop        [panel label d]
#          L (col 5)     legend stack (SM Grade on top, Genotype + N below)
#   Row 3: E (cols 1-5)  high-risk OR forest, full width [panel label e]
#
#   Panel letters follow first-citation order in prose (a=sm_grade_dist,
#   b=sm_components_dumbbell, c=clinical_history_dumbbell,
#   d=parietal_kras_rupture, e=hr_features_OR). Visual reading order a→e
#   is perfectly sequential left-to-right, top-to-bottom.
#     chain: p_A(sm_grade) + p_B(sm_comp) + p_C(clinical) + p_D(parietal) + p_E(hr) + L
#     tag:   "a",           "b",           "c",            "d",            "e",   ""
#
# Theme handling:
#   * A, B, C, D: theme_nature_panel() + theme(legend.position = "none").
#   * E (hr_features_OR forest): native theme_void() retained — see v6.2 note.
#
# Legend stack: SM Grade (top) + Genotype + N (bottom), rel_heights = c(1, 1.4).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork); library(cowplot)
})

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "panel_dims.R"))

load_panel <- function(token) {
  readRDS(file.path(panel_slot_dir(token), paste0(token, ".rds")))
}

# ---- Plot panels ------------------------------------------------------------
# A / B / C / D: apply theme_nature_panel() FIRST (complete theme), then
# strip legends so the override sticks. E: NO theme overlay — its
# producer-side theme_void() + annotate() layout is intentional.
# v6.4 (2026-05-21): blank horizontal gridlines on A/B/C/D AFTER
# theme_nature_panel() applies — the panel.grid.major line in
# theme_nature_panel() (via theme_classic %+replace% semantics)
# overrides the producer's blanking, so we re-blank here. All four
# panels have a categorical y-axis (KRAS-status strata for A,
# genotype groups for B, feature names for C/D), so per the canonical
# grid rule (FIGURE_FORMATTING_RULES.md §3.3) horizontal gridlines
# should not appear.
.strip_y_grid <- theme(panel.grid.major.y = element_blank(),
                       legend.position    = "none")

p_A <- load_panel("sm_grade_dist")             +
       theme_nature_panel() + .strip_y_grid
p_B <- load_panel("sm_components_dumbbell")    +
       theme_nature_panel() + .strip_y_grid
p_C <- load_panel("clinical_history_dumbbell") +
       theme_nature_panel() + .strip_y_grid +
       theme(axis.title.y = element_blank())   # match p_D which has y = NULL
p_D <- load_panel("parietal_kras_rupture")     +
       theme_nature_panel() + .strip_y_grid
p_E <- load_panel("hr_features_OR")            +
       theme(legend.position = "none")

# ---- Legends (extracted from FRESH copies that still show their guide) ------
.leg_smgrade <- get_legend(
  load_panel("sm_grade_dist") +
    theme_nature_panel() +
    theme(legend.position = "right",
          legend.margin   = margin(0, 0, 0, 0))
)
.leg_genotype <- get_legend(
  load_panel("sm_components_dumbbell") +
    theme_nature_panel() +
    theme(legend.position   = "right",
          legend.direction  = "vertical",
          legend.box        = "vertical",
          legend.margin     = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0)) +
    guides(colour = guide_legend(order = 1,
                                 override.aes = list(size = 3)),
           size   = guide_legend(order = 2))
)

# Stack: SM Grade on top, Genotype + N below. SM Grade is more compact
# (4 swatches) than Genotype + N together (2 colour + 3-break size), so
# the smaller rel_height for SM Grade keeps the proportions natural.
.legend_stack <- cowplot::plot_grid(
  .leg_smgrade, .leg_genotype,
  ncol = 1, align = "v", axis = "l",
  rel_heights = c(1, 1.4)
)
L <- wrap_elements(full = .legend_stack)

# ---- Compose using patchwork `design` so A/C/D lock to identical cells -----
# v6.3 bugfix: patchwork's plot_layout(design = ...) matches `+`-chain
# plots to design areas in ALPHABETICAL order of area letters (not the
# order in which letters first appear in the design string). v6.2's
# chain was A,B,C,D,L,E which made patchwork pair the 5th plot (the
# legend L) with the 5th alphabetical area (E, row 3 full-width) and
# the 6th plot (p_E forest) with the 6th area (L, row 2 col 5),
# crushing the forest into the tiny legend cell and stretching the
# legend across the bottom. Re-ordered to A,B,C,D,E,L so each plot
# lands in its intended area; tag_levels follows the same order.
design <- "
AAABB
CCDDL
EEEEE
"

composite <- p_A + p_B + p_C + p_D + p_E + L +
  plot_layout(design = design, heights = c(1.0, 1.0, 1.5)) +
  plot_annotation(
    # chain: sm_grade(A) + sm_comp(B) + clinical(C) + parietal(D) + hr(E) + L
    # sequential a→e; last "" blanks legend cell.
    tag_levels = list(c("a", "b", "c", "d", "e", ""))
  ) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family))

# ---- Save -------------------------------------------------------------------
panel_root <- here("results",
                   paste0("Figure", figure_number("gxp_associations")))
out_stem   <- file.path(panel_root,
                        sprintf("Fig%d_composite",
                                figure_number("gxp_associations")))

.W <- NM$width_in   # 7.20 in (183 mm)
.H <- 5.8           # 147 mm — under Nature's 170 mm cap

ggsave(paste0(out_stem, ".pdf"), composite,
       width = .W, height = .H, device = cairo_pdf, family = NM$font_family)
ggsave(paste0(out_stem, ".png"), composite,
       width = .W, height = .H, dpi = 300, type = "cairo")
# v6.49 (2026-05-21): emit SVG alongside PDF/PNG so Fig 2 is round-
# trippable into Adobe Illustrator with native text.
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_stem, ".svg"), composite,
         width = .W, height = .H, device = svglite::svglite, fix_text_size = FALSE)
}
cat(sprintf(
  "  ✓ composite: %s.{pdf,png,svg}  (%.2f x %.2f in)\n",
  out_stem, .W, .H
))
