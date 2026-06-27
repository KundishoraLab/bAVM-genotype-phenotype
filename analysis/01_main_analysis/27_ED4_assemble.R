# =============================================================================
# 27_ED4_assemble.R — Extended Data figure composer (new ED4).
# -----------------------------------------------------------------------------
# v6.21 (2026-05-20): 6-panel ed_cohort_heterogeneity split into TWO ED
# figures (this group + new ed_age_meta_variants in composer 45a).
# This group now carries the per-series heterogeneity story:
#
#   a — per_series_rate    JAMA stacked-bar of mutation rate per study
#                          (08_F1_cohort_variants.R)
#   b — per_series_vaf     per-study VAF distribution
#                          (08_F1_cohort_variants.R, violin + boxplot,
#                           Fig-1 polish 1-to-1: linewidth 0.4 / 0.3,
#                           fill alpha 0.2, jitter dots dropped)
#   c — pooled_rate_meta   REML pooled mutation-rate meta-analysis forest
#                          (11_ED4_power_forest_meta.R)
#
# Layout: 1 column × 3 rows. Each cell 6.60 × 2.03 in (inner); composer
# adds 0.3-in outer margin per side → final composite 7.20 × 6.69 in
# (Nature double-col x max height).
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("analysis", "pipeline", "helpers", "compose_figure.R"))
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa <- load_panel_assignments()

panel_root <- here("results", "ExtendedData", "ed_cohort_heterogeneity")

panel_rds <- list(
  per_series_rate   = file.path(panel_root, "per_series_rate.rds"),
  per_series_vaf    = file.path(panel_root, "per_series_vaf.rds"),
  pooled_rate_meta  = file.path(panel_root, "pooled_rate_meta.rds")
)

# 1-col x 3-row layout — panels read top-to-bottom in §1 citation order.
LAYOUT <- "A\nB\nC"

AREA_TO_TOKEN <- c(
  A = "per_series_rate",
  B = "per_series_vaf",
  C = "pooled_rate_meta"
)

# v6.25 (2026-05-20): rebalance heights so the pooled-rate forest
# (panel c) gets ~50% more vertical room. The forest has 8 study rows
# + 1 pooled row + axis + header — at uniform 2.03 in it was the most
# vertically cramped of the three panels. Panels A (stacked bar) and B
# (violin) tolerate a 20% squeeze.
# v6.41 (2026-05-28, Andy): widen inner cells 6.60 -> 6.90 in by trimming
# the composer's outer margin 0.30 -> 0.15 in. The forest in panel c uses
# table_forest_meta (single-row Rate "63.2 (50.6-74.3)") which is wider
# than the two-row layout that came before; the extra 0.30 in keeps the
# Rate string clear of the N column without overflowing Nature's 7.20-in
# (180-mm) double-column maximum (6.90 + 0.30 = 7.20).
ROW_HEIGHTS  <- c(1.62, 1.62, 2.85)   # 6.09 in inner total
COMPOSITE_W  <- 6.90
OUTER_MARGIN <- 0.15

compose_figure(
  panel_rds       = panel_rds,
  layout          = LAYOUT,
  area_to_token   = AREA_TO_TOKEN,
  out_stem        = file.path(panel_root,
                              sprintf("ED%02d_composite",
                                      .pa$group_number$EDFig[["ed_cohort_heterogeneity"]])),
  width_in        = COMPOSITE_W,
  height_in       = sum(ROW_HEIGHTS),
  panel_tag_size  = 8,
  heights         = ROW_HEIGHTS,
  outer_margin_in = OUTER_MARGIN
)
