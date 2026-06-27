# =============================================================================
# 28_ED5_assemble.R — Extended Data figure composer (new ED5).
# -----------------------------------------------------------------------------
# v6.21 (2026-05-20): new ED group spun out of ed_cohort_heterogeneity to
# carry the "outcome signal converges across cohorts" narrative paired
# with the rare-variant KM.
#
#   a — rupture_meta       REML pooled rupture-rate meta-analysis forest
#                          (11_ED4_power_forest_meta.R)
#   b — age_meta           REML pooled age-at-presentation meta-analysis
#                          forest (11_ED4_power_forest_meta.R)
#   c — rare_variants      Rare-variant (Other KRAS + BRAF V600E) vs.
#                          panel-negative KM (09_F1_km_age.R, save_km_panel
#                          pointer RDS reassembled by compose_figure)
#
# Layout: 1 column × 3 rows. Each cell 6.60 × 2.03 in (inner); composer
# adds 0.3 in outer margin per side → final composite 7.20 × 6.69 in
# (Nature double-col x max height).
#
# ED number is resolved at runtime by the citation-order resolver in
# panel_assignments.R / resolve_panels.R. Expected slot: ED5, with all
# subsequent ED figures (km_diagnostics, vaf_deep_dive, etc.) shifted
# +1 from their previous positions.
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("analysis", "pipeline", "helpers", "compose_figure.R"))
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa <- load_panel_assignments()

panel_root <- here("results", "ExtendedData", "ed_age_meta_variants")
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)

panel_rds <- list(
  rupture_meta  = file.path(panel_root, "rupture_meta.rds"),
  age_meta      = file.path(panel_root, "age_meta.rds")
)

# v6.22 (2026-05-20): rare_variants moved to ed_km_diagnostics. This
# group now carries the two phenotype meta-forests only. Each forest is
# tall enough to read comfortably in a 1-col x 2-row composite at half
# the Nature max height.
LAYOUT <- "A\nB"

AREA_TO_TOKEN <- c(
  A = "rupture_meta",
  B = "age_meta"
)

ROW_HEIGHTS  <- c(3.045, 3.045)   # 6.09 / 2 = 3.045 each
COMPOSITE_W  <- 6.60

compose_figure(
  panel_rds      = panel_rds,
  layout         = LAYOUT,
  area_to_token  = AREA_TO_TOKEN,
  out_stem       = file.path(panel_root,
                             sprintf("ED%02d_composite",
                                     .pa$group_number$EDFig[["ed_age_meta_variants"]])),
  width_in       = COMPOSITE_W,
  height_in      = sum(ROW_HEIGHTS),
  panel_tag_size = 8,
  heights        = ROW_HEIGHTS
)
