# 36_F4_assemble.R — Figure 4 composer (precision_medicine).
# -----------------------------------------------------------------------------
# Assembles Fig 4 from its component panels so it is reproducible like every
# other main figure — no more hand-assembly in Illustrator (which is what let
# the embedded waterfall go stale). Panels, in registry order:
#   a  biopsy_schematic     — Laura Roy illustration   (external PNG)
#   b  biopsy_angiograms     — DSA images               (external PNG)
#   c  biopsy_tapestation    — TapeStation gel          (external PNG)
#   d  biopsy_dpcr_waterfall  — CODE panel (this pipeline; read from RDS)
#   e  framework_diagram      — Laura Roy illustration  (external PNG)
#
# Panels a/b/c are illustration / wet-lab assets that no script can regenerate.
# They must be exported as SEPARATE images (they are currently fused inside the
# legacy Fig4_liquid_biopsy.png). Drop them at the PANELS paths below; until
# then the composer renders labelled placeholders so the layout is visible and
# the script runs cleanly in run_all.R.
#
# Output: results/Figure4/Fig04_precision_medicine.{png,pdf}
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(here); library(ggplot2); library(patchwork); library(cowplot)
  library(magick); library(grid); library(ggtext)
})
source(here("analysis", "helper_scripts", "utils.R"))

ext        <- here("manuscript", "external_figures")
panels_dir <- file.path(ext, "panels")   # drop separate a/b/c exports here

# Component panel image sources (a/b/c need separate exports; e exists today).
PANEL_IMG <- list(
  biopsy_schematic   = file.path(panels_dir, "biopsy_schematic.png"),
  biopsy_angiograms  = file.path(panels_dir, "biopsy_angiograms.png"),
  biopsy_tapestation = file.path(panels_dir, "biopsy_tapestation.png"),
  framework_diagram  = file.path(ext, "Fig4_precision_framework.png")
)
waterfall_rds <- here("results", "Figure4", "biopsy_dpcr_waterfall",
                      "biopsy_dpcr_waterfall.rds")

# A real materialised raster is >1 KB; a git-LFS pointer or missing file is not.
.is_real_image <- function(p) file.exists(p) && file.info(p)$size > 1024L

panel_from_image <- function(path, label) {
  if (.is_real_image(path)) {
    ggdraw() + draw_image(magick::image_read(path))
  } else {
    ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
               fill = "grey96", colour = "grey70", linetype = "dashed") +
      annotate("text", x = 0.5, y = 0.5, family = NM$font_family,
               size = NM$text$body_mm, colour = "grey40",
               label = paste0("[", label, "]\nexport separate asset →\n",
                              "manuscript/external_figures/panels/")) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void()
  }
}

p_a <- panel_from_image(PANEL_IMG$biopsy_schematic,   "biopsy_schematic")
p_b <- panel_from_image(PANEL_IMG$biopsy_angiograms,  "biopsy_angiograms")
p_c <- panel_from_image(PANEL_IMG$biopsy_tapestation, "biopsy_tapestation")
p_d <- if (file.exists(waterfall_rds)) readRDS(waterfall_rds) else
         panel_from_image("", "biopsy_dpcr_waterfall")
p_e <- panel_from_image(PANEL_IMG$framework_diagram,  "framework_diagram")

# Layout: vertical stack (default). Adjust `design`/`heights` to the intended
# arrangement once the real a/b/c assets and target layout are confirmed.
design <- "A\nB\nC\nD\nE"
composite <- p_a + p_b + p_c + p_d + p_e +
  plot_layout(design = design, heights = c(1, 1, 0.8, 1, 1.3)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e"))) &
  theme(plot.tag = element_text(size = NM$label_pt, face = "bold",
                                family = NM$font_family))

out_dir <- here("results", "Figure4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_stem <- file.path(out_dir, "Fig04_precision_medicine")
.W <- 7.20; .H <- 9.50
ggsave(paste0(out_stem, ".png"), composite, width = .W, height = .H, dpi = 200)
ggsave(paste0(out_stem, ".pdf"), composite, width = .W, height = .H,
       device = cairo_pdf, family = NM$font_family)

missing <- names(PANEL_IMG)[!vapply(PANEL_IMG, .is_real_image, logical(1))]
cat(sprintf("  ✓ Fig 4 composite: %s.{png,pdf}  (%.2f × %.2f in)\n",
            out_stem, .W, .H))
if (length(missing))
  message("  ⚠ placeholder panels (need separate exports): ",
          paste(missing, collapse = ", "))
