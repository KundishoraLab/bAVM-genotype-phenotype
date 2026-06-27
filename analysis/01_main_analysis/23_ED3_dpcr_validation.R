# 23_ED3_dpcr_validation.R — Extended Data figure: multiplex dPCR validation.
#
# Six per-assay partition-amplitude scatter panels (dpcr_scatter_g12d / g12v /
# g12c / g12a / v600e / tert). Each panel has two facets:
#   Left  — "Positive control": real QIAcuity partition-level RFU, the
#           Run 3 FFPE tumour well that lights this assay's channel.
#   Right — "NTC":              real non-template control (well H3), same
#           channel — the negative control per Andy (2026-05-26), replacing
#           the former synthetic-germline right facet.
#
# Input: upstream_path("scrna", "data", "dpcr_real", "edfig3_dpcr_real.csv")
#   De-identified extract written by avm-variant-detection/scripts/
#   09_figure_generation/edfig3_dpcr_validation/extract_real_dpcr.py from the
#   "UAB Priority Samples Run 3" QIAcuity exports.
#   Schema: assay, channel, variant_label, sample_role, partition_id,
#           fluorescence_rfu, threshold, call
#   sample_role ctrl|ntc ; call positive|negative.
#
# Channel->variant map (fixed assay design; confirmed from Run 20 genotyped
# controls): G=G12D, Y=G12V, O=BRAF V600E, R=G12C, C=G12A, Fr=TERT.
#
# Outputs: results/ExtendedData/ed_dpcr_validation/<key>/<key>.{png,pdf,rds}
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(here); library(ggtext)
})
source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "save_panel.R"))

out_root <- here("results", "ExtendedData", "ed_dpcr_validation")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(subdir, name, plot, w, h)
  save_panel_impl(file.path(out_root, subdir), name, plot, w, h, device = "cairo")

csv_path <- upstream_path("scrna", "data", "dpcr_real", "edfig3_dpcr_real.csv")
if (!file.exists(csv_path))
  stop("Real dPCR extract not found at ", csv_path,
       "\nRun avm-variant-detection/scripts/09_figure_generation/",
       "edfig3_dpcr_validation/extract_real_dpcr.py first.")

dpcr <- read_csv(csv_path, show_col_types = FALSE)

# ── Per-assay metadata ──────────────────────────────────────────────────────
# assay key (matches the extract + panel_registry), ggtext-markdown panel title
# (italic gene + superscript variant · dye channel), positive-partition colour,
# and output subdir / panel key. Colour is keyed by the assay's FLUOROPHORE
# (PAL_FLUOROPHORE = real emission-spectrum colours), so each panel reads as its
# dye channel rather than its variant tier — see utils.R for the rationale.
ASSAY_META <- tibble::tribble(
  ~assay,        ~title_md,                         ~mut_color,                  ~panel_key,
  "KRAS_G12D",   "*KRAS*<sup>G12D</sup> · FAM",     PAL_FLUOROPHORE[["FAM"]],    "dpcr_scatter_g12d",
  "KRAS_G12V",   "*KRAS*<sup>G12V</sup> · HEX",     PAL_FLUOROPHORE[["HEX"]],    "dpcr_scatter_g12v",
  "KRAS_G12C",   "*KRAS*<sup>G12C</sup> · ROX",     PAL_FLUOROPHORE[["ROX"]],    "dpcr_scatter_g12c",
  "KRAS_G12A",   "*KRAS*<sup>G12A</sup> · Cy5",     PAL_FLUOROPHORE[["Cy5"]],    "dpcr_scatter_g12a",
  "BRAF_V600E",  "*BRAF*<sup>V600E</sup> · TAMRA",  PAL_FLUOROPHORE[["TAMRA"]],  "dpcr_scatter_v600e",
  "TERT",        "*TERT* · Cy5.5",                  PAL_FLUOROPHORE[["Cy5.5"]],  "dpcr_scatter_tert"
)

call_levels <- c("negative", "positive")
call_size   <- c(negative = 0.35, positive = 0.70)
call_alpha  <- c(negative = 0.45, positive = 0.90)

# ── Build one panel ─────────────────────────────────────────────────────────
build_panel <- function(meta_row) {
  raw <- dpcr %>%
    filter(assay == meta_row$assay) %>%
    mutate(facet_cond = factor(
      if_else(sample_role == "ctrl", "Positive control", "NTC"),
      levels = c("Positive control", "NTC"))) %>%
    group_by(facet_cond) %>%
    mutate(part_x = row_number()) %>%
    ungroup()

  gate <- unique(raw$threshold)
  if (length(gate) != 1L)
    stop(sprintf("[build_panel] non-unique threshold for %s: %s",
                 meta_row$assay, paste(gate, collapse = ", ")), call. = FALSE)

  # Keep every positive; downsample the negative cloud for a light vector PDF.
  set.seed(42)
  neg <- raw %>% filter(call == "negative") %>%
    group_by(facet_cond) %>% slice_sample(n = 1500) %>% ungroup()
  sub <- bind_rows(raw %>% filter(call == "positive"), neg) %>%
    mutate(call = factor(call, levels = call_levels)) %>%
    arrange(call)  # negatives first, positives drawn on top

  pal <- c(negative = "#D0D0D0", positive = meta_row$mut_color)

  ggplot(sub, aes(part_x, fluorescence_rfu, colour = call,
                  size = call, alpha = call)) +
    geom_point(stroke = 0) +
    ref_hline(yintercept = gate, kind = "threshold") +
    facet_wrap(~ facet_cond, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = pal, breaks = call_levels, drop = FALSE, name = NULL) +
    scale_size_manual(values = call_size, guide = "none") +
    scale_alpha_manual(values = call_alpha, guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.08))) +
    scale_x_continuous(breaks = c(0, 10000, 20000),
                       labels = c("0", "10,000", "20,000"),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = meta_row$title_md, x = "Partition", y = "Fluorescence (RFU)") +
    # Per-panel colour legend stripped; composer builds one shared legend in
    # the right-margin gutter (positive = variant colour, negative = grey).
    guides(colour = "none") +
    theme_nature_panel() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position  = "none",
          plot.title = element_markdown(size = NM$body_pt, family = NM$font_family))
}

# ── Render six panels (dims match the Nature-spec 2x3 composite grid) ─────────
for (i in seq_len(nrow(ASSAY_META))) {
  meta_row <- ASSAY_META[i, ]
  p <- build_panel(meta_row)
  save_panel(meta_row$panel_key, meta_row$panel_key, p, w = 3.30, h = 2.03)
}

message("[ed_dpcr_validation] wrote 6 real-data panels to ", out_root)
