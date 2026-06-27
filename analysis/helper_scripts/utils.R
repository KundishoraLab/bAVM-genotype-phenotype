# utils.R — Shared helper functions for bAVM genotype-phenotype analysis
#
# Sourced by both analysis/00_data_prep/ and analysis/01_main_analysis/ scripts.
# Add reusable functions here to avoid duplication.

# ── Canonical random-number seed ─────────────────────────────────────────────
# Single source of truth for every set.seed() call in the manuscript pipeline.
# Audit 2026-05-12 (F13): previously each producer hardcoded its own seed
# (20260424L for risk-score, 20260429L for table1, 20260420 for vaf_outlier,
# 42 for cluster audits) and several Monte-Carlo Fisher calls
# (`fisher.test(..., simulate.p.value = TRUE)`) ran without any seed at all.
# Producers must now `set.seed(MASTER_SEED)` immediately before every
# Monte-Carlo or bootstrap operation, and pass MASTER_SEED through to
# helpers that need their own seed argument.
MASTER_SEED <- 20260424L

# ── Canonical location × genotype interaction fit ────────────────────────────
# Audit 2026-05-12 (F11 deep): both `13_F2_genotype_phenotype.R` (SuppTable08, broad
# 8-outcome family) and `17_ED8_anatomy.R` (ED Fig 8 panel C, focused 3-outcome
# family) need the per-(lobe, outcome) interaction P from
#
#     outcome ~ geno_binary * lobe
#
# fit on `df` rows that have a non-NA outcome, non-NA lobe indicator, and
# non-NA geno_binary. Previously each producer recomputed this independently,
# with subtly different sample filters and reference-level conventions, so
# the same parietal × rupture cell could ship two different raw P values
# AND two different beta signs across the manuscript. This helper centralises
# both, sets the geno_binary reference explicitly to "Panel-negative" so
# the interaction beta sign is stable (positive ⇒ variant-positive carries
# the elevated effect on parietal vs non-parietal), and returns a one-row
# tibble per (data, outcome, lobe).
fit_loc_geno_interaction <- function(data,
                                     outcome_var,
                                     loc_var,
                                     outcome_type = c("binary", "continuous"),
                                     min_n        = 40L,
                                     geno_var     = "geno_binary",
                                     geno_ref     = "Panel-negative") {
  outcome_type <- match.arg(outcome_type)
  stopifnot(geno_var    %in% names(data),
            outcome_var %in% names(data),
            loc_var     %in% names(data))

  s <- data %>%
    dplyr::filter(!is.na(.data[[outcome_var]]),
                  !is.na(.data[[loc_var]]),
                  !is.na(.data[[geno_var]]))

  .na_result <- function(method)
    tibble::tibble(beta = NA_real_, se = NA_real_, p = NA_real_,
                   n = nrow(s), method = method)

  if (nrow(s) < min_n)                          return(.na_result("skipped_n_below_min"))
  if (dplyr::n_distinct(s[[loc_var]])  < 2)     return(.na_result("skipped_loc_constant"))
  if (dplyr::n_distinct(s[[geno_var]]) < 2)     return(.na_result("skipped_geno_constant"))

  # Explicit reference level so the interaction beta sign is stable across
  # producers and across re-runs. "Panel-negative" is the canonical factor
  # level (renamed from legacy "Genotype-negative" 2026-05-16) and the
  # conventional comparator in §4 prose ("variant-positive lesions vs.
  # panel-negative").
  s[[geno_var]] <- stats::relevel(factor(s[[geno_var]]), ref = geno_ref)

  fam <- if (outcome_type == "binary") stats::binomial() else stats::gaussian()
  fml <- stats::as.formula(paste(outcome_var, "~", geno_var, "*", loc_var))
  m   <- tryCatch(stats::glm(fml, data = s, family = fam),
                  error = function(e) NULL)
  if (is.null(m)) return(.na_result("glm_fit_failed"))

  cs <- summary(m)$coefficients
  ir <- grep(":", rownames(cs))
  if (length(ir) == 0L) return(.na_result("no_interaction_term"))

  tibble::tibble(
    beta   = cs[ir[1], 1],
    se     = cs[ir[1], 2],
    p      = cs[ir[1], 4],
    n      = nrow(s),
    method = "glm_geno_interaction"
  )
}

# ── Cross-submodule path helpers ─────────────────────────────────────────────
# Manuscript producers consume scRNA + spatial outputs directly from the
# sibling submodule clones (per STRUCTURE.md R3) instead of via committed
# mirrors. Single source of truth, never stale.
#
# upstream_path("scrna",   "data/results", "per_cell_results.tsv")
#   -> ../avm-variant-detection/data/results/per_cell_results.tsv
# upstream_path("spatial", "results/07_cell_communication", "lr_comm_top_hits_ctscaled.tsv")
#   -> ../avm-spatial-tx/results/07_cell_communication/lr_comm_top_hits_ctscaled.tsv
#
# The clone topology is the brain-avm parent monorepo (recursive submodule
# clone). A standalone genotype-phenotype clone WITHOUT the siblings will
# fail loudly at read time — by design. Don't silently substitute stale
# committed snapshots.
upstream_path <- function(repo, ...) {
  base <- switch(
    repo,
    scrna   = "avm-variant-detection",
    spatial = "avm-spatial-tx",
    stop("unknown upstream repo: ", repo, " (use 'scrna' or 'spatial')")
  )
  file.path(here::here(), "..", base, ...)
}

# ── Packages commonly needed ─────────────────────────────────────────────────
load_common_packages <- function() {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
    library(janitor)
    library(knitr)
    library(kableExtra)
  })
}

# ── AVM panel allowlist (mirror of scRNA repo's scripts/utils/avm_panel.R) ──
# Canonical (gene, variant) tuples for somatic AVM driver mutations
# (Nikolaev 2018 + Al-Olabi 2018 hotspot set). Used by every scRNA consumer
# to defensively drop off-panel calls — V600V_syn synonymous −1 bp neighbour
# of BRAF V600E, TEK R849W, raw chr* positional noise, opportunistic
# ClinVar/COSMIC scans — that may slip through the upstream scRNA repo
# preprocessor. Run-time enforcement here means the canonical inventory is
# reproducible from any TSV path, not just the one-shot preprocessor output.
AVM_PANEL_GENE_VARIANT <- c(
  "KRAS G12D", "KRAS G12V", "KRAS G12A", "KRAS G12C", "KRAS Q61H",
  "BRAF V600E"
)

filter_avm_panel <- function(df, mutation_col = "mutation", verbose = TRUE) {
  stopifnot(mutation_col %in% names(df))
  keep <- df[[mutation_col]] %in% AVM_PANEL_GENE_VARIANT
  if (verbose && sum(!keep) > 0) {
    dropped_tally <- as.data.frame(table(df[[mutation_col]][!keep]))
    message(sprintf("filter_avm_panel: %d → %d rows (%d off-panel dropped)",
                    nrow(df), sum(keep), sum(!keep)))
    for (i in seq_len(nrow(dropped_tally))) {
      message(sprintf("  dropped: %s  (n=%d)",
                      as.character(dropped_tally$Var1[i]),
                      dropped_tally$Freq[i]))
    }
  }
  df[keep, , drop = FALSE]
}

# ── Genotype display-label helper ────────────────────────────────────────────
# Producers carry the raw factor level "Negative" (variant_group / geno_binary
# / geno_variant). The canonical display label everywhere in the manuscript
# is "Panel-negative". Use this helper at the render layer to swap the level
# string without touching the upstream factor:
#
#   ggsurvplot(..., legend.labs = display_levels(dat$variant_group))
#   ggplot(...) + scale_x_discrete(labels = display_levels)
#
# Promoted from 09_F1_km_age.R (2026-05-16) so every ED producer can share one
# definition instead of redefining it locally.
display_levels <- function(f) {
  sub("^Negative$", "Panel-negative", levels(f))
}

# Convenience: apply the rename to a factor's levels in-place (returns a new
# factor) so producers can call it before plotting categorical aesthetics
# that don't accept a `labels` formatter (e.g. facet strips, dumbbell groups).
relabel_geno_factor <- function(f) {
  if (!is.factor(f)) f <- factor(f)
  levels(f) <- display_levels(f)
  f
}

# ── Paper palette (consistent with survey project) ───────────────────────────
PAPER_PALETTE <- c("#000000", "#E69F00", "#56B4E9", "#CC79A7")

# ── Colorblind-safe categorical palette (Wong/Okabe-Ito 8-colour) ───────────
# 2026-05-20 migration: replaced the prior IBM Carbon CB family with the
# Wong palette (Okabe & Ito 2008; popularised in Wong 2011 *Nature Methods*
# 8:441). The eight tokens are colorblind-safe under deuteranopia,
# protanopia, and tritanopia simulations, and reproduce cleanly in both
# RGB and CMYK without out-of-gamut clipping. Hex codes locked to the
# canonical RGB values supplied by the journal palette:
#
#   Black           #000000   (  0,   0,   0)
#   Orange          #E69F00   (230, 159,   0)
#   Sky blue        #56B4E9   ( 86, 180, 233)
#   Bluish green    #009E73   (  0, 158, 115)
#   Yellow          #F0E442   (240, 228,  66)
#   Blue            #0072B2   (  0, 114, 178)
#   Vermillion      #D55E00   (213,  94,   0)
#   Reddish purple  #CC79A7   (204, 121, 167)
#
# Defined here as W_* constants so downstream palettes reference them by
# semantic name. Greys (#737373 medium / #5F5F5F dark / #B0B0B0 light)
# remain reserved for null/reference strata (Panel-negative, score 0,
# Literature source) because the Wong palette intentionally omits a
# neutral grey.
W_BLACK         <- "#000000"
W_ORANGE        <- "#E69F00"
W_SKYBLUE       <- "#56B4E9"
W_BLUISHGREEN   <- "#009E73"
W_YELLOW        <- "#F0E442"
W_BLUE          <- "#0072B2"
W_VERMILLION    <- "#D55E00"
W_REDDISHPURPLE <- "#CC79A7"

# TIERED COLOUR CONTRACT — colours in Tier A are reserved for the four
# named variants and MUST NOT be reused for any other category in any
# figure that displays variants. Tier B (binary outcome) is reserved
# for Variant-positive vs Panel-negative. Tier C (risk score) uses
# distinct hues outside Tier A/B. Tier D (study, sample sub-strata)
# uses the remaining Wong palette. Grey is the only colour permitted
# to recur across tiers (Panel-negative, score 0, Literature, etc.).
#
# COLOUR-BLIND RULE — the Wong palette is internally validated; the
# only additional constraint is to avoid using both Vermillion and
# Reddish purple in the same categorical scale when their hues
# converge under tritanopia. Acceptable pairings inside one panel:
# Blue + Orange, Blue + Vermillion, Bluish green + Reddish purple,
# Sky blue + Vermillion, Black + any.

# ── Tier A: per-variant (locked) ──────────────────────────────────────────
# Variant-tier factors (variant_group, geno_variant) still carry the bare
# "Negative" level internally; render-layer helpers relabel it as
# "Panel-negative" in legends. Both keys point at the same grey so palette
# lookup works pre- and post-relabel.
PAL_VARIANT <- c(
  "KRAS G12D"      = W_BLUE,           # #0072B2
  "KRAS G12V"      = W_VERMILLION,     # #D55E00
  "Other KRAS"     = W_ORANGE,         # #E69F00
  "BRAF"           = W_REDDISHPURPLE,  # #CC79A7
  "Negative"       = "#737373",        # factor-level key (variant_group)
  "Panel-negative" = "#737373"         # canonical display label
)

# ── Tier B: binary outcome (Mut+ vs Panel-neg) ────────────────────────────
# Variant-positive is Bluish green (Wong's CB-safe green), NOT blue —
# blue is reserved for KRAS G12D in Tier A.
PAL_BINARY <- c(
  "Variant-positive" = W_BLUISHGREEN, # #009E73
  "Panel-negative"    = "#737373"      # medium grey (null tier)
)

# Detailed variant (for landscape / oncoprint — rare variants split out).
# Extends Tier A; lighter shades for rarer same-family variants are mixed
# with white to keep the bar landscape reading as variant-family ramps.
.mix_white <- function(hex, frac = 0.45) {
  rgb_to_hex <- function(r, g, b) sprintf("#%02X%02X%02X", r, g, b)
  rgb_in <- col2rgb(hex)[, 1]
  rgb_to_hex(
    as.integer(round(rgb_in[1] + (255 - rgb_in[1]) * frac)),
    as.integer(round(rgb_in[2] + (255 - rgb_in[2]) * frac)),
    as.integer(round(rgb_in[3] + (255 - rgb_in[3]) * frac))
  )
}
PAL_DETAILED <- c(
  "KRAS G12D"      = W_BLUE,
  "KRAS G12V"      = W_VERMILLION,
  "KRAS G12C"      = W_ORANGE,
  "KRAS G12A"      = .mix_white(W_ORANGE, 0.45),
  "KRAS Q61H"      = W_SKYBLUE,
  "KRAS dup"       = W_BLUISHGREEN,
  "BRAF V600E"     = W_REDDISHPURPLE,
  "BRAF Q636X"     = .mix_white(W_REDDISHPURPLE, 0.45),
  "Negative"       = "#B0B0B0",
  "Panel-negative" = "#B0B0B0",
  "Unassigned variant-positive" = "#999999",
  "No tissue"      = "#E8E8E8"
)

# ── dPCR dye-channel emission colours (ED Fig 3 + Fig 4d only) ──────────────
# These two figures colour positive partitions by the assay's FLUOROPHORE
# emission spectrum (which optical channel lit up), NOT by the variant tier —
# the panels are about multiplex dye channels, so the colour communicates the
# dye, not the mutation. This is the ONE sanctioned departure from the Tier-A
# variant contract above and applies ONLY to the dPCR partition scatters; the
# variant tiers remain locked for every variant-bearing figure (Fig 1b, etc.).
# Per Jonah (2026-05-27): match real emission maxima — warm dyes ordered by λ
# (HEX 556 < TAMRA 580 < ROX 608). NB: this warm ramp is not colour-blind safe,
# but colour is constant within each panel (positive vs grey negative) and every
# panel is titled with its dye, so channel identity never rests on colour alone.
PAL_FLUOROPHORE <- c(
  "FAM"   = "#00A651",  # ~517 nm  green
  "HEX"   = "#F2A900",  # ~556 nm  yellow-orange
  "TAMRA" = "#F57C00",  # ~580 nm  orange
  "ROX"   = "#E8401C",  # ~608 nm  red-orange
  "Cy5"   = "#C8102E",  # ~668 nm  far-red / crimson
  "Cy5.5" = "#7D0633"   # ~694 nm  infrared / deep crimson
)

# SM grade (blue sequential — light to dark). Sequential not categorical,
# perceptually uniform; kept on the existing Brewer Blues ramp.
PAL_SM <- c(
  "I"   = "#D1E5F0",
  "II"  = "#92C5DE",
  "III" = "#4393C3",
  "IV"  = "#2166AC",
  "V"   = "#053061"
)

# ── Tier D: study palette (8 contributing series) ─────────────────────────
# Wong gives 8 distinct hues. Map 1-to-1 across the 8 series. Cross-figure
# reuse with PAL_BINARY's Bluish green (BCH = Mut+ green) is acceptable
# because study composition panels never co-occur with binary-outcome
# panels in the same figure.
PAL_STUDY <- c(
  "BCH"      = W_BLUISHGREEN,
  "UAB"      = W_SKYBLUE,
  "CHOP"     = W_BLUE,
  "Nikolaev" = W_REDDISHPURPLE,
  "Priemer"  = W_VERMILLION,
  "Hong"     = W_ORANGE,
  "Goss"     = W_BLACK,
  "Gao"      = W_YELLOW
)

# ── Tier D: sample type palette ───────────────────────────────────────────
PAL_SAMPLE <- c(
  "FFPE"          = W_SKYBLUE,
  "Fresh/Frozen"  = W_REDDISHPURPLE,
  "Literature"    = "#737373"
)

# Heatmap gradients (use for scale_fill_gradient)
PAL_HEAT_LOW  <- "#F7F7F7"
PAL_HEAT_HIGH <- W_BLUE

# ── Tier F: anatomy-locked palette ────────────────────────────────────────
# Anatomy panels use a Yellow -> Vermillion sequential ramp (warm hues
# outside Tier A's KRAS Blue / BRAF Reddish-purple), with Reddish purple
# as the highlight hue for a called-out lobe. Both endpoints are Wong
# colours so the ramp stays colorblind-safe.
PAL_ANATOMY <- list(
  LOW       = "#FFF5CC",       # very light yellow (tint of W_YELLOW)
  HIGH      = W_VERMILLION,    # warm anchor
  ACCENT    = W_VERMILLION,
  HIGHLIGHT = W_REDDISHPURPLE  # called-out lobe (parietal x rupture)
)

# EC subtype palette — arteriovenous axis (red artery -> purple capillary
# -> blue vein). Biologically themed gradient; intentionally NOT mapped
# to Wong because the colour encoding carries anatomy, not category.
EC_SUBTYPE_COLORS <- c(
  "Large artery"         = "#8B1E1E",
  "Artery"               = "#C0392B",
  "Arteriole"            = "#E74C3C",
  "Capillary"            = "#9B59B6",
  "Angiogenic capillary" = "#6C3483",
  "Venule"               = "#2471A3",
  "Vein"                 = "#1A5276",
  "Large vein"           = "#0E2F43",
  "EndoMT"               = "#7F8C8D",
  "Stem-to-EC"           = "#95A5A6",
  "Proliferating cell"   = "#BDC3C7",
  "Mitochondrial"        = "#D5DBDB",
  "Lymphatic"            = "#1ABC9C"
)

# Single-cell variant palette keyed on `variant_status` labels in
# avm-variant-detection. Uses Tier A hues so sc UMAPs match Fig 1B.
VARIANT_SC_COLORS <- c(
  "G12D"     = W_BLUE,
  "G12V"     = W_VERMILLION,
  "G12C"     = W_ORANGE,
  "G12A"     = .mix_white(W_ORANGE, 0.45),
  "Q61H"     = W_SKYBLUE,
  "V600E"    = W_REDDISHPURPLE,
  "Wildtype" = "grey85"
)

# ── KM-context palettes ────────────────────────────────────────────────────
# Same Tier A hues for KRAS arms; the Panel-negative stratum uses a
# slightly DARKER grey (#5F5F5F vs PAL_*'s #737373) because KM curves
# render as thin lines on grey-93 panel grids and a darker negative grey
# holds contrast in line-only contexts.

# 3-arm KM palette (G12D / G12V / Panel-negative).
PAL_KM <- c(
  "KRAS G12D"      = W_BLUE,
  "KRAS G12V"      = W_VERMILLION,
  "Negative"       = "#5F5F5F",
  "Panel-negative" = "#5F5F5F"
)

# Rare-variant companion palette (Other KRAS + BRAF V600E vs negative).
PAL_RARE <- c(
  "Other KRAS"     = W_ORANGE,
  "BRAF V600E"     = W_REDDISHPURPLE,
  "Negative"       = "#5F5F5F",
  "Panel-negative" = "#5F5F5F"
)

# Binary KM / density palette — Variant-positive uses Bluish green.
PAL_BINARY_KM <- c(
  "Variant-positive" = W_BLUISHGREEN,
  "Panel-negative"    = "#5F5F5F"
)

# ── KM x-axis display cutoff ───────────────────────────────────────────────
# Display-only upper bound (years). The Kaplan-Meier fit uses the full age
# range; only the rendered window is truncated.
KM_AGE_XLIM_MAX <- 60

# ── Tier E: clinical event / alarm token ───────────────────────────────────
# Vermillion as the universal alarm convention. CB-safe under all three
# common dichromacy types (Wong-validated).
PAL_EVENT <- c(
  "Event"    = W_VERMILLION,
  "No event" = "#737373"
)

# ── Tier C: risk-score palette (integer score 0/1/2) ───────────────────────
# Score 0 reuses grey because it represents the null/lowest-risk tier.
# Score 1 = Sky blue, score 2 = Reddish purple. Categorical contrast,
# NOT a sequential ramp.
PAL_SCORE <- c(
  "0" = "#737373",
  "1" = W_SKYBLUE,
  "2" = W_REDDISHPURPLE
)

# ── Reference / threshold line styling ─────────────────────────────────────
# Uniform appearance for vlines / hlines marking null references (OR=1,
# β=0), p-value thresholds, secondary thresholds (Bonferroni), and
# data-driven median lines. Use ref_vline() / ref_hline() with a `kind`
# argument; pass extra args (e.g. `data =`, `aes(...)`) via `...`.
#
# Conventions:
#   null      OR=1 / β=0 / x=0 reference. Dashed, grey60, lw 0.5.
#   threshold p-value cutoff (e.g. -log10(0.05)). Same look as null.
#   bonf      secondary / Bonferroni threshold. Dotted, grey50, lw 0.5.
#   median    median markers in densities (colour aes-mapped per group).
#             Dashed, lw 0.8 — chunkier so it cuts through the density curve.
#
# Hard QC cutoffs that carry semantic colour (e.g. red Q30 floor in
# 17_ed_scrna_phred_qc.R) intentionally bypass these tokens.
REF_LINE <- list(
  # v6.48 (2026-05-21): linewidths dropped to match Fig 1 line weight (0.4)
  # for native Nature-spec rendering. Prior 1.0–1.6 values were sized for
  # the legacy 14-in producer canvases; at 7.20-in Nature double-col print
  # they read as heavy dotted scaffolding overpowering the data ink.
  null      = list(linetype = "dashed",  colour = "grey60", linewidth = 0.4),
  threshold = list(linetype = "dashed",  colour = "grey60", linewidth = 0.4),
  bonf      = list(linetype = "longdash", colour = "grey35", linewidth = 0.5),
  median    = list(linetype = "dashed",                     linewidth = 0.5)
)

ref_vline <- function(xintercept = NULL, kind = "null", ...) {
  spec <- REF_LINE[[kind]]
  args <- c(spec, list(...))
  if (!is.null(xintercept)) args$xintercept <- xintercept
  do.call(geom_vline, args)
}

ref_hline <- function(yintercept = NULL, kind = "null", ...) {
  spec <- REF_LINE[[kind]]
  args <- c(spec, list(...))
  if (!is.null(yintercept)) args$yintercept <- yintercept
  do.call(geom_hline, args)
}

# ── Typography contract ─────────────────────────────────────────────────────
# Source-of-truth font sizes AND font family for every figure in the
# manuscript. Two scales, one for each render path:
#
#   * theme_avm() / TYPO        → standalone producer PNG/PDF (review path)
#   * theme_avm_native() / NATIVE → submission composite at journal print size
#
# Both themes lock `base_family = "Arial"` (Nature Medicine figure-guideline
# preference) so the typeface is identical across every panel regardless of
# the cairo device default. Bumped 2026-05-16 from base 14/7 to 18/9 in the
# uniform-typography pass — TYPO/NATIVE annotation tokens were scaled
# proportionally (×18/14 ≈ ×1.29) so inline counts, panel tags, and
# JAMA-forest annotations track the larger axis text.
#
# theme_avm() — universal panel theme, base_size = 18, family = Arial
#   Standalone producer files are saved at native ggsave size (no resize),
#   so axis title = 18pt and axis text = ~14.4pt in the final PNG/PDF.
#
# theme_avm_composite() — same base_size, tighter margins for grid packing.
#   Use when a panel is embedded into a multi-panel composite. The reduced
#   plot.margin keeps inter-panel gutters tasteful.
#
# Inline geom_text() / annotate("text") sizes go through TYPO so every
# in-panel annotation (counts, percents, p-values, legend stuffing) reads
# at the same physical size across figures. ggplot::geom_text size is in
# mm; pt-equivalents below assume the standard size.pt = 72/25.4 conversion.
#
# Rule of thumb when adding a new annotation:
#   - "Label" (axis-label sized; counts on bars, OR text in forests):
#       use TYPO$geom_text_label  (~10pt — matches base_size 14 axis text)
#   - "Small" (dense risk-table-like, fine print):
#       use TYPO$geom_text_small  (~8pt — for tightly packed labels)
#   - "Tiny" (last-resort, only when small still doesn't fit):
#       use TYPO$geom_text_tiny   (~6.8pt — at the NatMed floor)
#
# Update the contract here, not in producers, so a single edit re-flows
# every figure consistently.
TYPO <- list(
  # Per-panel tag in patchwork composites (a / b / c / …). Bold throughout.
  # Scales with base_size: ~22pt for a base_size 18 figure.
  panel_tag        = 22,
  # ggplot geom_text sizes (mm). Approx pt = mm × 72 / 25.4 ≈ mm × 2.83.
  geom_text_label  = 4.5,    # ~ 12.7pt
  geom_text_small  = 3.6,    # ~ 10.2pt
  geom_text_tiny   = 3.1,    # ~  8.8pt — last-resort dense annotations only
  # Legend-key squares; lines unit so it scales with theme.
  legend_key_size  = 1.15,
  # ─ Inline annotation hierarchy for hand-rolled forests / tables ──
  # Use these in annotate("text", size = TYPO$annot_*) when building
  # JAMA-style stacked text headers + body rows by hand. Tweak here
  # to reflow every hand-rolled annotated panel. Standard forests
  # (table_forest / table_forest_meta) derive their internal sizes
  # from base_size automatically and don't reference these tokens.
  annot_header_bold = 5.5,   # ~ 15.6pt  — column / row headers (bold)
  annot_body_bold   = 5.3,   # ~ 15.0pt  — bold predictor / label rows
  annot_body_plain  = 4.8,   # ~ 13.6pt  — body row OR (95% CI) text
  annot_axis_text   = 4.2    # ~ 11.9pt  — tick value labels
)

# Lock the font family that every theme inherits. Set as a single constant
# so a future house-style change (e.g. switching to Helvetica or a journal
# requiring Times) touches exactly one line. systemfonts (used by ggplot2
# ≥3.4 and the cairo backends) resolves "Arial" via the OS font lookup; on
# macOS this maps to the installed Arial / Arial.ttf, on Linux build hosts
# it falls back to a Helvetica-class metric-compatible alias.
# v6.50 (2026-05-21): Arial replaces Helvetica as the manuscript-wide
# font family. Nature accepts both; Arial is more universally
# installed on Windows reviewer machines and avoids missing-font
# warnings when reviewers open the PDFs / SVGs in Adobe Illustrator.
AVM_FONT_FAMILY <- "Arial"

# Helper: convert ggplot::element_text size in pt → mm for geom_text. Useful
# when an annotation should track the panel's axis-text size automatically
# (e.g. count text on top of a bar).
typo_pt_to_mm <- function(pt_size) pt_size / ggplot2::.pt

# ── Native typography contract for submission composites ─────────────────────
# Producer scripts save STANDALONE panels at base_size 14 (TYPO above) for
# internal review. The submission composite — the thing Nature Medicine
# actually receives — is rendered separately by analysis/pipeline/helpers/compose_figure.R
# at the journal's exact print footprint, with NO scale-down at the publisher
# stage. That requires a smaller native base_size so axis text lands inside
# Nature Medicine's 5–7pt spec at the on-page width.
#
#   NATIVE$widths_in   — Nature Medicine column targets (single / 1.5 / double)
#   NATIVE$base_size   — ggplot base_size that gives ~5.5pt axis-text and
#                        ~7pt axis-title at native composite render (no resize)
#   NATIVE$panel_tag   — bold a/b/c letter size in pt (matches NatMed style)
#   NATIVE$geom_text_* — mm sizes for inline annotations at native scale
#
# Compose helpers consume this. Producer scripts that keep their inline
# patchwork (i.e. anything not yet migrated to compose_figure) should remain
# on TYPO + theme_avm() — switching is opt-in figure by figure.
NATIVE <- list(
  widths_in        = list(single = 3.50, one_half = 4.72, double = 7.20),
  base_size        = 9,    # axis-title; axis-text ≈ 7.2pt (top of NatMed 5–7pt range)
  panel_tag        = 12,   # bold a/b/c letter, slightly larger than axis-title
  geom_text_label  = 2.5,  # ~ 7.1pt — inline counts / OR labels
  geom_text_small  = 2.05, # ~ 5.8pt — risk-table fine print
  geom_text_tiny   = 1.80, # ~ 5.1pt — last-resort dense annotations
  legend_key_size  = 0.70  # lines unit (smaller than TYPO$legend_key_size)
)

# theme_avm_native(base_size = NATIVE$base_size) — drop-in theme for panels
# rendered inside compose_figure() composites at native column width. Same
# theme contract as theme_avm() (classic + grey93 grid + bold strip + tight
# margins) but calibrated for ~5–7pt rendered text at print size.
# ── theme_nature_panel — composer-side canonical theme ───────────────────────
# 2026-05-20 v6: enforces Nature display-item font sizes (5-7 pt body, 8 pt
# panel-tag bold lowercase) at the SINGLE chokepoint where every panel in
# a composite passes through. Producers save UNTHEMED ggplot objects to
# RDS; the composer reads them and applies this theme via `&` so every
# panel inherits identical font sizes regardless of its cell footprint.
#
#   axis.title  = NM$body_pt (7 pt)        — uniform across panels
#   axis.text   = NM$tick_pt (6 pt)        — one size down per spec
#   legend.text = NM$tick_pt
#   plot.tag    = NM$label_pt (8 pt)       — bold panel letter
#   font.family = NM$font_family           — Arial / Helvetica (Nature spec)
#
# Use INSIDE compose_figure() via `&` to override per-panel themes.
# Producers MAY still apply theme_avm() for review-sized standalone PDFs
# (those are not the final submission artifacts) but should NOT bake
# theme into the panel RDS object that the composer will read.
NM <- list(
  width_in    = 7.20,   # 183 mm double-column
  height_in   = 6.69,   # 170 mm max
  body_pt     = 7,      # axis title / in-plot text
  tick_pt     = 6,      # axis text / legend text
  label_pt    = 8,      # bold lowercase panel tag
  # v6.50 (2026-05-21): Arial replaces Helvetica (Windows-stock,
  # matches AVM_FONT_FAMILY above; Nature accepts both).
  font_family = "Arial",
  # v6.50 — single source of truth for hand-rolled geom_text /
  # annotate("text", ...) sizes. Values are in ggplot2 mm (size = X is
  # text height in mm). pt ≈ mm × 2.83.
  #   body_mm  2.0 mm  ≈ 5.7 pt   — in-cell labels, dense annotations
  #   small_mm 1.8 mm  ≈ 5.1 pt   — dense column-header / micro-labels
  #   tiny_mm  1.6 mm  ≈ 4.5 pt   — last-resort below Nature's 5-pt
  #                                  floor; use sparingly
  # Producers should NEVER use literal mm values; reference NM$text$*
  # instead so a future Nature spec change touches one place.
  text        = list(body_mm = 2.0, small_mm = 1.8, tiny_mm = 1.6)
)

theme_nature_panel <- function() {
  theme_classic(base_size = NM$body_pt, base_family = NM$font_family) +
    theme(
      axis.text       = element_text(size = NM$tick_pt,  colour = "black"),
      axis.title      = element_text(size = NM$body_pt,  colour = "black"),
      legend.text     = element_text(size = NM$tick_pt),
      legend.title    = element_text(size = NM$body_pt,  face = "bold"),
      legend.key.size = unit(0.28, "cm"),
      strip.text      = element_text(size = NM$body_pt,  face = "bold"),
      panel.grid.major = element_line(colour = "grey93"),
      plot.margin     = margin(2, 4, 2, 2, "pt"),
      # Preserve the patchwork-stamped tag against theme_classic's
      # %+replace% semantics (theme_classic sets plot.tag = NULL).
      plot.tag        = element_text(size = NM$label_pt, face = "bold",
                                     family = NM$font_family)
    )
}

theme_avm_native <- function(base_size = NATIVE$base_size) {
  theme_classic(base_size = base_size, base_family = AVM_FONT_FAMILY) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(color = "grey40", size = base_size - 1),
      strip.text       = element_text(face = "bold", size = base_size),
      # 2026-05-20 uniform font hierarchy across the figure:
      #   axis.title       = base_size       (e.g., 9 pt)
      #   axis.text        = base_size - 1   (e.g., 8 pt — "one size down" from labels)
      # Ensures x/y-label sizes are consistent across panels and tick
      # labels are uniformly one step smaller, regardless of any
      # per-panel theme overrides the producer applied earlier.
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.position  = "bottom",
      legend.key.size  = unit(NATIVE$legend_key_size, "lines"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1),
      panel.grid.major = element_line(colour = "grey93"),
      plot.margin      = margin(3, 4, 3, 4),
      # Preserve panel-tag visibility. When this theme is applied via
      # patchwork's `&` in compose_figure(), theme_classic's %+replace%
      # semantics set plot.tag = NULL on each child plot, which kills the
      # panel-letter rendered by plot_annotation(tag_levels = "a"). Set
      # plot.tag explicitly here so the tag survives the override.
      plot.tag         = element_text(face = "bold", size = NATIVE$panel_tag)
    )
}

# ── KM-panel save helper ─────────────────────────────────────────────────────
# survminer::ggsurvplot(..., risk.table = TRUE) returns an object with
# $plot (the curve ggplot) and $table (the at-risk-table ggplot). Producers
# historically combined the two via patchwork into a single object and
# saved that as RDS — but compose_figure() can't re-theme the combined
# patchwork at native size (wrap_elements rasterises it as a single grob
# and curve/table text sizes are baked in at the producer's standalone
# dimensions).
#
# save_km_panel() writes THREE RDS files per panel:
#
#   <name>__curve.rds   the $plot ggplot, untouched
#   <name>__table.rds   the $table ggplot, untouched
#   <name>.rds          a small `list(kind = "km_panel", ...)` pointer
#
# compose_figure detects the pointer, reads curve and table separately,
# re-themes each at the composite's native base_size (with the table at
# base_size - 1 so it doesn't dominate), and reassembles at the requested
# `ratio` height split. Standalone PDF/PNG are still written for review.
#
# Use for any survminer ggsurvplot with `risk.table = TRUE`. For
# ggsurvplot_facet (returns a single ggplot, no risk table), use save_trio.
save_km_panel <- function(dir, name, ggsurv, w, h, ratio = c(3, 1)) {
  stopifnot(is.list(ggsurv),
            inherits(ggsurv$plot,  "ggplot"),
            inherits(ggsurv$table, "ggplot"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  composite <- ggsurv$plot / ggsurv$table +
               patchwork::plot_layout(heights = ratio)
  ggsave(file.path(dir, paste0(name, ".pdf")), composite,
         width = w, height = h, device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), composite,
         width = w, height = h, dpi = 300, type = "cairo")
  saveRDS(ggsurv$plot,  file.path(dir, paste0(name, "__curve.rds")))
  saveRDS(ggsurv$table, file.path(dir, paste0(name, "__table.rds")))
  saveRDS(
    list(
      kind  = "km_panel",
      curve = paste0(name, "__curve.rds"),
      table = paste0(name, "__table.rds"),
      ratio = ratio
    ),
    file.path(dir, paste0(name, ".rds"))
  )
}

# ── Custom ggplot theme ──────────────────────────────────────────────────────
# Universal theme. base_size 18 + base_family "Arial" is the contract;
# bumped from 14 → 18 in the 2026-05-16 uniform-typography pass so the
# standalone review files read comfortably at on-screen and slide sizes
# while still meeting NatMed's print-text floor after composite resize.
# Font family locked via AVM_FONT_FAMILY so the typeface no longer depends
# on whatever the cairo device picks as its default.
theme_avm <- function(base_size = 18) {
  theme_classic(base_size = base_size, base_family = AVM_FONT_FAMILY) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(color = "grey40", size = base_size - 1),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.key.size  = unit(TYPO$legend_key_size, "lines"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      panel.grid.major = element_line(colour = "grey93"),
      plot.margin      = margin(6, 8, 6, 8)
    )
}

# ── Composite-context theme ─────────────────────────────────────────────────
# Same typography contract as theme_avm(); tighter plot.margin for packing
# multiple panels into a Fig*_composite.{png,pdf} grid without uncomfortable
# whitespace gutters between panels. Use this at the composite-assembly
# step (e.g. inside patchwork::wrap_plots(...)) so standalone per-panel
# saves can keep theme_avm()'s default margins.
theme_avm_composite <- function(base_size = 14) {
  theme_avm(base_size = base_size) +
    theme(
      plot.margin      = margin(4, 6, 4, 6)
    )
}

# ── Table styling helper ─────────────────────────────────────────────────────
style_table <- function(k, ...) {
  k %>% kable_styling(
    font_size = 12,
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = TRUE, ...
  )
}

# ── JAMA-style horizontal stacked bar with alignment guides ──────────────────
# data: tibble with columns [group, category, pct] where pct sums to 100 per group
# palette: named color vector keyed by category
# n_per_group: named integer vector keyed by group (for y-axis labels)
# min_label_pct: minimum pct to show label inside segment
# ── VAF × age regression scatter — Fig 1 Panel G canonical design ────────────
# Returns a ggplot styled 1-to-1 with Fig 1 Panel G (vaf_age_scatter):
#   geom_point(alpha = 0.3, size = 0.6)
#   geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.2)
#   PAL_KM colour + fill, theme_nature_panel + no_grid
#   5-tick VAF axis (0/2/4/6/8) with .00 accuracy formatter
# `vaf_axis = "y"` matches Fig 1 (VAF on y, age on x). `vaf_axis = "x"`
# is the inverted ED07 vaf_age_rupture_scatter case (VAF on x, age on y).
# Both forms emit identical line/alpha/colour grammar; only the breaks/
# limits attach to the appropriate axis.
vaf_age_scatter_panel <- function(data,
                                  x_var, y_var,
                                  color_var = "variant_group",
                                  palette   = NULL,
                                  vaf_axis  = c("y", "x"),
                                  x_lab     = NULL,
                                  y_lab     = NULL,
                                  vaf_lim   = c(0, 8),
                                  vaf_breaks = c(0, 2, 4, 6, 8),
                                  point_alpha = 0.3,
                                  point_size  = 0.6,
                                  smooth_lw   = 0.5,
                                  smooth_alpha = 0.2) {
  vaf_axis <- match.arg(vaf_axis)
  if (is.null(palette)) palette <- PAL_KM
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]],
                        color = .data[[color_var]])) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_smooth(aes(fill = .data[[color_var]]),
                method = "lm", se = TRUE,
                linewidth = smooth_lw, alpha = smooth_alpha) +
    scale_color_manual(values = palette, name = "Genotype") +
    scale_fill_manual(values  = palette, guide = "none") +
    labs(x = x_lab, y = y_lab) +
    theme_nature_panel() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  vaf_scale <- scale_y_continuous(breaks = vaf_breaks,
                                  minor_breaks = NULL,
                                  labels = scales::label_number(accuracy = 0.01),
                                  limits = vaf_lim)
  if (vaf_axis == "x") {
    vaf_scale <- scale_x_continuous(breaks = vaf_breaks,
                                    minor_breaks = NULL,
                                    labels = scales::label_number(accuracy = 0.01),
                                    limits = vaf_lim)
  }
  p + vaf_scale
}

jama_stacked_bar <- function(data, group_col = "group", cat_col = "category",
                             pct_col = "pct", palette, n_per_group = NULL,
                             legend_title = "", min_label_pct = 4,
                             bar_height = 0.35, base_size = 14,
                             border_col = "white",
                             show_guide_segments = TRUE,
                             n_wrap = FALSE) {
  # show_guide_segments: when FALSE, suppresses the dotted inter-row
  # connection segments that visually link adjacent group boundaries.
  # Callers with many groups + categories (e.g. ED04 per_series_rate
  # with 8 studies × 5 categories) find the resulting dotted scaffolding
  # visually taxing; passing FALSE returns a cleaner stacked-bar.
  # border_col controls the linewidth-0.5 border ggplot draws around each
  # bar segment. Default "white" reproduces the JAMA-style inter-segment
  # dividers used by Fig 1/2. Pass NA to drop the border entirely — useful
  # when the leftmost segment must sit flush against the y-axis line (e.g.
  # Fig 4 Panel E in a narrow composite cell, where the ~0.25-pt right
  # half of that border paints over the fill at x=0 and reads as a gap
  # between the axis and the first bar).

  data <- data %>%
    arrange(!!sym(group_col), !!sym(cat_col)) %>%
    group_by(!!sym(group_col)) %>%
    mutate(
      xmax = cumsum(!!sym(pct_col)),
      xmin = xmax - !!sym(pct_col),
      xmid = (xmin + xmax) / 2
    ) %>%
    ungroup()

  groups <- rev(unique(data[[group_col]]))
  y_map <- setNames(seq_along(groups), groups)
  # Force name-based lookup. If data[[group_col]] is a factor whose levels
  # differ from rev(unique(...)) order (e.g. ED Fig 4 panel a after the
  # v6.40 fct_relevel(ED_STUDY_ORDER)), indexing y_map[<factor>] would use
  # the factor's integer CODES positionally — drawing each group's bar at
  # another group's y-position while the y-axis labels stayed correct, so
  # labels and bars silently flipped. as.character() restores name lookup.
  data$y <- y_map[as.character(data[[group_col]])]

  if (!is.null(n_per_group)) {
    # v6.26 (2026-05-20): cohort name + n on the SAME line (was 2-line
    # newline-separated). Saves vertical space inside each bar row and
    # reads as "BCH (n = 110)" instead of "BCH" stacked above "(n = 110)".
    # v6.39 (2026-05-21): n_wrap = TRUE restores the 2-line wrap for
    # panels with long category names (e.g. ED08 Panel C rupture_categories
    # with strata like "KRAS-NonParietal").
    .sep <- if (isTRUE(n_wrap)) "\n" else " "
    y_labels <- paste0(names(y_map), .sep,
                       "(n = ", n_per_group[names(y_map)], ")")
  } else {
    y_labels <- names(y_map)
  }

  label_data <- data %>% filter(!!sym(pct_col) >= min_label_pct)

  p <- ggplot(data) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = y - bar_height, ymax = y + bar_height,
                  fill = !!sym(cat_col)),
              color = border_col,
              linewidth = if (is.na(border_col)) 0 else 0.5) +
    geom_text(data = label_data,
              aes(x = xmid, y = y, label = sprintf("%.1f", !!sym(pct_col))),
              size = base_size / 3.5, color = "black") +
    scale_fill_manual(values = palette, name = legend_title) +
    scale_y_continuous(breaks = unname(y_map), labels = y_labels) +
    scale_x_continuous(limits = c(0, 100.5), breaks = seq(0, 100, 10),
                       expand = c(0, 0)) +
    labs(x = "Patients (%)", y = NULL) +
    theme_avm(base_size = base_size) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )

  if (isTRUE(show_guide_segments) && length(groups) >= 2) {
    guide_segs <- tibble()
    for (i in seq_len(length(groups) - 1)) {
      g_top <- groups[i]; g_bot <- groups[i + 1]
      d_top <- data %>% filter(!!sym(group_col) == g_top) %>% arrange(!!sym(cat_col))
      d_bot <- data %>% filter(!!sym(group_col) == g_bot) %>% arrange(!!sym(cat_col))
      n_bounds <- min(nrow(d_top), nrow(d_bot)) - 1
      if (n_bounds > 0) {
        guide_segs <- bind_rows(guide_segs, tibble(
          x = d_top$xmax[seq_len(n_bounds)],
          xend = d_bot$xmax[seq_len(n_bounds)],
          y = y_map[g_top] - bar_height,
          yend = y_map[g_bot] + bar_height
        ))
      }
    }
    if (nrow(guide_segs) > 0) {
      p <- p + geom_segment(data = guide_segs,
        aes(x = x, xend = xend, y = y, yend = yend),
        linetype = "dotted", color = "grey50", linewidth = 0.5,
        inherit.aes = FALSE)
    }
  }

  p
}

# ── Dumbbell / lollipop prevalence plot ──────────────────────────────────────
# data: tibble with columns [feature, group, prevalence, n]
# palette: named color vector keyed by group
# size_range: range for point size mapping (min, max)
dumbbell_prevalence <- function(data, feature_col = "feature", group_col = "group",
                                prev_col = "prevalence", n_col = "n",
                                lo_col = NULL, hi_col = NULL,
                                palette, size_range = c(2.5, 5),
                                size_breaks = NULL, size_limits = NULL,
                                base_size = 14) {
  # `size_breaks` / `size_limits` let the caller force identical N-scale
  # legends across multiple dumbbells in a composite. With `guides =
  # "collect"`, patchwork only merges legends whose scales match
  # exactly — so when two dumbbell panels are stacked side-by-side
  # (e.g. Fig 3 sm_components + clinical_history), the producer should
  # compute the union of n values across both data frames and pass the
  # same `size_breaks` + `size_limits` to both calls.
  #
  # `lo_col` / `hi_col` are optional. When both are supplied, a
  # horizontal 95% CI segment is rendered behind each dot (same colour
  # mapping as the dot). If either is NULL no error bars are drawn —
  # preserves the no-CI behaviour for callers that don't pass them.
  size_scale <- if (is.null(size_breaks)) {
    scale_size_continuous(range = size_range, name = "N",
                          breaks = scales::breaks_pretty(3),
                          limits = size_limits)
  } else {
    scale_size_continuous(range = size_range, name = "N",
                          breaks = size_breaks,
                          limits = size_limits)
  }

  # Vertical dodge so the two genotype groups don't overlap within a
  # feature row. width=0.6 shifts each group by ±0.30 row-units above /
  # below the feature's y tick (tightened from ±0.45 so adjacent feature
  # rows don't visually bleed into each other while still keeping a
  # clear gap between the two within-row groups).
  dodge <- position_dodge(width = 0.6)

  # No connecting line: with the groups dodged into separate y bands
  # within each feature row, a "dumbbell" line at the un-dodged
  # centerline would float between (not through) the two dots and read
  # as an erroneous orphan. The vertical dodge alone pairs the groups
  # within each feature.
  p <- ggplot(data, aes(x = !!sym(prev_col), y = !!sym(feature_col)))

  if (!is.null(lo_col) && !is.null(hi_col)) {
    # Horizontal 95% CI bar, dodged with the dot so it sits on the same
    # vertical line as its group. Width is the cap span in y-axis
    # (categorical) units; 0.12 keeps the perpendicular caps inside the
    # dodged half-band. ggplot2 >=3.5 prefers `width` over `height`.
    p <- p + geom_errorbarh(
      aes(xmin = !!sym(lo_col), xmax = !!sym(hi_col),
          colour = !!sym(group_col)),
      position = dodge,
      width = 0.12, linewidth = 0.6, alpha = 0.85, show.legend = FALSE
    )
  }

  p <- p +
    geom_point(aes(colour = !!sym(group_col), size = !!sym(n_col)),
               position = dodge, alpha = 0.85) +
    scale_color_manual(values = palette, name = "Genotype") +
    size_scale +
    scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Prevalence (%)", y = NULL) +
    guides(color = guide_legend(order = 1, override.aes = list(size = 4)),
           size  = guide_legend(order = 2)) +
    theme_avm(base_size = base_size) +
    theme(panel.grid.major.y = element_blank())

  p
}

# ── Table-forest plot (JAMA style) ───────────────────────────────────────────
# Adapted from UKB NfL Cardiovascular project (Keller et al.)
# Layout: left-side text table | center forest plot | right-side estimate text
#
# Arguments:
#   data       — tibble with columns: label, est, lo, hi (and optional extras)
#   cols       — list of list(col=, header=) for left-side text columns
#   null_value — reference line (1 for OR/HR, 0 for beta)
#   log_scale  — if TRUE, x-axis is log-transformed
#   axis_ticks — numeric vector of axis tick values (original scale)
#   point_col  — fill color for data points
#   pooled_col — column flagging pooled row
#   pooled_flag— value identifying pooled rows
#   size_col   — column to map to point size (converted to numeric)
#   size_range — range for point size scaling
#   title, subtitle — plot title/subtitle
#   x_lab      — axis label below ticks
#   est_col_header — header for right-side estimate column
#   est_fmt    — sprintf format for estimate text
#   base_size  — text base size
table_forest <- function(data,
                         cols = list(),
                         null_value = 1,
                         log_scale = TRUE,
                         axis_ticks = NULL,
                         point_col = "#2166AC",
                         pooled_col = NULL,
                         pooled_flag = "pooled",
                         size_col = NULL,
                         size_range = c(2, 6),
                         title = NULL,
                         subtitle = NULL,
                         x_lab = "Odds Ratio (95% CI)",
                         est_col_header = NULL,
                         est_fmt = "%.2f (%.2f\u2013%.2f)",
                         # est_col_side: "left" (default; OR becomes the last
                         # column of the left-side text table, aligned with
                         # N / FDR P etc.) or "right" (OR text renders to the
                         # right of the forest bars). The repo convention is
                         # "left" — used by Fig 2C/E and every ED forest. The
                         # default was flipped from "right" → "left" after a
                         # round of ED panels rendered with the OR floating in
                         # otherwise-empty white space to the right of the
                         # forest bars; making "left" the default means new
                         # producers don't have to opt in to the convention.
                         est_col_side = c("left", "right"),
                         # est_col_two_rows: when TRUE, the estimate string
                         #   (e.g. "28.6 (11.3-52.2)") is rendered on a sub-row
                         #   directly beneath the predictor/N row instead of
                         #   competing for horizontal space alongside them.
                         #   Header gets the same treatment ("Rate (95% CI)"
                         #   sits below "Predictor / N"). Used in composites
                         #   where forest cells are squeezed below ~7 in wide
                         #   and the est-column text would otherwise collide.
                         est_col_two_rows = FALSE,
                         base_size = 14,
                         # ── ELab visual options (Keller/avm-survey style) ──
                         # group_col: column whose values define row groupings;
                         #   rows whose group_col value matches group_bg_flag
                         #   are rendered against a light-grey band (used for
                         #   the "Primary / Secondary" split in the avm-survey
                         #   forest plot). Default NULL = no band.
                         group_col = NULL,
                         group_bg_flag = NULL,
                         group_bg_fill = "#f5f5f5",
                         group_bg_alpha = 0.5,
                         # left_label / right_label: directional text placed
                         #   under the x-axis (e.g. "<- Favours control",
                         #   "Favours treatment ->"). Default NULL = omitted.
                         left_label = NULL,
                         right_label = NULL,
                         # show_bottom_rule: when FALSE, suppresses the
                         #   horizontal rule between the last data row and
                         #   the x-axis. Default TRUE preserves the JAMA-
                         #   table appearance. Fig 2's high-risk OR forest
                         #   uses FALSE so the panel reads as a clean forest
                         #   without an extra hairline above the OR axis.
                         show_bottom_rule = TRUE) {

  # ── Transform to log scale if requested ─────────────────────────────────────
  if (log_scale) {
    data <- data %>% mutate(
      x_est = log(.data$est),
      x_lo  = log(.data$lo),
      x_hi  = log(.data$hi)
    )
    null_x <- log(null_value)
    if (is.null(axis_ticks)) axis_ticks <- c(0.25, 0.5, 1, 2, 4)
    tick_x <- log(axis_ticks)
    tick_labels <- as.character(axis_ticks)
  } else {
    data <- data %>% mutate(x_est = .data$est, x_lo = .data$lo, x_hi = .data$hi)
    null_x <- null_value
    if (is.null(axis_ticks)) axis_ticks <- pretty(c(data$lo, data$hi), n = 5)
    tick_x <- axis_ticks
    tick_labels <- as.character(axis_ticks)
  }

  # ── Auto-extend axis so every point estimate lies within the forest range ──
  # Rationale: CI endpoints can be arbitrarily wide under sparse data and are
  # still best handled by clipping + arrowheads, but a *point estimate* that
  # falls outside the axis produces a hollow marker at the axis edge that
  # visually collides with the CI's arrowhead. Extending the axis to cover all
  # point estimates (with a small pad) keeps every box cleanly inside the
  # range. Ticks are auto-added in the user's tick cadence (doubling for log
  # OR axes, pretty breaks otherwise).
  est_finite <- data$x_est[is.finite(data$x_est)]
  if (length(est_finite) > 0) {
    est_lo <- min(est_finite)
    est_hi <- max(est_finite)
    if (log_scale) {
      # Extend by doubling/halving the outermost tick (i) until every
      # median fits inside the range AND (ii) one extra step beyond the
      # most extreme median, so a clipped CI arrowhead always has
      # horizontal headroom to render. Without the extra step, a median
      # that sits within half a log-step of the edge tick would clip its
      # arrow into the axis line. tick-step on a log-2 axis = log(2),
      # so we require >= half a tick (i.e. est_hi within log(sqrt(2)) of
      # the edge → extend).
      half_log_step <- log(2) / 2  # log(sqrt(2)) ≈ 0.347
      max_iter <- 20L
      i <- 0L
      while ((est_hi + half_log_step) > max(tick_x) + 1e-9 && i < max_iter) {
        next_val <- max(axis_ticks) * 2
        axis_ticks <- c(axis_ticks, next_val)
        tick_x <- c(tick_x, log(next_val))
        tick_labels <- c(tick_labels, formatC(next_val, format = "fg"))
        i <- i + 1L
      }
      i <- 0L
      while ((est_lo - half_log_step) < min(tick_x) - 1e-9 && i < max_iter) {
        next_val <- min(axis_ticks) / 2
        axis_ticks <- c(next_val, axis_ticks)
        tick_x <- c(log(next_val), tick_x)
        tick_labels <- c(formatC(next_val, format = "fg"), tick_labels)
        i <- i + 1L
      }
    } else {
      # Linear scale: widen pretty ticks to cover all point estimates,
      # PLUS one extra step of headroom beyond the extreme median so
      # clipped CI arrowheads always have horizontal room to render.
      step <- if (length(tick_x) >= 2) diff(range(tick_x)) / (length(tick_x) - 1) else 1
      pad  <- max(step, (max(tick_x) - min(tick_x)) * 0.05)
      need_extend_hi <- est_hi >= max(tick_x) - 1e-9
      need_extend_lo <- est_lo <= min(tick_x) + 1e-9
      if (need_extend_hi || need_extend_lo) {
        lo_lim <- if (need_extend_lo) est_lo - pad else min(tick_x)
        hi_lim <- if (need_extend_hi) est_hi + pad else max(tick_x)
        new_range <- range(c(tick_x, lo_lim, hi_lim))
        axis_ticks <- pretty(new_range, n = length(axis_ticks))
        tick_x <- axis_ticks
        tick_labels <- as.character(axis_ticks)
      }
    }
  }

  # ── Default estimate column header ──────────────────────────────────────────
  if (is.null(est_col_header)) {
    est_col_header <- if (log_scale) "OR (95% CI)" else "Estimate (95% CI)"
  }

  # ── Identify pooled row ─────────────────────────────────────────────────────
  if (!is.null(pooled_col) && pooled_col %in% names(data)) {
    data$is_pooled <- data[[pooled_col]] == pooled_flag
  } else {
    data$is_pooled <- FALSE
  }

  # ── Convert size_col to numeric ─────────────────────────────────────────────
  if (!is.null(size_col) && size_col %in% names(data)) {
    data$size_val <- as.numeric(as.character(data[[size_col]]))
  }

  # ── Y positions (bottom-up: first row at top) ──────────────────────────────
  # In two-row mode each study occupies a vertically taller block (top
  # sub-row + bottom sub-row + breathing room), so we multiply y_pos by
  # row_spacing > 1 to give every study its own visual band and stop the
  # bottom sub-row from colliding with the next study's top sub-row.
  n_rows      <- nrow(data)
  row_spacing <- if (isTRUE(est_col_two_rows)) 2.2 else 1.0
  data$y_pos  <- (rev(seq_len(n_rows)) - 1) * row_spacing
  if (any(data$is_pooled)) {
    pooled_y <- data$y_pos[data$is_pooled]
    # Lift non-pooled rows above the pooled row by 0.7 * row_spacing in
    # two-row mode (was 0.5) so the gap between BCH's CI sub-row and the
    # pooled "Pooled" label is visibly larger than the inter-study gap.
    pooled_lift <- if (isTRUE(est_col_two_rows)) 0.7 else 0.5
    data$y_pos[!data$is_pooled & data$y_pos > pooled_y] <-
      data$y_pos[!data$is_pooled & data$y_pos > pooled_y] + pooled_lift * row_spacing
  }

  # ── X layout: LEFT table | CENTER forest | (optional) RIGHT estimate text ─
  est_col_side <- match.arg(est_col_side)
  forest_min   <- min(tick_x)
  forest_max   <- max(tick_x)
  forest_range <- forest_max - forest_min

  # Left table region. When the OR column is on the left, one extra slot is
  # packed into the same strip, so widen the strip proportionally to keep
  # long Predictor labels clear of the right-aligned text columns (we saw
  # "Venous Outflow Stenosis" / N = "184" collide at the 2.0x width).
  #
  # Wide-axis guard: for log-scale forests that span >2 log10 decades (e.g.
  # Fig 4A, axis 0.25–1024), a forest_range-proportional table balloons to
  # ~11 log-units wide and visually pushes the forest into the right ~25% of
  # the cell. When log_scale && forest_range > 2, switch to a fixed per-slot
  # width so the table packs tight regardless of axis span. Fig 2 forests
  # stay on the original formula (their axes span ~1.2 log-units).
  #
  # Linear-scale forests (log_scale = FALSE) live in raw data units (e.g.
  # 0–100 for rate-%, 0–60 for years), so forest_range > 2 is essentially
  # always true and the fixed per-slot constant (3.5) is far too tight in
  # raw units — it crams every column into the same x position. For those
  # we use a fraction-of-forest_range formula tuned to leave the forest
  # ~60% of the panel and the table ~30% (linear_left/linear_right below).
  n_extra      <- length(cols)
  n_col_slots  <- if (est_col_side == "left") n_extra + 1 else n_extra
  # Per-column slot width (log10 units). Constants below were hand-tuned at
  # base_size 11; physical text width scales linearly with base_size, so we
  # rescale the multipliers by (base_size / 11) so the table packs the same
  # PROPORTION of the panel cell at any base_size.
  bs_scale     <- base_size / 11
  # In single-row mode the est column renders the full "%.2f (%.2f–%.2f)"
  # string (~17 chars) on one line vs. just "(%.2f–%.2f)" (~10 chars) in
  # two-row mode, so the est slot needs more horizontal room when est sits
  # LEFT of the forest. We don't widen table_width itself — total axis span
  # grows with it, so text width in axis units scales the same way and the
  # gap stays proportionally cramped. Instead, shift only the est slot
  # rightward by one full col_spacing below, after table_width is set.
  one_row_left <- !isTRUE(est_col_two_rows) && est_col_side == "left"
  table_width  <- if (log_scale && forest_range > 2) {
    (n_col_slots + 1) * 3.5 * bs_scale
  } else if (!log_scale) {
    # Linear-scale axes (raw rate-%, years, slope-units): take the larger of
    # (a) the log-scale wide-axis fixed-per-slot width, which dominates when
    #     forest_range is small (few raw units) — without it, a slope axis
    #     spanning −6…2 leaves only ~6 raw-units of table strip and the
    #     Predictor/N columns collide.
    # (b) a forest_range-proportional width, which dominates when forest_range
    #     is large (rate-% 0–100, age 0–60), so the forest doesn't get
    #     squeezed into a small fraction of the panel.
    max(
      (n_col_slots + 1) * 3.5,
      forest_range * (if (est_col_side == "left") 0.70 else 0.30)
    ) * bs_scale
  } else {
    forest_range * (if (est_col_side == "left") 3.0 else 2.0) * bs_scale
  }
  # Cap the table-to-forest gap so it doesn't balloon with wide axes.
  table_gap    <- min(forest_range * 0.25, 0.4)
  x_table_left <- forest_min - table_gap - table_width
  TX_LABEL     <- x_table_left

  if (n_col_slots > 0) {
    avail        <- table_width * 0.92
    # Predictor column gets 2.5 slot widths because its text is typically
    # the widest in the table (e.g. "Pooled (REML)", "KRAS G12V (|z| ≤ 2)"
    # or "Venous Outflow Stenosis"); the numeric / FDR / OR columns are
    # narrow. The extra padding keeps long predictor labels clear of the
    # right-aligned N column even at half-width composite cells.
    pred_slots   <- 2.5
    col_spacing  <- avail / (n_col_slots + pred_slots)
    slot_x       <- TX_LABEL + col_spacing * (pred_slots - 1 + seq_len(n_col_slots))
    # Single-row left-side est: shift the est slot one col_spacing right so
    # the "%.1f (%.1f–%.1f)" string clears the adjacent N column. The slot
    # ends ~1 col_spacing shy of forest_min (table_gap reserves ≤0.4 raw
    # units), so the shift stays inside the table strip.
    if (one_row_left && n_col_slots >= 2) {
      slot_x[n_col_slots] <- slot_x[n_col_slots] + col_spacing
    }
  } else {
    slot_x <- numeric(0)
  }

  if (est_col_side == "left") {
    tx_positions <- if (n_extra > 0) slot_x[seq_len(n_extra)] else numeric(0)
    TX_EST       <- slot_x[n_col_slots]
    x_right      <- forest_max + forest_range * 0.12
  } else {
    tx_positions <- slot_x
    TX_EST       <- forest_max + forest_range * 0.12
    x_right      <- TX_EST + forest_range * 0.85
  }

  # ── Key Y positions ────────────────────────────────────────────────────────
  # Two-row mode: each data row + the header gets a sub-row 0.5 axis-units
  # below its primary y-position, on which the est-column string is rendered.
  # We lift the header-band and drop the axis-band by sub_gap so the sub-rows
  # don't crowd their neighbours.
  sub_gap  <- if (isTRUE(est_col_two_rows)) 1.0 else 0
  # In two-row mode at composite scale, annotate text spans ~0.7 axis units
  # vertically. Push header_y up by hdr_pad and the rule down by rule_offset
  # so neither the "(95% CI)" sub-row nor the first data row's top text
  # crosses the rule line.
  # One-row mode pads bumped 2026-05-28 (ED Fig 4 panel c): each text row is
  # ~0.5 axis-units tall at base_size 8 in a 2.85-in panel, so the previous
  # 0.6 / 0.35 hdr/rule pads put the top rule right on the first row's text
  # top edge. Lift hdr_pad to 1.0 → rule_y sits ~0.65 au above the top row,
  # giving a clean separation that matches the two-row case proportionally.
  hdr_pad     <- if (isTRUE(est_col_two_rows)) 1.5 else 1.0
  rule_offset <- if (isTRUE(est_col_two_rows)) 0.6 else 0.35
  header_y <- max(data$y_pos) + hdr_pad + sub_gap
  axis_y   <- min(data$y_pos) - 1.0 - sub_gap
  rule_y   <- header_y - rule_offset - sub_gap

  # ── Build plot ──────────────────────────────────────────────────────────────
  dr <- data %>% filter(!is_pooled)
  pl <- data %>% filter(is_pooled)

  p <- ggplot()

  # ── ELab-style group shading band (opt-in via group_col + group_bg_flag) ─
  # Rendered first so all other geoms sit on top. Typical use: flag the
  # secondary-outcome block in a multi-row forest so readers see the
  # primary/secondary boundary at a glance (see avm-survey Figure 3).
  if (!is.null(group_col) && !is.null(group_bg_flag) &&
      group_col %in% names(data)) {
    band_rows <- data[data[[group_col]] == group_bg_flag, , drop = FALSE]
    if (nrow(band_rows) > 0) {
      p <- p + annotate(
        "rect",
        xmin = TX_LABEL - 0.3, xmax = x_right,
        ymin = min(band_rows$y_pos) - 0.4,
        ymax = max(band_rows$y_pos) + 0.4,
        fill = group_bg_fill, alpha = group_bg_alpha
      )
    }
  }

  # 2026-05-19 NM pass: dashed (not dotted) reference line at the null
  # value, slightly thicker and darker so it's actually visible at
  # native composite scale. Spans from the axis line up to the top of
  # the data rows (just below the header rule).
  p <- p +
    annotate("segment", x = null_x, xend = null_x,
             y = axis_y + 0.1, yend = rule_y,
             linetype = "dashed", colour = "grey45", linewidth = 0.35)

  # Axis ticks + labels. tick_text_y_pad is the gap between the axis line
  # and the tick label centre. Text bounding boxes span ~0.7 axis units in
  # two-row mode at composite scale, so a 0.08-unit pad isn't enough — the
  # tick labels touch the bottom rule. Use 0.55 in two-row mode.
  # 2026-05-28 (ED Fig 4 panel c): tick text top edge sits ~0.25 au above
  # its centre at this scale, so a 0.08 pad left the labels straddling the
  # axis line (axis_y + 0.10). Drop the labels by 0.40 au so the gap to the
  # axis line is ~0.25 au — visible but compact.
  tick_text_y_pad <- if (isTRUE(est_col_two_rows)) 0.55 else 0.40
  # 2026-05-19 NM pass: pad each tick label with a thin space on either
  # side so adjacent log-spaced labels (e.g. "0.5", "1", "2") read as
  # distinct words instead of running together. U+2009 THIN SPACE is
  # narrower than a normal space; one on each side adds ~0.15 axis
  # units of visual gap at base_size=14 without offsetting the centre.
  tick_labels_padded <- paste0("  ", tick_labels, "  ")
  for (i in seq_along(tick_x)) {
    p <- p +
      annotate("segment", x = tick_x[i], xend = tick_x[i],
               y = axis_y + 0.05, yend = axis_y + 0.18,
               colour = "grey60", linewidth = 0.3) +
      annotate("text", x = tick_x[i], y = axis_y - tick_text_y_pad,
               label = tick_labels_padded[i], size = base_size / 4,
               colour = "black", hjust = 0.5)
  }
  p <- p + annotate("segment", x = forest_min, xend = forest_max,
                     y = axis_y + 0.1, yend = axis_y + 0.1,
                     colour = "grey60", linewidth = 0.4)
  # x-axis title sits below the tick labels. 1.4 units below axis_y in
  # two-row mode (was 0.85) to leave a clear gap from the tick labels.
  # 2026-05-28: paired with the tick_text_y_pad bump above; the title needs
  # to drop with the ticks or it lands on top of them.
  axis_title_y_pad <- if (isTRUE(est_col_two_rows)) 1.4 else 1.30
  p <- p + annotate("text", x = mean(c(forest_min, forest_max)),
                     y = axis_y - axis_title_y_pad, label = x_lab,
                     size = base_size / 3.8, colour = "black", hjust = 0.5)

  # ── ELab-style directional labels (opt-in via left_label / right_label) ──
  # Placed beneath the x-axis on the null-line sides, e.g.
  # "<- Favours no therapy" / "Favours therapy ->". Y-position is tracked in
  # `dir_label_y` so the final scale_y_continuous() call (below) can extend
  # the lower bound to keep the labels inside the plot region.
  dir_label_y <- NA_real_
  if (!is.null(left_label) || !is.null(right_label)) {
    dir_label_y <- axis_y - 0.85
    if (!is.null(left_label)) {
      p <- p + annotate(
        "text",
        x = mean(c(forest_min, null_x)), y = dir_label_y,
        label = left_label, size = base_size / 4.5,
        colour = "black", hjust = 0.5
      )
    }
    if (!is.null(right_label)) {
      p <- p + annotate(
        "text",
        x = mean(c(null_x, forest_max)), y = dir_label_y,
        label = right_label, size = base_size / 4.5,
        colour = "black", hjust = 0.5
      )
    }
  }

  # ── Data row CIs + sized points ───────────────────────────────────────────
  # Clipping behavior (standard forest-plot convention, e.g. Cochrane/metafor):
  #   • CI endpoints outside [forest_min, forest_max] are clipped to the axis
  #     and marked with an arrowhead on the clipped side.
  #   • Point estimates outside the axis are clamped to the axis edge and
  #     rendered as an open (hollow) marker to indicate the estimate itself
  #     lies beyond the displayed range. The numeric OR (right-side column)
  #     still shows the true value.
  # This preserves faithful representation without silently dropping markers
  # or misaligning the point with its error bar.
  arrow_style <- grid::arrow(length = unit(0.10, "inches"),
                             ends = "last", type = "closed")

  draw_forest_row <- function(p, rows, err_height, err_lw, pt_size_default,
                              pt_shape = 15, pt_col, pt_fill = NA,
                              pt_border = NULL, use_size_col = FALSE) {
    if (nrow(rows) == 0) return(p)
    rows <- rows %>% mutate(
      x_lo_clip = pmax(x_lo, forest_min),
      x_hi_clip = pmin(x_hi, forest_max),
      x_est_clip = pmax(pmin(x_est, forest_max), forest_min),
      lo_clipped = x_lo < forest_min,
      hi_clipped = x_hi > forest_max,
      est_clipped = x_est < forest_min | x_est > forest_max
    )

    # Base horizontal bar (between the two clipped endpoints)
    p <- p + geom_segment(
      data = rows,
      aes(y = y_pos, yend = y_pos, x = x_lo_clip, xend = x_hi_clip),
      linewidth = err_lw, colour = pt_col,
      inherit.aes = FALSE
    )
    # Whisker caps on unclipped sides
    cap_lo <- rows %>% filter(!lo_clipped)
    cap_hi <- rows %>% filter(!hi_clipped)
    if (nrow(cap_lo) > 0) {
      p <- p + geom_segment(
        data = cap_lo,
        aes(x = x_lo_clip, xend = x_lo_clip,
            y = y_pos - err_height / 2, yend = y_pos + err_height / 2),
        linewidth = err_lw, colour = pt_col,
        inherit.aes = FALSE
      )
    }
    if (nrow(cap_hi) > 0) {
      p <- p + geom_segment(
        data = cap_hi,
        aes(x = x_hi_clip, xend = x_hi_clip,
            y = y_pos - err_height / 2, yend = y_pos + err_height / 2),
        linewidth = err_lw, colour = pt_col,
        inherit.aes = FALSE
      )
    }
    # Arrowheads on clipped sides (point inward → outward along axis)
    arr_lo <- rows %>% filter(lo_clipped)
    arr_hi <- rows %>% filter(hi_clipped)
    if (nrow(arr_lo) > 0) {
      p <- p + geom_segment(
        data = arr_lo,
        aes(y = y_pos, yend = y_pos,
            x = forest_min + forest_range * 0.04, xend = forest_min),
        linewidth = err_lw, colour = pt_col, arrow = arrow_style,
        inherit.aes = FALSE
      )
    }
    if (nrow(arr_hi) > 0) {
      p <- p + geom_segment(
        data = arr_hi,
        aes(y = y_pos, yend = y_pos,
            x = forest_max - forest_range * 0.04, xend = forest_max),
        linewidth = err_lw, colour = pt_col, arrow = arrow_style,
        inherit.aes = FALSE
      )
    }

    # Points (clamped). Clamped points use hollow marker (shape + 0 = filled,
    # so we override with open variants when clipped).
    pts_inside  <- rows %>% filter(!est_clipped)
    pts_clipped <- rows %>% filter(est_clipped)

    if (use_size_col && "size_val" %in% names(rows)) {
      if (nrow(pts_inside) > 0) {
        p <- p + geom_point(
          data = pts_inside,
          aes(y = y_pos, x = x_est_clip, size = size_val),
          shape = pt_shape, colour = pt_col, fill = pt_col,
          show.legend = FALSE, inherit.aes = FALSE
        )
      }
      if (nrow(pts_clipped) > 0) {
        # Hollow marker to flag: point estimate is beyond axis range
        p <- p + geom_point(
          data = pts_clipped,
          aes(y = y_pos, x = x_est_clip, size = size_val),
          shape = 0, colour = pt_col, stroke = 1.1,
          show.legend = FALSE, inherit.aes = FALSE
        )
      }
      p <- p + scale_size_continuous(range = size_range)
    } else {
      if (nrow(pts_inside) > 0) {
        p <- p + geom_point(
          data = pts_inside,
          aes(y = y_pos, x = x_est_clip),
          size = pt_size_default, shape = pt_shape,
          colour = ifelse(is.null(pt_border), pt_col, pt_border),
          fill = if (!is.na(pt_fill)) pt_fill else pt_col,
          inherit.aes = FALSE
        )
      }
      if (nrow(pts_clipped) > 0) {
        p <- p + geom_point(
          data = pts_clipped,
          aes(y = y_pos, x = x_est_clip),
          size = pt_size_default, shape = 0,
          colour = pt_col, stroke = 1.1,
          inherit.aes = FALSE
        )
      }
    }
    p
  }

  # In two-row mode, the dot/CI bar should sit vertically centered between
  # the predictor+N row (top sub-row) and the rate-value row (bottom sub-row).
  # Shift the y_pos in the data frame passed to dot rendering by -sub_gap/2;
  # text rendering further down still uses the original data$y_pos.
  dot_y_offset <- if (isTRUE(est_col_two_rows)) -sub_gap / 2 else 0
  dr_dot <- if (dot_y_offset != 0) dr %>% mutate(y_pos = y_pos + dot_y_offset) else dr

  # Per-study rows
  p <- draw_forest_row(
    p, dr_dot,
    err_height = 0.18, err_lw = 0.8,
    pt_size_default = 3.5, pt_shape = 15, pt_col = point_col,
    use_size_col = !is.null(size_col)
  )

  # Pooled row (diamond) — keep diamond; if clipped (unlikely for pooled),
  # fall back to hollow diamond.
  if (nrow(pl) > 0) {
    pl <- pl %>% mutate(
      x_lo_clip = pmax(x_lo, forest_min),
      x_hi_clip = pmin(x_hi, forest_max),
      x_est_clip = pmax(pmin(x_est, forest_max), forest_min),
      lo_clipped = x_lo < forest_min,
      hi_clipped = x_hi > forest_max,
      est_clipped = x_est < forest_min | x_est > forest_max,
      y_pos = y_pos + dot_y_offset
    )
    p <- p +
      geom_segment(
        data = pl,
        aes(y = y_pos, yend = y_pos, x = x_lo_clip, xend = x_hi_clip),
        linewidth = 1.0, colour = "black", inherit.aes = FALSE
      )
    # Pooled arrows if clipped
    arr_lo_p <- pl %>% filter(lo_clipped)
    arr_hi_p <- pl %>% filter(hi_clipped)
    if (nrow(arr_lo_p) > 0) {
      p <- p + geom_segment(
        data = arr_lo_p,
        aes(y = y_pos, yend = y_pos,
            x = forest_min + forest_range * 0.04, xend = forest_min),
        linewidth = 1.0, colour = "black", arrow = arrow_style,
        inherit.aes = FALSE
      )
    }
    if (nrow(arr_hi_p) > 0) {
      p <- p + geom_segment(
        data = arr_hi_p,
        aes(y = y_pos, yend = y_pos,
            x = forest_max - forest_range * 0.04, xend = forest_max),
        linewidth = 1.0, colour = "black", arrow = arrow_style,
        inherit.aes = FALSE
      )
    }
    pooled_inside  <- pl %>% filter(!est_clipped)
    pooled_clipped <- pl %>% filter(est_clipped)
    if (nrow(pooled_inside) > 0) {
      p <- p + geom_point(
        data = pooled_inside,
        aes(y = y_pos, x = x_est_clip),
        size = 5, shape = 23, fill = point_col, colour = "black",
        inherit.aes = FALSE
      )
    }
    if (nrow(pooled_clipped) > 0) {
      p <- p + geom_point(
        data = pooled_clipped,
        aes(y = y_pos, x = x_est_clip),
        size = 5, shape = 5, colour = "black", stroke = 1.1,
        inherit.aes = FALSE
      )
    }
  }

  # ── Left table: label column ───────────────────────────────────────────────
  # In two-row mode, parenthetical qualifiers in the label (typically the
  # pooled row's "Pooled (REML)") would push the bold label into the N
  # column. Split the label at the first " (" boundary and render the head
  # ("Pooled") on the top sub-row at full weight, with the parenthetical
  # ("(REML)") on the bottom sub-row at lighter weight, mirroring the est
  # column's point/CI split.
  label_face <- ifelse(data$is_pooled, "bold", "plain")
  display_label <- data$label
  sub_label     <- rep("", nrow(data))
  if (isTRUE(est_col_two_rows)) {
    has_parens   <- grepl(" \\(.*\\)$", data$label, perl = TRUE)
    display_label[has_parens] <- sub(" \\(.*\\)$", "", data$label[has_parens],
                                     perl = TRUE)
    sub_label[has_parens]     <- regmatches(data$label[has_parens],
                                            regexpr(" \\(.*\\)$",
                                                    data$label[has_parens],
                                                    perl = TRUE))
    sub_label[has_parens]     <- sub("^ ", "", sub_label[has_parens])
  }
  p <- p +
    annotate("text", x = TX_LABEL, y = data$y_pos,
             label = display_label, hjust = 0,
             size = base_size / 3.5, fontface = label_face, colour = "black")
  if (any(nzchar(sub_label))) {
    sub_rows <- which(nzchar(sub_label))
    p <- p +
      annotate("text", x = TX_LABEL, y = data$y_pos[sub_rows] - sub_gap,
               label = sub_label[sub_rows], hjust = 0,
               size = base_size / 3.8, fontface = "plain", colour = "black")
  }

  # ── Left table: additional columns ─────────────────────────────────────────
  for (j in seq_along(cols)) {
    col_data <- as.character(data[[cols[[j]]$col]])
    col_data[is.na(col_data) | col_data == "NA"] <- "\u2014"
    p <- p +
      annotate("text", x = tx_positions[j], y = data$y_pos,
               label = col_data, hjust = 1,
               size = base_size / 4, colour = "black")
  }

  # ── Estimate (95% CI) text ─────────────────────────────────────────────────
  # hjust: 0 (left-aligned) when OR sits to the right of the forest; 1
  # (right-aligned) when OR is the last column of the left table so its
  # text ends flush with neighbouring table columns.
  #
  # In two-row mode the est-column is split: the point estimate ("28.6")
  # stays on the predictor/N row and only the CI portion ("(11.3-52.2)")
  # drops to the sub-row. We split est_fmt at the first " (" or " [" so the
  # producer can keep its existing one-string format.
  est_h    <- if (est_col_side == "left") 1 else 0
  .split_fmt <- function(s) {
    m <- regmatches(s, regexpr(" \\(.*\\)$| \\[.*\\]$", s, perl = TRUE))
    if (length(m) == 1L && nzchar(m))
      list(point = sub(" \\(.*\\)$| \\[.*\\]$", "", s, perl = TRUE),
           ci    = sub("^ ", "", m))
    else list(point = s, ci = NULL)
  }
  if (isTRUE(est_col_two_rows)) {
    est_split   <- .split_fmt(est_fmt)
    point_text  <- sprintf(est_split$point, data$est)
    ci_text     <- if (!is.null(est_split$ci))
      sprintf(est_split$ci, data$lo, data$hi) else rep("", nrow(data))
    p <- p +
      annotate("text", x = TX_EST, y = data$y_pos,
               label = point_text, hjust = est_h,
               size = base_size / 3.8, colour = "black",
               fontface = ifelse(data$is_pooled, "bold", "plain"))
    if (any(nzchar(ci_text))) {
      # Center the CI text under the rate value's visual center, rather
      # than right-aligning it at TX_EST (which makes the wider CI string
      # extend further left than the rate it belongs to). Estimate the
      # rate's half-width using the longest formatted point and shift x.
      #
      # The shift is expressed as a fraction of the inter-column spacing
      # rather than in raw axis units. The previous formula
      # (max_pt_chars * 0.45 * base_size/11) worked when the axis was
      # linear (rate %, age years, raw counts) — there 0.45 axis units
      # ≈ one character of rendered text. Under a log10 axis (OR / HR
      # forests with span ≤ 2 log10-units) the axis units are an order
      # of magnitude smaller per char, so the same formula shifted the
      # CI text out of the OR column entirely (visually landing under
      # the N column instead). Tying the shift to col_spacing keeps the
      # CI under the point estimate regardless of axis scale.
      max_pt_chars <- max(nchar(point_text), na.rm = TRUE)
      half_pt_w    <- if (n_col_slots > 0) {
        # ~7% of inter-column gap per character, capped at half a column
        # so the CI never crosses into the previous text column.
        min(max_pt_chars * col_spacing * 0.07, col_spacing * 0.45)
      } else {
        max_pt_chars * 0.45 * (base_size / 11)
      }
      ci_x         <- if (est_col_side == "left") TX_EST - half_pt_w
                      else TX_EST + half_pt_w
      p <- p +
        annotate("text", x = ci_x, y = data$y_pos - sub_gap,
                 label = ci_text, hjust = 0.5,
                 size = base_size / 3.8, colour = "black",
                 fontface = "plain")
    }
  } else {
    or_text <- sprintf(est_fmt, data$est, data$lo, data$hi)
    p <- p +
      annotate("text", x = TX_EST, y = data$y_pos,
               label = or_text, hjust = est_h,
               size = base_size / 3.8, colour = "black",
               fontface = ifelse(data$is_pooled, "bold", "plain"))
  }

  # ── Column headers ──────────────────────────────────────────────────────────
  p <- p +
    annotate("text", x = TX_LABEL, y = header_y, label = "Predictor",
             hjust = 0, size = base_size / 3.8, fontface = "bold")
  for (j in seq_along(cols)) {
    p <- p +
      annotate("text", x = tx_positions[j], y = header_y,
               label = cols[[j]]$header, hjust = 1,
               size = base_size / 3.8, fontface = "bold")
  }
  if (isTRUE(est_col_two_rows)) {
    # Drop the "(95% CI)" sub-line on the header — the parenthetical CI
    # format is self-evident from the data rows below, and removing the
    # sub-header simplifies the top of the plot.
    hdr_split <- .split_fmt(est_col_header)
    p <- p +
      annotate("text", x = TX_EST, y = header_y,
               label = hdr_split$point, hjust = est_h,
               size = base_size / 3.8, fontface = "bold")
  } else {
    p <- p +
      annotate("text", x = TX_EST, y = header_y,
               label = est_col_header, hjust = est_h,
               size = base_size / 3.8, fontface = "bold")
  }

  # ── Horizontal rules ───────────────────────────────────────────────────────
  x_left <- TX_LABEL - 0.05
  p <- p +
    annotate("segment", x = x_left, xend = x_right,
             y = rule_y, yend = rule_y,
             colour = "grey50", linewidth = 0.4)
  # Bottom rule (between last data row and x-axis) is opt-out via the
  # show_bottom_rule arg so panels that want a cleaner forest without
  # the extra hairline above the axis can suppress it.
  if (isTRUE(show_bottom_rule)) {
    p <- p +
      annotate("segment", x = x_left, xend = x_right,
               y = axis_y + 0.35, yend = axis_y + 0.35,
               colour = "grey50", linewidth = 0.4)
  }

  # ── Title / subtitle ───────────────────────────────────────────────────────
  # In two-row mode, annotate text spans ~0.7 axis-units, so the previous
  # 0.6-unit title-to-header gap had the title baseline overlapping the
  # Predictor header. Push title up by 1.4 units in two-row mode.
  if (!is.null(title)) {
    title_pad <- if (isTRUE(est_col_two_rows)) 1.4 else 0.6
    title_y <- header_y + ifelse(!is.null(subtitle), 1.0, title_pad)
    p <- p + annotate("text", x = TX_LABEL, y = title_y,
                       label = title, hjust = 0,
                       size = base_size / 2.8, fontface = "bold")
  }
  if (!is.null(subtitle)) {
    p <- p + annotate("text", x = TX_LABEL, y = header_y + 0.55,
                       label = subtitle, hjust = 0,
                       size = base_size / 3.4, colour = "black")
  }

  # ── Final theme ─────────────────────────────────────────────────────────────
  # If directional labels were drawn, extend the y-limit so they don't clip.
  # 2026-05-28: expanded from 1.05 -> 1.70 to keep the dropped axis title
  # (axis_y - 1.30) inside the plot area.
  y_lo_pad <- if (isTRUE(est_col_two_rows)) 1.8 else 1.70
  y_lo <- if (!is.na(dir_label_y)) dir_label_y - 0.25 else axis_y - y_lo_pad
  # Top buffer above the header row. 1.6 units reserves room for the
  # optional title/subtitle; when neither is passed we collapse the
  # buffer so short (2-3 row) forests fill their allocated vertical
  # space instead of floating in the middle of the panel.
  top_buffer <- if (!is.null(title) || !is.null(subtitle)) {
    if (isTRUE(est_col_two_rows)) 2.2 else 1.6
  } else 0.25
  p <- p +
    scale_x_continuous(limits = c(x_left, x_right), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_lo, header_y + top_buffer)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size) +
    theme(plot.margin = margin(8, 14, 8, 14))

  p
}

# ── Meta-analysis forest plot wrapper ─────────────────────────────────────────
# For per-study meta-analysis with pooled diamond row; dot size ∝ N.
# data must have: label, est, lo, hi, type ("study"/"pooled")
# size_col defaults to "n_label" (character N converted to numeric internally)
table_forest_meta <- function(data,
                              cols = list(),
                              null_value = 1,
                              log_scale = TRUE,
                              axis_ticks = NULL,
                              point_col = "#2166AC",
                              size_range = c(2, 7),
                              title = NULL,
                              subtitle = NULL,
                              x_lab = "Odds Ratio (95% CI)",
                              est_col_header = NULL,
                              est_fmt = "%.2f (%.2f\u2013%.2f)",
                              est_col_side = c("left", "right"),
                              est_col_two_rows = FALSE,
                              base_size = 14,
                              group_col = NULL, group_bg_flag = NULL,
                              group_bg_fill = "#f5f5f5",
                              group_bg_alpha = 0.5,
                              left_label = NULL, right_label = NULL,
                              show_bottom_rule = TRUE) {
  # v6.25 (2026-05-20): show_bottom_rule forwarded to table_forest so
  # the meta forest can drop the extra horizontal rule below the last
  # study row (matches the table_forest_general signature; previously
  # only the _general wrapper exposed this knob).
  table_forest(data, cols = cols, null_value = null_value,
               log_scale = log_scale, axis_ticks = axis_ticks,
               point_col = point_col,
               pooled_col = "type", pooled_flag = "pooled",
               size_col = "n_label", size_range = size_range,
               title = title, subtitle = subtitle,
               x_lab = x_lab, est_col_header = est_col_header,
               est_fmt = est_fmt, est_col_side = est_col_side,
               est_col_two_rows = est_col_two_rows,
               base_size = base_size,
               group_col = group_col, group_bg_flag = group_bg_flag,
               group_bg_fill = group_bg_fill,
               group_bg_alpha = group_bg_alpha,
               left_label = left_label, right_label = right_label,
               show_bottom_rule = show_bottom_rule)
}

# ── General OR forest plot wrapper ────────────────────────────────────────────
# For predictor-level OR forests (no pooled row); dot size ∝ N.
# data must have: label, est, lo, hi, n_label
table_forest_general <- function(data,
                                 cols = list(),
                                 null_value = 1,
                                 log_scale = TRUE,
                                 axis_ticks = NULL,
                                 point_col = "#2166AC",
                                 size_range = c(2, 7),
                                 title = NULL,
                                 subtitle = NULL,
                                 x_lab = "Odds Ratio (95% CI)",
                                 est_col_header = NULL,
                                 est_fmt = "%.2f (%.2f\u2013%.2f)",
                                 est_col_side = c("left", "right"),
                                 est_col_two_rows = FALSE,
                                 base_size = 14,
                                 group_col = NULL, group_bg_flag = NULL,
                                 group_bg_fill = "#f5f5f5",
                                 group_bg_alpha = 0.5,
                                 left_label = NULL, right_label = NULL,
                                 show_bottom_rule = TRUE) {
  table_forest(data, cols = cols, null_value = null_value,
               log_scale = log_scale, axis_ticks = axis_ticks,
               point_col = point_col,
               pooled_col = NULL, pooled_flag = NULL,
               size_col = "n_label", size_range = size_range,
               title = title, subtitle = subtitle,
               x_lab = x_lab, est_col_header = est_col_header,
               est_fmt = est_fmt, est_col_side = est_col_side,
               est_col_two_rows = est_col_two_rows,
               base_size = base_size,
               group_col = group_col, group_bg_flag = group_bg_flag,
               group_bg_fill = group_bg_fill,
               group_bg_alpha = group_bg_alpha,
               left_label = left_label, right_label = right_label,
               show_bottom_rule = show_bottom_rule)
}

# \u2500\u2500 Compact-cell forest wrappers \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# Convenience wrappers for the half-/third-width composite cell case.
# Setting est_col_two_rows = TRUE in either of these flips on the full bundle
# of layout adjustments tuned in 2026-04-25 for ED Fig 4:
#   * Each study row splits into a top sub-row (Predictor + N + point
#     estimate) and a bottom sub-row (CI in parentheses, centered under
#     the point estimate). Forest dot/CI bar sits vertically between them.
#   * Header reads as "Predictor / N / Rate" (or Mean / OR / etc.); the
#     "(95% CI)" sub-line is dropped \u2014 the CI format is self-evident from
#     the data rows below.
#   * Predictor labels with a " (...)" suffix (e.g. "Pooled (REML)") split
#     so the parenthetical sits on the bottom sub-row at lighter weight.
#   * Row spacing 1.6\u00d7 wider; header / rule / axis-title pads scaled up so
#     no two text rows overlap at composite scale.
#
# Use these from any producer whose forest panel will be composed into a
# narrow cell:
#   p <- forest_compact_meta(data, x_lab = "...", title = "...")
#   p <- forest_compact_general(data, x_lab = "...", title = "...")
#
# Standalone callers that don't hit a cell-width crunch should keep using
# table_forest_meta / table_forest_general (default est_col_two_rows = FALSE).

forest_compact_meta <- function(...) {
  table_forest_meta(..., est_col_two_rows = TRUE)
}

forest_compact_general <- function(...) {
  table_forest_general(..., est_col_two_rows = TRUE)
}

# ── Single-cell UMAP theme ───────────────────────────────────────────────────
# Matches theme_avm() for typography (bold title, grey-40 subtitle) but
# strips everything that isn't meaningful for a UMAP: tick labels, gridlines,
# axis text. Axis arrows in the corner signal "UMAP 1 / UMAP 2" without
# consuming a full tick bar. base_size follows theme_avm's default so
# single- and multi-panel composites stay typographically consistent.
theme_umap_avm <- function(base_size = 14) {
  # Axis-blanking theme for UMAPs. We DO NOT draw the axis arrows here
  # any longer — they're rendered as data-coordinate segments inside
  # `umap_origin_labels()` so their geometry stays consistent across
  # panels of differing patchwork-cell aspect (theme-drawn arrows scale
  # with the panel area, which produced size/position drift between
  # paired UMAPs in the Fig 4 composite).
  #
  # `aspect.ratio = 1` forces every UMAP panel to render square in
  # display units regardless of the patchwork cell aspect — so two
  # paired UMAPs (e.g. atlas + variant overlay) always look the same
  # shape. Cells whose aspect doesn't match get whitespace padding
  # rather than a stretched plot.
  theme_avm(base_size = base_size) +
    theme(
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.title        = element_blank(),
      axis.line         = element_blank(),
      panel.grid        = element_blank(),
      aspect.ratio      = 1
    )
}

# ── umap_origin_labels ──────────────────────────────────────────────────────
# Draws the UMAP axis arrows AND their "UMAP 1" / "UMAP 2" labels as
# data-coordinate annotate() layers anchored at the bottom-left corner
# of the plot. Use on plots already themed with theme_umap_avm() (which
# blanks ggplot's own axis line + titles).
#
# Why data-coord arrows instead of theme axis-line arrows: the theme's
# `axis.line` runs the full length of the panel's axis, so its arrow
# tip lands at a different physical position whenever patchwork hands
# the panel a different-shape cell. By drawing both arrows as
# `annotate("segment", arrow = ...)` with length = `arrow_frac` × data
# range, two panels sharing the same UMAP coordinate space (e.g. atlas
# vs variant overlay) render arrows of identical data-extent, so they
# match across the composite.
#
# `arrow_frac` defaults to 0.18 — i.e. each axis arrow spans ~18% of
# the data range. Tweak per-call for unusually compact data (smaller
# UMAPs may want 0.25) or to taste.
umap_origin_labels <- function(p, x_lab = "UMAP 1", y_lab = "UMAP 2",
                               base_size = 14, colour = "grey40",
                               arrow_frac = 0.18,
                               arrow_head = 0.18) {
  bld <- ggplot_build(p)
  xr  <- bld$layout$panel_params[[1]]$x.range
  yr  <- bld$layout$panel_params[[1]]$y.range
  xs  <- diff(xr); ys <- diff(yr)
  ax_len_x <- xs * arrow_frac
  ax_len_y <- ys * arrow_frac
  text_size <- (base_size - 2) / ggplot2::.pt
  ox <- xr[1]
  oy <- yr[1]

  p +
    # X-axis arrow: tail at (ox, oy), head at (ox + ax_len_x, oy)
    annotate("segment",
             x = ox, xend = ox + ax_len_x,
             y = oy, yend = oy,
             arrow = grid::arrow(length = unit(arrow_head, "cm"),
                                 ends = "last", type = "closed"),
             colour = colour, linewidth = 0.5) +
    # Y-axis arrow: tail at (ox, oy), head at (ox, oy + ax_len_y)
    annotate("segment",
             x = ox, xend = ox,
             y = oy, yend = oy + ax_len_y,
             arrow = grid::arrow(length = unit(arrow_head, "cm"),
                                 ends = "last", type = "closed"),
             colour = colour, linewidth = 0.5) +
    # Labels just below / left of the arrow tails
    annotate("text",
             x = ox, y = oy - ys * 0.04,
             label = x_lab, hjust = 0, vjust = 1,
             size = text_size, colour = colour) +
    annotate("text",
             x = ox - xs * 0.04, y = oy,
             label = y_lab, hjust = 0, vjust = 0, angle = 90,
             size = text_size, colour = colour) +
    # Re-pass the build-resolved x/y range to coord_cartesian. Reading
    # via ggplot_build means we pick up any upstream coord_cartesian
    # (e.g. xy_lims for per-patient UMAP grids) and re-apply it here
    # with clip = "off"; otherwise this last coord call would wipe
    # callers' custom xlim/ylim.
    coord_cartesian(xlim = xr, ylim = yr, clip = "off") +
    theme(plot.margin = margin(6, 8, 12, 14))
}

# ── umap_background_highlight ───────────────────────────────────────────────
# Two-layer UMAP: faint grey backdrop of all cells + coloured highlight of
# a subset. Used by Fig 4d (variant overlay on atlas) and by the per-patient
# facet UMAPs in the Extended Data. Inputs are two tibbles with columns
# `umap_1`, `umap_2`, plus `color_col` on the highlight tibble.
#
# Arguments:
#   background_df — data frame of non-highlighted cells (rendered grey)
#   highlight_df  — data frame of highlighted cells (rendered coloured)
#   color_col     — column in highlight_df to map to colour
#   palette       — named colour vector (defaults to VARIANT_SC_COLORS)
#   bg_color      — background point colour (default grey92)
#   bg_alpha      — background alpha (default 0.35)
#   bg_size       — background point size (default 0.25)
#   fg_size       — foreground point size (default 1.4)
#   fg_alpha      — foreground alpha (default 0.9)
#   facet_col     — optional column to facet by (NULL = single panel)
#   legend_title  — legend heading
umap_background_highlight <- function(background_df, highlight_df,
                                      color_col,
                                      palette = VARIANT_SC_COLORS,
                                      bg_color = "grey92",
                                      bg_alpha = 0.35,
                                      bg_size  = 0.25,
                                      fg_size  = 1.4,
                                      fg_alpha = 0.9,
                                      facet_col = NULL,
                                      legend_title = NULL,
                                      base_size = 14) {
  p <- ggplot() +
    geom_point(data = background_df,
               aes(x = umap_1, y = umap_2),
               colour = bg_color, alpha = bg_alpha, size = bg_size,
               stroke = 0) +
    geom_point(data = highlight_df,
               aes(x = umap_1, y = umap_2,
                   colour = .data[[color_col]]),
               size = fg_size, alpha = fg_alpha, stroke = 0) +
    scale_colour_manual(values = palette, name = legend_title,
                        na.value = "grey60") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_umap_avm(base_size = base_size)

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)))
  }
  p
}

# ── umap_by_category ─────────────────────────────────────────────────────────
# Single-layer UMAP colouring every cell by a categorical variable. Used for
# Fig 4c (atlas coloured by EC subtype) and for the cell-type reference panel
# in per-patient eFigs. Input df needs columns `umap_1`, `umap_2`, `color_col`.
umap_by_category <- function(df, color_col,
                             palette = EC_SUBTYPE_COLORS,
                             point_size = 0.35,
                             alpha = 0.75,
                             legend_title = NULL,
                             legend_ncol = 2,
                             base_size = 14) {
  ggplot(df, aes(x = umap_1, y = umap_2,
                 colour = .data[[color_col]])) +
    geom_point(size = point_size, alpha = alpha, stroke = 0) +
    scale_colour_manual(values = palette, name = legend_title,
                        na.value = "grey70") +
    guides(colour = guide_legend(
      ncol = legend_ncol,
      override.aes = list(size = 2.8, alpha = 1)
    )) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_umap_avm(base_size = base_size)
}

# ── vascular_axis_barplot ────────────────────────────────────────────────────
# Ordered horizontal bar of a per-subtype statistic along the arteriovenous
# axis (artery -> vein), with a FDR-significance star overlay. Used by Fig 4f
# (per-subtype mutant-cell fraction). Input df must have columns:
#   `cell_type` — categorical, must match EC_SUBTYPE_COLORS keys
#   value_col   — numeric column to plot (e.g. pct_mutant, fold_enrichment)
#   sig_col     — optional column with pre-computed significance label
# Order defaults to artery -> vein; off-axis types are appended at the bottom.
vascular_axis_barplot <- function(df,
                                  value_col = "pct_mutant",
                                  sig_col = NULL,
                                  axis_order = c("Large artery", "Artery",
                                                 "Arteriole",
                                                 "Capillary",
                                                 "Angiogenic capillary",
                                                 "Venule", "Vein",
                                                 "Large vein",
                                                 "EndoMT", "Stem-to-EC",
                                                 "Proliferating cell",
                                                 "Mitochondrial",
                                                 "Lymphatic"),
                                  palette = EC_SUBTYPE_COLORS,
                                  x_lab = "Mutant cells (% of cell type)",
                                  base_size = 14) {
  plot_df <- df %>%
    filter(.data$cell_type %in% axis_order) %>%
    mutate(cell_type = factor(.data$cell_type,
                              levels = rev(axis_order[axis_order %in% .data$cell_type])))

  p <- ggplot(plot_df,
              aes(x = .data[[value_col]], y = .data$cell_type,
                  fill = .data$cell_type)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = x_lab, y = NULL) +
    theme_avm(base_size = base_size) +
    theme(panel.grid.major.y = element_blank())

  if (!is.null(sig_col) && sig_col %in% names(plot_df)) {
    p <- p + geom_text(aes(label = .data[[sig_col]]),
                       hjust = -0.25, size = base_size / 3.2,
                       colour = "grey25")
  }
  p
}

# ── vascular_axis_lollipop ───────────────────────────────────────────────────
# Lollipop variant of vascular_axis_barplot. Segment length encodes the
# statistic (e.g. pct_mutant), dot size encodes sample support (e.g.
# n_mutant). Used where the raw bar lengths compress low-count subtypes
# that are nonetheless significant — the dot size surfaces the underlying
# cell-count driving each bar.
#
# Input df must have columns:
#   cell_type   — categorical (must match axis_order entries)
#   value_col   — numeric length of the lollipop stick
#   size_col    — numeric count driving point size
#   sig_col     — optional pre-built text to paste at the tip of each dot
vascular_axis_lollipop <- function(df,
                                   value_col = "pct_mutant",
                                   size_col  = "n_mutant",
                                   sig_col   = NULL,
                                   direction_col = NULL,
                                   axis_order = c("Large artery", "Artery",
                                                  "Arteriole",
                                                  "Capillary",
                                                  "Angiogenic capillary",
                                                  "Venule", "Vein",
                                                  "Large vein",
                                                  "EndoMT", "Stem-to-EC",
                                                  "Proliferating cell",
                                                  "Mitochondrial",
                                                  "Lymphatic"),
                                   palette = EC_SUBTYPE_COLORS,
                                   x_lab = "Mutant cells (% of cell type)",
                                   size_title = "Mutant cells (n)",
                                   size_range = c(2, 9),
                                   base_size = 14) {
  plot_df <- df %>%
    filter(.data$cell_type %in% axis_order) %>%
    mutate(cell_type = factor(.data$cell_type,
                              levels = rev(axis_order[axis_order %in% .data$cell_type])))

  # Direction-of-effect encoding (audit 05): when `direction_col` is supplied,
  # depleted rows (BH-FDR < 0.05 AND fold < 1) render as hollow rings, while
  # every other row renders as a filled disc. This surfaces the significantly
  # depleted Capillary subtype (fold = 0.44, BH-q = 0.003) which would
  # otherwise be visually indistinguishable from the ns rows. Enrichment is
  # already conveyed by stick length × dot size (long stick + large dot = a
  # subtype where many cells are mutant at high fraction), so the shape cue
  # is reserved for depletion. No inline FDR text — encoding is described in
  # the figure caption.
  if (!is.null(direction_col) && direction_col %in% names(plot_df)) {
    plot_df$.is_depleted <- factor(
      ifelse(plot_df[[direction_col]] == "depleted",
             "Depleted (BH-FDR < 0.05)",
             "Other"),
      levels = c("Other", "Depleted (BH-FDR < 0.05)")
    )
  }

  p <- ggplot(plot_df,
              aes(x = .data[[value_col]], y = .data$cell_type,
                  colour = .data$cell_type)) +
    geom_segment(aes(x = 0, xend = .data[[value_col]],
                     y = .data$cell_type, yend = .data$cell_type),
                 linewidth = 0.7, alpha = 0.9)

  if (!is.null(direction_col) && direction_col %in% names(plot_df)) {
    p <- p + geom_point(aes(size = .data[[size_col]],
                            shape = .data$.is_depleted),
                        stroke = 1.1) +
      scale_shape_manual(
        values = c("Other" = 16,
                   "Depleted (BH-FDR < 0.05)" = 1),
        name   = NULL,
        drop   = FALSE,
        guide  = guide_legend(override.aes = list(size = 4, colour = "grey25"))
      )
  } else {
    p <- p + geom_point(aes(size = .data[[size_col]]))
  }

  p <- p + scale_colour_manual(values = palette, guide = "none") +
    scale_size_continuous(range = size_range, name = size_title,
                          breaks = scales::breaks_pretty(4)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = x_lab, y = NULL) +
    theme_avm(base_size = base_size) +
    theme(panel.grid.major.y = element_blank())

  if (!is.null(sig_col) && sig_col %in% names(plot_df)) {
    p <- p + geom_text(aes(label = .data[[sig_col]]),
                       hjust = -0.25, size = base_size / 3.2,
                       colour = "grey25", show.legend = FALSE,
                       nudge_x = 0.05)
  }
  p
}

# ── null_vs_avm_forest ───────────────────────────────────────────────────────
# Thin wrapper around table_forest_general() that renders the TL (temporal-
# lobe) null vs AVM enrichment OR per EC subtype. Input must be a tibble
# with one row per (cell_type, mutation) pair and columns:
#   cell_type  — EC subtype label
#   mutation   — variant label (e.g. "KRAS G12D")
#   odds_ratio, ci_lo, ci_hi — from tl_null_per_subtype.R Fisher test
#   p_adj      — FDR-adjusted p (used to build the right-side p-value column)
#   avm_n_cells — sample-size badge on forest dot
null_vs_avm_forest <- function(df,
                               mutation = NULL,
                               point_col = "#CC0066",
                               x_lab = "AVM vs temporal-lobe odds ratio (log scale, 95% CI)",
                               axis_ticks = c(0.25, 1, 4, 16, 64, 256),
                               base_size = 14) {
  d <- df
  if (!is.null(mutation)) {
    d <- d %>% filter(.data$mutation == !!mutation)
  }
  forest_df <- d %>%
    mutate(
      label = sprintf("%s  [%s]", .data$cell_type, .data$mutation),
      est   = pmax(.data$odds_ratio, 1e-3),
      lo    = pmax(.data$ci_lo, 1e-3),
      hi    = pmax(.data$ci_hi, 1e-3),
      n_label = as.character(.data$avm_n_cells),
      p_label = ifelse(is.na(.data$p_adj), "",
                       sprintf("FDR %s", format.pval(.data$p_adj,
                                                     digits = 2, eps = 1e-4)))
    ) %>%
    arrange(desc(.data$est)) %>%
    select(.data$label, .data$est, .data$lo, .data$hi,
           .data$n_label, .data$p_label)

  table_forest_general(
    forest_df,
    cols = list(list(col = "p_label", header = "FDR p-adj")),
    null_value = 1,
    log_scale = TRUE,
    axis_ticks = axis_ticks,
    point_col = point_col,
    x_lab = x_lab,
    est_col_header = "OR (95% CI)",
    left_label  = "\u2190 Under-represented vs TL",
    right_label = "Enriched vs TL \u2192",
    base_size = base_size
  )
}

# ── SupplementaryTable workbook helpers ──────────────────────────────────────
# write_supp_table_sheets() lets multiple scripts contribute sheets to the same
# SuppTable xlsx without clobbering each other. If the file does not yet exist
# it is created; if it does, any sheet with a matching name is removed first
# (so reruns are idempotent) and the new sheets are appended. This is the
# "live render" pattern used to assemble SuppTable05 (script 12/14 + 15) and
# SuppTable10 (script 13 + 14).
#
# Args:
#   path   : absolute or here()-style path to the target .xlsx.
#   sheets : named list where names are sheet names and values are data frames.
write_supp_table_sheets <- function(path, sheets) {
  stopifnot(is.list(sheets), !is.null(names(sheets)),
            all(nzchar(names(sheets))))

  # 2026-04-26: openxlsx 4.2.8.1 has a sheetId round-trip bug — every
  # loadWorkbook -> saveWorkbook cycle corrupts the workbook so the next
  # loadWorkbook fails with "subscript out of bounds" inside an internal
  # gsub over wb$workbook$sheets. Producers 13 and 14 both contribute
  # sheets to SuppTable10 via this helper, and the second producer trips
  # the bug whenever the first has already written.
  #
  # Workaround: never call openxlsx::loadWorkbook. Instead read existing
  # sheets via readxl (robust), merge with the producer's new sheets
  # (new wins on name collision so re-runs stay idempotent), and write
  # the whole file fresh via openxlsx::write.xlsx. write.xlsx is the
  # high-level entry point that builds a clean workbook from a named
  # list, side-stepping the load/save round-trip path entirely.
  merged <- list()
  if (file.exists(path)) {
    existing_names <- tryCatch(readxl::excel_sheets(path),
                               error = function(e) character(0))
    for (nm in existing_names) {
      if (nm %in% names(sheets)) next  # producer is rewriting this sheet
      df <- tryCatch(
        as.data.frame(readxl::read_excel(path, sheet = nm)),
        error = function(e) NULL
      )
      if (!is.null(df)) merged[[nm]] <- df
    }
  }
  for (nm in names(sheets)) merged[[nm]] <- sheets[[nm]]

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(merged, path, overwrite = TRUE)
  invisible(path)
}

# ── Resolver-prefix lookup ───────────────────────────────────────────────────
# Given a registry panel token, return the human-navigation filename prefix
# from the current resolver assignment in panel_assignments.rds:
#
#   "Fig2A"   for a panel of a main-figure group at manuscript Fig 2 panel a
#   "Fig5"    for a single-panel main figure (no panel letter)
#   "ED01A"   for an ED panel at manuscript ED 1, panel a
#   "ED03"    for a single-panel ED at manuscript ED 3
#
# Numbering: main figures use single-digit (Fig1..Fig9 only — we have 6),
# ED uses 2-digit zero-padded (ED01..ED24+) so files sort correctly in
# Finder / git diff. Panel letters are uppercased.
#
# Returns NULL when the token is not in the registry or panel_assignments.rds
# is missing — callers should fall back to writing token-only filenames in
# that case.
#
# Used by analysis/pipeline/helpers/sync_panel_prefixes.R and (optionally) by producer
# save helpers if they want to write the prefixed copy directly. Composer
# scripts read <token>.rds so they don't need to look up the prefix.
compute_panel_prefix <- function(token, assignments = NULL) {
  pa <- if (is.null(assignments)) {
    tryCatch({
      source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"),
             local = TRUE)
      load_panel_assignments()
    }, error = function(e) NULL)
  } else assignments
  if (is.null(pa)) return(NULL)

  # Find the group containing this panel token.
  group <- NULL
  for (grp in names(pa$panel_letter)) {
    if (token %in% names(pa$panel_letter[[grp]])) {
      group <- grp; break
    }
  }
  if (is.null(group)) return(NULL)

  letter <- pa$panel_letter[[group]][[token]]

  # Find the track (Fig vs EDFig) the group belongs to.
  for (track in names(pa$group_number)) {
    if (group %in% names(pa$group_number[[track]])) {
      n   <- pa$group_number[[track]][[group]]
      pre <- if (track == "EDFig") "ED" else "Fig"
      n_s <- if (track == "EDFig") sprintf("%02d", n) else sprintf("%d", n)
      # Single-panel groups: drop the letter so the prefix reads cleanly
      # (e.g. ED02_missingness_heatmap, not ED02A_missingness_heatmap).
      n_panels_in_group <- length(pa$panel_letter[[group]])
      let <- if (n_panels_in_group > 1L &&
                 length(letter) > 0L && nzchar(letter)) {
               toupper(letter)
             } else ""
      return(paste0(pre, n_s, let))
    }
  }
  NULL
}

# ── Prose figure tag ─────────────────────────────────────────────────────────
# Given one or more registry panel tokens, return a human-readable figure
# citation built from the current resolver assignment, e.g.
#   "Fig. 1d"                          single panel of a main figure
#   "Fig. 1e,g"                        two+ panels of the same main figure
#   "Extended Data Fig. 6a"            single ED panel
#   "Extended Data Fig. 10c,d"         two+ panels of the same ED figure
#   "Fig. 5"                           single-panel main figure (letter omitted)
#   "Fig. 2e; Extended Data Fig. 10b"  tokens spanning different figures
#
# Mirrors compute_panel_prefix() token->group->track resolution but emits the
# prose form used in manuscript citations (panel_registry.R syntax notes).
# Stats-file headers call this so the figure numbers track the resolver
# instead of going stale when prose re-letters panels. Returns NA_character_
# if any token is absent from the registry / cache (caller can fall back).
panel_prose_tag <- function(tokens, assignments = NULL) {
  pa <- if (is.null(assignments)) {
    tryCatch({
      source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"),
             local = TRUE)
      load_panel_assignments()
    }, error = function(e) NULL)
  } else assignments
  if (is.null(pa)) return(NA_character_)

  resolve_one <- function(token) {
    group <- NULL
    for (grp in names(pa$panel_letter)) {
      if (token %in% names(pa$panel_letter[[grp]])) { group <- grp; break }
    }
    if (is.null(group)) return(NULL)
    letter   <- pa$panel_letter[[group]][[token]]
    n_panels <- length(pa$panel_letter[[group]])
    for (track in names(pa$group_number)) {
      if (group %in% names(pa$group_number[[track]])) {
        n   <- pa$group_number[[track]][[group]]
        pre <- if (track == "EDFig") "Extended Data Fig." else "Fig."
        let <- if (n_panels > 1L && length(letter) > 0L && nzchar(letter)) {
                 letter
               } else ""
        return(list(pre = pre, n = n, let = let))
      }
    }
    NULL
  }

  parts <- lapply(tokens, resolve_one)
  if (any(vapply(parts, is.null, logical(1)))) return(NA_character_)

  # Collapse tokens sharing the same prefix+number into one citation with
  # comma-joined letters; tokens from different figures join with "; ".
  keys <- vapply(parts, function(p) paste(p$pre, p$n), character(1))
  out <- character(0)
  for (k in unique(keys)) {
    gp   <- parts[keys == k]
    pre  <- gp[[1]]$pre
    n    <- gp[[1]]$n
    lets <- vapply(gp, function(p) p$let, character(1))
    lets <- lets[nzchar(lets)]
    out  <- c(out, if (length(lets)) {
      sprintf("%s %d%s", pre, n, paste(lets, collapse = ","))
    } else {
      sprintf("%s %d", pre, n)
    })
  }
  paste(out, collapse = "; ")
}

# ── Producer-side disk-slot resolution ───────────────────────────────────────
# Given a registry panel token, return the absolute on-disk directory the
# producer should write that panel into, derived from the current
# panel_assignments cache:
#
#   Main fig token at Fig N panel L  -> results/Figure<N>/panel_<L>/
#   Main fig token in single-panel grp -> results/Figure<N>/         (flat)
#   ED   token at ED N panel L        -> results/ExtendedData/<grp>/panel_<L>/
#   ED   token in single-panel grp    -> results/ExtendedData/<grp>/  (flat)
#
# Panel-letter and group lookups are identical to compute_panel_prefix();
# this helper just maps them to a directory rather than a filename prefix.
# Returns NULL if the token isn't in the registry or the cache is missing.
#
# Use this from any producer whose tokens may move between groups —
# producers bind the disk path through the cache, so registry edits
# propagate automatically without producer-side hardcoding.
panel_slot_dir <- function(token, assignments = NULL) {
  pa <- if (is.null(assignments)) {
    tryCatch({
      source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"),
             local = TRUE)
      load_panel_assignments()
    }, error = function(e) NULL)
  } else assignments
  if (is.null(pa)) return(NULL)

  group <- NULL
  for (grp in names(pa$panel_letter)) {
    if (token %in% names(pa$panel_letter[[grp]])) {
      group <- grp; break
    }
  }
  if (is.null(group)) return(NULL)

  letter <- pa$panel_letter[[group]][[token]]
  n_panels <- length(pa$panel_letter[[group]])

  for (track in names(pa$group_number)) {
    if (group %in% names(pa$group_number[[track]])) {
      n <- pa$group_number[[track]][[group]]
      base <- if (track == "EDFig") {
        here::here("results", "ExtendedData", group)
      } else {
        here::here("results", sprintf("Figure%d", n))
      }
      if (n_panels > 1L && nzchar(letter)) {
        return(file.path(base, sprintf("panel_%s", toupper(letter))))
      }
      return(base)
    }
  }
  NULL
}

# ── Angioarchitecture-adjusted rupture Cox (Supplementary Table 4) ───────────
# Tidy multivariable Cox for the rupture endpoint, adjusting for the three
# non-collinear angioarchitectural covariates (Spetzler–Martin nidus size,
# deep venous drainage, high-risk feature count). Genotype is releveled to
# panel-negative as the reference so the genotype rows read directly as
# "KRAS G12D/G12V vs. panel-negative" (HR > 1 = earlier rupture), matching the
# §1 prose. Cohort + model definition are kept identical to cox_angio in
# analysis/01_main_analysis/09_F1_km_age.R (which feeds the in-text
# rupt_cox_angio_* manifest keys); the G12D HR here equals
# stats$fig3$rupt_cox_angio_hr_g12d (1.68) by construction. The VIF column is
# the car::vif GVIF^(1/(2*Df)) squared (== ordinary VIF for 1-df terms); the
# ~1.0 values substantiate the "non-collinear" wording. Returns a tibble with
# model_n / model_events attributes.
rupture_angio_cox_table <- function(df) {
  stopifnot(all(c("variant_group", "rupture_category", "age", "sm_size_num",
                  "sm_drainage_num", "n_high_risk_num") %in% names(df)))
  d <- df
  d$km2_event <- ifelse(d$rupture_category == "Ruptured at presentation", 1L,
                 ifelse(d$rupture_category == "Never ruptured",           0L,
                        NA_integer_))
  d$km2_group <- factor(as.character(d$variant_group),
                        levels = c("KRAS G12D", "KRAS G12V", "Negative"))
  keep <- !is.na(d$km2_event) & !is.na(d$age) & !is.na(d$km2_group) &
          !is.na(d$sm_size_num) & !is.na(d$sm_drainage_num) &
          !is.na(d$n_high_risk_num)
  d <- d[keep, , drop = FALSE]
  # Panel-negative reference -> genotype coefficients are the vs-negative HRs.
  d$km2_group <- stats::relevel(d$km2_group, ref = "Negative")

  fit <- survival::coxph(
    survival::Surv(age, km2_event) ~ km2_group + sm_size_num +
      sm_drainage_num + n_high_risk_num, data = d)
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  # coxph has no intercept -> car::vif warns but the GVIFs are valid here.
  v <- suppressWarnings(car::vif(fit))
  vif_df <- if (is.matrix(v)) {
    data.frame(vterm = rownames(v), VIF = (v[, "GVIF^(1/(2*Df))"])^2,
               stringsAsFactors = FALSE)
  } else {
    data.frame(vterm = names(v), VIF = unname(v), stringsAsFactors = FALSE)
  }

  labels <- c(
    "km2_groupKRAS G12D" = "KRAS G12D (vs. panel-negative)",
    "km2_groupKRAS G12V" = "KRAS G12V (vs. panel-negative)",
    "sm_size_num"        = "Spetzler–Martin nidus size",
    "sm_drainage_num"    = "Deep venous drainage",
    "n_high_risk_num"    = "High-risk feature count")
  vmap <- c("km2_groupKRAS G12D" = "km2_group",
            "km2_groupKRAS G12V" = "km2_group",
            "sm_size_num" = "sm_size_num",
            "sm_drainage_num" = "sm_drainage_num",
            "n_high_risk_num" = "n_high_risk_num")

  out <- data.frame(
    Term      = unname(labels[tt$term]),
    HR        = round(tt$estimate, 2),
    `CI low`  = round(tt$conf.low, 2),
    `CI high` = round(tt$conf.high, 2),
    P         = signif(tt$p.value, 2),
    VIF       = round(vif_df$VIF[match(unname(vmap[tt$term]), vif_df$vterm)], 2),
    check.names = FALSE, stringsAsFactors = FALSE)
  attr(out, "model_n")      <- fit$n
  attr(out, "model_events") <- fit$nevent
  out
}

# ── Paths (relative to project root) ─────────────────────────────────────────
paths <- list(
  raw_data      = here::here("data", "raw"),
  processed_data = here::here("data", "processed"),
  helpers       = here::here("analysis", "helper_scripts"),
  processing    = here::here("analysis", "00_data_prep"),
  manuscript    = here::here("analysis", "01_main_analysis")
)
