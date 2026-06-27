# 07_F1_consort_flow.R — Figure 1A: CONSORT-style cohort flow diagram
#
# Live-rendered, JAMA-styled successor to the previous Graphviz/DiagrammeR
# implementation. All counts (per-site, per-published-series, pre/post-
# exclusion, dPCR vs literature, mut+/mut-/pending, per-variant breakdown)
# are pulled from the live data + manifest so cohort refreshes propagate
# without manual edits.
#
# Tier structure (refined v2 of the BioRender / Python / TikZ prototypes):
#
#   Tier 1  | 3 institutional cohorts |   | 5 published series |
#                       \                    /
#   Tier 2          [ N pre-exclusion bAVM tissues ]  ──→ (excluded — dashed)
#                                |
#   Tier 3                [ N sporadic bAVMs ]         ──→ (genotype pending)
#                          /                \
#   Tier 4    [ dPCR genotyping ]   [ Patient-level meta-analysis ]
#                          \                /
#   Tier 5         [ Genotype-phenotype harmonization ]
#                          /                \
#   Tier 6     [ Variant positive ]   [ Panel negative ]
#
# Two side-arms are now CONSORT-standard: an exclusion side-arm off Tier 2
# (dashed border, dashed arrow — these cases leave the analysis pipeline)
# and a pending side-arm off Tier 3 (solid border, solid arrow — these
# cases stay in the cohort but bypass the genotyping/harmonization arms).
# The variant breakdown that used to live inside the Variant-positive box
# was removed because it duplicates Fig 1B (variant_landscape).
#
# Visual conventions (JAMA Network–style):
#   * Sans-serif (Arial), 11–13 pt; bold for box headers, regular for body
#   * Uniform white-fill / dark-gray stroke for every flow node — no
#     coral/amber/blue accents on the outcome boxes; categorical
#     differences come from labels, not color
#   * Excluded side-arm: dashed border + dashed arrow (lighter penwidth)
#   * All other arrows uniform: solid, single penwidth, top-to-bottom on
#     the main spine; horizontal-only on side-arms (no diagonal routing).
#
# Output: results/Figure1/Fig1A_consort_flow.{pdf,png}
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(here)
  library(DiagrammeR)
  library(dplyr)
})

source(here::here("analysis", "helper_scripts", "utils.R"))   # panel_slot_dir

# Canonical slot: Figure1/panel_A/ matches panel_slot_dir("consort_flow")
# so sync_panel_prefixes' zeroth-pass placement check leaves the file in
# place. Previous location (results/Figure1/) was mis-classified as
# "registry-misplaced" and swept on every build.
output_dir <- panel_slot_dir("consort_flow")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Pull live counts
# ─────────────────────────────────────────────────────────────────────────────
df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))

n_post_excl <- nrow(df)  # sporadic bAVMs (analysis-ready)

.les <- table(df$study_clean)
.get <- function(s) as.integer(.les[[s]] %||% 0L)

n_bch        <- .get("BCH")
n_chop       <- .get("CHOP")
n_uab        <- .get("UAB")
n_inst       <- n_bch + n_chop + n_uab
n_nikolaev   <- .get("Nikolaev")
n_priemer    <- .get("Priemer")
n_hong       <- .get("Hong")
n_goss       <- .get("Goss")
n_gao        <- .get("Gao")
n_pub        <- n_nikolaev + n_priemer + n_hong + n_goss + n_gao

# CONSORT integrity check: Tier 1 sum must equal Sporadic.
if ((n_inst + n_pub) != n_post_excl) {
  warning(sprintf(
    "CONSORT reconciliation: inst %d + pub %d = %d ≠ Sporadic %d. Investigate.",
    n_inst, n_pub, n_inst + n_pub, n_post_excl
  ))
}

# Genotyping arms — every analysis-ready row has a mutation call after the
# May 2026 master ingest (n_pending = 0). After the CHOP46-orphan study
# imputation in 01_clean_master.R, every row has study_clean ∈ {BCH,CHOP,
# UAB} for dPCR or {Nikolaev,Priemer,Hong,Goss,Gao} for meta; the
# Tier-4 residual category is no longer needed.
df_genotyped <- df %>% filter(!is.na(mutation_positive))
n_pending    <- nrow(df) - nrow(df_genotyped)

n_dpcr_arm   <- df_genotyped %>% filter(study_clean %in% c("BCH", "CHOP", "UAB")) %>% nrow()
n_meta_arm   <- df_genotyped %>%
  filter(study_clean %in% c("Nikolaev", "Priemer", "Hong", "Goss", "Gao")) %>% nrow()
n_arm_unassigned <- nrow(df_genotyped) - n_dpcr_arm - n_meta_arm
if (n_arm_unassigned != 0L) {
  warning(sprintf(
    "Tier-4 arm assignment incomplete: %d analysis-ready row(s) have study_clean ∉ {institutional, published}. Expected 0 after CHOP46-orphan imputation in 01_clean_master.R.",
    n_arm_unassigned
  ))
}

# Outcome counts
n_mut_pos    <- sum(df$mutation_positive == TRUE, na.rm = TRUE)
n_mut_neg    <- sum(df$mutation_positive == FALSE, na.rm = TRUE)

# Per-variant breakdown intentionally omitted from this panel — Fig 1B
# (variant_landscape) carries that information, and including it inside
# the Variant-positive box made the bottom row visually unbalanced
# against Panel-negative.

# ─────────────────────────────────────────────────────────────────────────────
# 2. Reusable render helper (rsvg → PDF + PNG, fallback to HTML widget)
# ─────────────────────────────────────────────────────────────────────────────
render_consort <- function(dot_string, pdf_path, png_path,
                           width_in = 8, height_in = 11, dpi = 300) {
  g <- grViz(dot_string)

  tryCatch({
    if (!requireNamespace("DiagrammeRsvg", quietly = TRUE))
      stop("DiagrammeRsvg not installed")
    if (!requireNamespace("rsvg", quietly = TRUE))
      stop("rsvg not installed")

    svg_text <- DiagrammeRsvg::export_svg(g)
    svg_raw  <- charToRaw(svg_text)

    rsvg::rsvg_pdf(svg_raw, pdf_path,
                   width = width_in * 72, height = height_in * 72)
    rsvg::rsvg_png(svg_raw, png_path,
                   width = width_in * dpi, height = height_in * dpi)

    cat(sprintf("  Saved: %s\n  Saved: %s\n", pdf_path, png_path))
  }, error = function(e) {
    html_path <- sub("\\.(pdf|png)$", ".html", pdf_path)
    htmlwidgets::saveWidget(g, html_path, selfcontained = TRUE)
    cat(sprintf("  Saved HTML (install DiagrammeRsvg + rsvg for PDF/PNG): %s\n",
                html_path))
    cat(sprintf("  Error: %s\n", e$message))
  })

  invisible(g)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. DOT specification — JAMA-styled CONSORT
# ─────────────────────────────────────────────────────────────────────────────
# Style choices, all uniform across nodes:
#   * fontname  = "Arial"     — sans-serif throughout
#   * fontsize  = 22 default; 13 for box header line, 11 for body
#   * fillcolor = "#FFFFFF"  — pure white fill
#   * color     = "#3A3A3A"  — dark gray border
#   * fontcolor = "#1A1A1A"  — near-black text
#   * penwidth  = 1.0         — single consistent stroke weight
#
# Sizing rule: every main-spine node has an explicit (width × height)
# applied via fixedsize=true; all boxes within the same tier share the
# same dimensions so the diagram reads as a regular grid. Side-arm
# nodes (excl, pending) are intentionally sized differently from their
# parent tier so the eye distinguishes the spine from the horns at a
# glance. Per-tier dimensions are tuned to fit the densest box in that
# tier (e.g. Tier 1 published-series listing has 7 rows, so Tier 1
# height = 1.8 in; the institutional box pads visually with vertical
# whitespace).
#
# Bolding: ONLY box-header lines (the first <TR>) and the n-total row
# get <B>; per-site / per-series body rows are regular. All <TD> use
# ALIGN="CENTER" so the label sits centered horizontally within the
# box; Graphviz centers the HTML table vertically within fixed-size
# nodes, so single-line / two-line labels float in the middle rather
# than top-aligning.

consort_dot <- sprintf('
digraph CONSORT {

  graph [
    rankdir  = TB
    splines  = ortho
    # v6.4 (2026-05-21): bump nodesep + ranksep so sibling boxes no
    # longer collide horizontally ("All bAVM tissues assembled" was
    # overlapping the published-series box; "Multiplex dPCR somatic
    # tissue genotyping" was overlapping "Patient-level meta-analysis").
    nodesep  = 0.6
    ranksep  = 0.55
    bgcolor  = "transparent"
    fontname = "Arial"
  ]

  /* See R-level comments above for sizing/padding/colour rationale. */
  node [
    shape     = box
    style     = "filled,rounded"
    fillcolor = "#FFFFFF"
    color     = "#3A3A3A"
    fontname  = "Arial"
    fontcolor = "#1A1A1A"
    fontsize  = 22
    penwidth  = 1.0
    margin    = "0.0,0.0"
  ]

  edge [
    color     = "#3A3A3A"
    arrowsize = 0.7
    penwidth  = 1.0
  ]

  /* Tier 1: Source cohorts. v6.12 (2026-05-20): bumped 2.6 x 2.5 ->
     3.2 x 3.0 to stop the bold two-line header ("3 institutional / cohorts",
     "5 published series") from clipping at the box edge in the rendered
     PNG. Both boxes share dimensions + {rank = same} → vertically aligned. */
  inst [width=3.2, height=4.0, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">3 institutional<BR/>cohorts</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER">CHOP, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER">BCH, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER">UAB, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  pub [width=3.2, height=4.0, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">5 published<BR/>series</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER">Nikolaev, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER">Priemer, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER">Hong, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER">Goss, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER">Gao, <B>n = %d</B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  {rank = same; inst; pub}

  /* Tier 2: sporadic bAVM cohort (inst + pub, no exclusion step). */
  sporadic [width=3.1, height=1.2, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">Sporadic bAVMs</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  /* Tier 3 horn previously held a "Genotype pending" side-arm. Retired
     2026-05-19 — every sporadic bAVM now has a mutation call after the
     18 May 26 master ingest (n_pending = 0). */

  /* Tier 4: two analysis arms. v6.18 (2026-05-20): width 2.73 -> 3.0,
     height 1.47 -> 1.8 — labels still clipped at v6.17 dims; bumped
     more aggressively. dpcr and meta scale together. */
  dpcr [width=3.0, height=1.8, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">Multiplex dPCR<BR/>somatic tissue<BR/>genotyping</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  meta [width=3.0, height=1.8, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">Patient-level<BR/>meta-analysis</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  {rank = same; dpcr; meta}

  /* Tier 5: harmonization. v6.18: width 3.12 -> 3.7 (now ~+42%%
     vs original), height 1.0 -> 1.3 — bottom-row labels needed
     more breathing room than the v6.17 +20%% width bump provided. */
  harm [width=3.7, height=1.3, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">Genotype–phenotype</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">harmonization</FONT></B></TD></TR>
</TABLE>>]

  /* Tier 6: outcomes. v6.18: width 2.86 -> 3.3, height 1.0 -> 1.3.
     Bottom-row width/height match the rest of the v6.18 bump. */
  mutpos [width=3.3, height=1.3, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">Variant positive</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  mutneg [width=3.3, height=1.3, fixedsize=true, label = <
<TABLE BORDER="0" CELLPADDING="6" CELLSPACING="0">
  <TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="24">Panel negative</FONT></B></TD></TR>
  <TR><TD ALIGN="CENTER"><B>n = %d</B></TD></TR>
</TABLE>>]

  {rank = same; mutpos; mutneg}

  /* Edges: inst and pub feed directly into sporadic cohort. */
  inst       -> sporadic
  pub        -> sporadic
  sporadic   -> dpcr
  sporadic   -> meta
  dpcr       -> harm
  meta       -> harm
  harm       -> mutpos
  harm       -> mutneg
}',
  # Tier 1 institutional
  n_chop, n_bch, n_uab, n_inst,
  # Tier 1 published
  n_nikolaev, n_priemer, n_hong, n_goss, n_gao, n_pub,
  # Tier 2 sporadic
  n_post_excl,
  # Tier 3 arms
  n_dpcr_arm, n_meta_arm,
  # Tier 5 outcomes
  n_mut_pos, n_mut_neg
)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Render
# ─────────────────────────────────────────────────────────────────────────────
g <- render_consort(
  consort_dot,
  pdf_path = file.path(output_dir, "Fig1A_consort_flow.pdf"),
  png_path = file.path(output_dir, "Fig1A_consort_flow.png"),
  width_in = 9, height_in = 16
)

# NB: 26_F1_assemble.R reads Fig1A_consort_flow.png directly
# via rasterGrob, so no per-panel RDS slot is written here.

cat("══ Fig1A CONSORT flow diagram complete ══\n")
cat(sprintf("  Tier 1 institutional: BCH %d, CHOP %d, UAB %d = %d\n",
            n_bch, n_chop, n_uab, n_inst))
cat(sprintf("  Tier 1 published:     %d\n", n_pub))
cat(sprintf("  Tier 1 sum:           %d\n", n_inst + n_pub))
cat(sprintf("  sporadic bAVMs:       %d\n", n_post_excl))
cat(sprintf("  dPCR arm:             %d\n", n_dpcr_arm))
cat(sprintf("  meta-analysis arm:    %d\n", n_meta_arm))
cat(sprintf("  mut+ / mut- / pending: %d / %d / %d\n",
            n_mut_pos, n_mut_neg, n_pending))
