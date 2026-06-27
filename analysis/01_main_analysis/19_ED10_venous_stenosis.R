# 19_ED10_venous_stenosis.R — Theme C (Hale v2)
# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 (was Figure 5 in producer nomenclature): venous-stenosis teaser.
# Output layout follows the same convention as Figs 1/2/3:
#   • results/Figure4/panel_{A..F}/  → individual panel PDF + PNG
#   • results/Figure4/Fig4_composite.{png,pdf} → composite
# The p_5A/p_5B/p_5C ggplot variable names are retained because they're
# purely internal and renaming them adds noise without value.
#
# Andy's framing: *"Put together 1 single figure of these data ... focused on
# framing the venous stenosis genotype-phenotype association ... a 'teaser'
# for a larger omics paper."*  The venous-stenosis signal is the single
# positive phenotype finding that the MDE audit (Extended Data Fig. 20) flags as
# underpowered on the increase side — the observed OR is large (~8) but
# the CI is enormous because of a zero-cell in the Panel-negative arm.
# This figure makes that caveat visually explicit.
#
# Panel design (Apr 2026 revision):
#   4A  Forest plot — venous-stenosis OR from simple 2x2 and from Firth
#       penalised logistic, side-by-side. Annotation calls out "0 events in
#       Panel-negative arm" explicitly so the reader sees the signal
#       comes from the asymmetry, not a mutation-driven rate increase.
#   4B  2x2 waffle / mosaic of genotype x stenosis counts. Empty quadrant
#       for (Negative, stenosis) is the point of the panel.
#   4C  Wälchli et al. brain-EC UMAP coloured by EC subtype — reference
#       map that grounds the coordinates used in panels D and F.
#   4D  Same UMAP with wildtype cells rendered as a faint backdrop and
#       mutant cells highlighted by detected variant (KRAS G12V / Q61H /
#       BRAF V600E). Collapses two previous subplots into one visual.
#   4E  AVM-vs-temporal-lobe null-distribution forest: per-subtype OR for
#       per-read mutant detection in AVM relative to the TL null, showing
#       the venous-compartment signal is real and not a background-rate
#       artefact.
#   4F  Per-subtype mutant-cell fraction along the arteriovenous axis
#       (phase3_vascular_axis_distribution.tsv): Layer 2 q-values
#       (subtype-vs-rest Fisher within AVM, non-self-referential) flag
#       the venule/vein enrichment and capillary depletion. See
#       brain-avm/TODO_2026-05-13.md §A.2-prime for the three-layer
#       design (L1: specificity vs global TL; L2: subtype-vs-rest
#       enrichment; L3: per-subtype TL with Jeffreys smoothing).
#
# The per-study stenosis-ascertainment bar (former panel 4c/4e) remains as
# a standalone Extended Data figure saved below — not wired into the main
# composite.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(here)
  library(patchwork); library(logistf); library(scales); library(stringr)
  library(lme4); library(purrr); library(tibble)
})
source(here("analysis", "helper_scripts", "utils.R"))

# Wrap long panel titles to a sensible character width so they don't get
# truncated when patchwork lays out 3 panels per row at 14 in width
# (each panel ≈ 4.5 in). 30 chars fits two lines comfortably for the
# longest titles without hyphenation.
# Generic title wrapper retained for any panel that genuinely needs
# auto-wrapping. Fig 4 titles below mostly read fine on one line and
# are passed through as raw strings; only Panel E carries an explicit
# "\n" break (placed before "reference" so the wrap doesn't bisect
# "mutant vs").
.wrap_title <- function(s) stringr::str_wrap(s, width = 30)

# Phase 2: composite panel tags are derived from the registry-driven
# resolver cache (written by run_all.R step [1.5]) rather than hardcoded.
# Tags attach via plot_annotation(tag_levels = ...) on the composite so
# each of the six panels (a-f) carries the letter the registry assigns.
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.fig4_tags <- tryCatch({
  a <- load_panel_assignments()
  list(stenosis_waffle = panel_letter("stenosis_waffle", a))
}, error = function(e) {
  message("  [panel_assignments cache missing] falling back to default")
  list(stenosis_waffle = "a")
})

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
# 2026-05-20: panel_dims provides save_panel_native(token, plot) for
# tokens registered in FIGURE_LAYOUTS. Used for stenosis_waffle below;
# stenosis_ascertainment still uses the legacy local save_panel().
source(here("analysis", "pipeline", "helpers", "panel_dims.R"))
# All call sites pass through the cairo device for unicode glyph fidelity.
save_panel <- function(dir, name, plot, w, h)
  save_panel_impl(dir, name, plot, w, h, device = "cairo")

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))

# ── Analysis cohort ─────────────────────────────────────────────────────────
vs <- df %>%
  filter(!is.na(venous_outflow_stenosis_num), !is.na(geno_binary)) %>%
  mutate(geno = factor(geno_binary,
                       levels = c("Panel-negative", "Variant-positive")),
         stenosis = factor(venous_outflow_stenosis_num,
                           levels = c(0, 1),
                           labels = c("No stenosis", "Stenosis")))

tab_2x2 <- with(vs, table(geno, stenosis))
cat("── Figure 4: venous outflow stenosis ──\n")
cat(sprintf("Total n = %d (Mut+ %d, Neg %d)\n",
            nrow(vs), sum(vs$geno == "Variant-positive"),
            sum(vs$geno == "Panel-negative")))
print(tab_2x2)

# ─────────────────────────────────────────────────────────────────────────────
# 4A — Forest of stenosis OR: simple 2x2 (0.5 continuity) vs Firth
# ─────────────────────────────────────────────────────────────────────────────
a <- tab_2x2["Variant-positive", "Stenosis"]
b <- tab_2x2["Variant-positive", "No stenosis"]
c <- tab_2x2["Panel-negative", "Stenosis"]
d <- tab_2x2["Panel-negative", "No stenosis"]

# Simple 2x2 with 0.5 continuity correction (the "naive" OR used in MDE)
ac <- a + 0.5; bc <- b + 0.5; cc <- c + 0.5; dc <- d + 0.5
or_simple  <- (ac * dc) / (bc * cc)
se_simple  <- sqrt(1/ac + 1/bc + 1/cc + 1/dc)
lo_simple  <- exp(log(or_simple) - 1.96 * se_simple)
hi_simple  <- exp(log(or_simple) + 1.96 * se_simple)

# Firth-penalised logistic, sample-type-adjusted (matches the model used in
# Supp Table 3 high-risk-feature panel so the headline OR in Fig 4 and the
# supplementary table are identical). Adjustment for sample_type_clean is
# motivated by the ascertainment caveat: stenosis calls originated exclusively
# from one institution whose sample_type composition differs from the other
# contributing series.
#
# AUDIT 2026-05-12 — ACKNOWLEDGED structural limitation: all positive
# venous-outflow-stenosis calls originate from BCH (stenosis_n_bch == total
# positives, see stats dump below). study_clean cannot be added to this Firth
# model without inducing perfect separation, so the institutional axis is not
# adjusted for. The OR/CI/P emitted here is, by construction, a "BCH-only
# subseries" estimate. §5 prose must label this as "consistent but
# ascertainment-confounded" rather than as supportive non-null evidence.
vs_adj <- vs %>% filter(!is.na(sample_type_clean))
firth_fit <- logistf(venous_outflow_stenosis_num ~ geno + sample_type_clean,
                     data = vs_adj)
firth_i   <- match("genoVariant-positive", names(coef(firth_fit)))
or_firth  <- exp(as.numeric(coef(firth_fit)[firth_i]))
lo_firth  <- exp(as.numeric(firth_fit$ci.lower[firth_i]))
hi_firth  <- exp(as.numeric(firth_fit$ci.upper[firth_i]))
p_firth   <- as.numeric(firth_fit$prob[firth_i])

# 2026-04-30: dropped forest_df + p_5A hand-rolled forest. The panel
# was retired from Fig 4 on 2026-04-27 (duplicate of one row of
# Fig 3 hr_features_OR), so the ~120 lines that built it were dead
# code. Firth + simple-2x2 OR statistics still flow into the prose
# stats and the schema dump below \u2014 those don't need a ggplot.
# If the panel is ever revived, build via table_forest() (utils.R)
# rather than re-rolling annotate() geometry.


# ─────────────────────────────────────────────────────────────────────────────
# 4B — 2x2 waffle: each cell is a grid of dots, one dot = one patient.
# Purpose: make the empty (Negative, Stenosis) cell impossible to miss.
# ─────────────────────────────────────────────────────────────────────────────
counts <- vs %>% count(geno, stenosis, name = "n_cases")

build_waffle <- function(n, ncol = 10) {
  if (n == 0) return(tibble(x = numeric(0), y = numeric(0)))
  idx <- seq_len(n) - 1
  tibble(x = idx %% ncol, y = idx %/% ncol)
}

waffle_df <- counts %>%
  mutate(dots = lapply(n_cases, build_waffle)) %>%
  tidyr::unnest(dots) %>%
  mutate(
    geno     = factor(ifelse(geno == "Panel-negative",
                             "Panel-neg", "Mut+"),
                      levels = c("Panel-neg", "Mut+")),
    stenosis = factor(stenosis,
                      levels = c("No stenosis", "Stenosis"))
  )

cell_labels <- counts %>%
  mutate(
    # Must use the SAME factor levels as waffle_df — otherwise
    # facet_grid takes the union of levels across both data sources
    # and renders 4 columns instead of 2.
    geno  = factor(ifelse(geno == "Panel-negative", "Panel-neg", "Mut+"),
                   levels = c("Panel-neg", "Mut+")),
    stenosis = factor(stenosis, levels = c("No stenosis", "Stenosis")),
    label = sprintf("n = %d", n_cases)
  )

# 10-wide waffle: 10 columns, rows grow with count. Max cell is 136/10 = 14
# rows. Letting the panel stretch (no coord_fixed) fills the cell — dots
# become vertically elongated but the 1-dot-per-patient reading still holds
# and the facet strips no longer clip at composite width.
p_5B <- ggplot(waffle_df, aes(x = x, y = y,
                              color = stenosis, shape = stenosis)) +
  geom_point(size = 1.9) +
  facet_grid(stenosis ~ geno, switch = "y") +
  geom_text(data = cell_labels,
            aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.6,
            size = TYPO$geom_text_label, color = "grey25") +
  # Tier A reserved for variant identity — "Stenosis" is an abnormal-finding
  # category, mapped to the Tier E event/alarm token (red).
  scale_color_manual(values = c("No stenosis" = "grey75",
                                "Stenosis"    = unname(PAL_EVENT[["Event"]]))) +
  scale_shape_manual(values = c("No stenosis" = 16, "Stenosis" = 17)) +
  scale_y_reverse(expand = expansion(mult = c(0.08, 0.04))) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  # 2026-05-19 NM pass: title dropped (moves to caption).
  labs(title = NULL,
       x = NULL, y = NULL, color = NULL, shape = NULL) +
  guides(color = "none", shape = "none") +
  # 2026-05-20: switched theme_avm() -> theme_avm_native() so the
  # waffle saves at NATIVE base_size 9 with the producer's
  # element_blank() overrides applied AFTER (and therefore preserved).
  # save_panel_native() does not re-theme so axis.text stays hidden.
  theme_avm_native() +
  theme(axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        panel.grid       = element_blank(),
        strip.text.x      = element_text(face = "bold"),
        strip.text.y.left = element_text(face = "bold", angle = 0),
        aspect.ratio     = 1)

# ─────────────────────────────────────────────────────────────────────────────
# 4e-sister (kept in Extended Data) — Per-study stenosis-ascertainment bar.
# ─────────────────────────────────────────────────────────────────────────────
study_df <- df %>%
  filter(!is.na(venous_outflow_stenosis_num)) %>%
  count(study, stenosis_present = venous_outflow_stenosis_num == 1,
        name = "n_cases") %>%
  tidyr::pivot_wider(names_from = stenosis_present, values_from = n_cases,
                     values_fill = 0L) %>%
  rename(stenosis = `TRUE`, no_stenosis = `FALSE`) %>%
  mutate(total = stenosis + no_stenosis) %>%
  arrange(desc(total)) %>%
  mutate(study = factor(study, levels = rev(study)))

study_plot_df <- study_df %>%
  tidyr::pivot_longer(c(stenosis, no_stenosis),
                      names_to = "status", values_to = "n_cases") %>%
  mutate(status = factor(status,
                         levels = c("no_stenosis", "stenosis"),
                         labels = c("No stenosis", "Stenosis")))

p_5C_study <- ggplot(study_plot_df,
               aes(x = n_cases, y = study, fill = status)) +
  geom_col(width = 0.7) +
  # Tier A reserved for variant identity — "Stenosis" maps to Tier E alarm (red).
  scale_fill_manual(values = c("No stenosis" = "grey75",
                               "Stenosis"    = unname(PAL_EVENT[["Event"]]))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  # v6.48 (2026-05-21): wrap long study names on the y-axis so they fit
  # inside the 3.30-in cell width. The previous one-line labels were
  # overflowing into the bar area.
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  # v6.48: x-axis title wrapped to two lines so it doesn't overflow.
  labs(x = stringr::str_wrap(
         "Number of patients with imaging assessed for venous outflow stenosis",
         width = 32),
       y = NULL, fill = NULL) +
  theme_nature_panel() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ─────────────────────────────────────────────────────────────────────────────
# Composite & save
#
# 2-D grid (matches Figs 1/2/3 aspect-ratio conventions; prior single-column
# 9x20 stack was overflowing the Nat Med page area):
#
#   ┌─────────────────────────────────────────────┐
#   │              A  forest (full width)         │   row 1
#   ├──────────────┬───────────────┬──────────────┤
#   │   B waffle   │  C UMAP sub   │  D UMAP var  │   row 2
#   ├──────────────┴───────────────┴──────────────┤
#   │            E  comp stacked bar              │   row 3
#   ├─────────────────────────────────────────────┤
#   │            F  subtype lollipop              │   row 4
#   └─────────────────────────────────────────────┘
#
# Design-string approach keeps panels as live ggplots (no PNG stitching),
# so tag letters and guides still resolve correctly.
# ─────────────────────────────────────────────────────────────────────────────
# 2026-05-19 scRNA retirement: original Fig 4 (6 panels, 4 scRNA-derived)
# dissolved. Survivors are stenosis_waffle (Phase 3 merges it into the
# gxp_associations group as the last panel of the new Fig 2) and
# stenosis_ascertainment (ED panel inside ed_null_audit). Output paths
# are registry-resolved via panel_slot_dir() so the resolver places
# them in the correct Figure<N>/panel_<letter>/ subdir without producer
# changes when the group structure shifts.
# v6.46 (2026-05-21): stenosis_waffle was retired from Fig 2 in v3; the
# resolver no longer maps the token to a slot dir, so panel_slot_dir
# returns NULL and dir.create blows up the producer before the
# stenosis_ascertainment save below. Gate the legacy lookup so the rest
# of the producer can still run.
waffle_dir <- tryCatch(panel_slot_dir("stenosis_waffle"),
                       error = function(e) NULL)
if (!is.null(waffle_dir))
  dir.create(waffle_dir, recursive = TRUE, showWarnings = FALSE)
# 2026-05-20 v3: stenosis_waffle retired from Fig 2 (the
# venous-outflow-stenosis finding is fully represented as a panel of the
# hr_features_OR forest). The p_5B ggplot is kept above so the
# stenosis-ascertainment ED panel and Supplementary Table workflows
# downstream can still inspect it if needed, but it no longer ships as
# a Fig 2 panel and is not saved.
# (formerly: save_panel_native("stenosis_waffle", p_5B))

efig_asc_dir <- panel_slot_dir("stenosis_ascertainment")
dir.create(efig_asc_dir, recursive = TRUE, showWarnings = FALSE)
save_panel(efig_asc_dir, "stenosis_ascertainment", p_5C_study, 3.30, 3.04)

# Composite assembly is now centralised in 34_ED10_assemble.R
# The legacy inline
# block below is gated behind if(FALSE){...} as historical documentation
# of the layout intent; canonical Fig 4 composite path stays at
# results/Figure4/Fig4_composite.{pdf,png}.
cat("\n── Fig 4 composite: deferred to 34_ED10_assemble.R ──\n")


# ── Stats dump ──────────────────────────────────────────────────────────────
stats_dir <- here("results", "stats")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)



writeLines(c(
  sprintf(paste0("# Venous outflow stenosis — supporting statistics ",
                 "(effect in %s high-risk-feature forest; ",
                 "per-study ascertainment in %s)"),
          panel_prose_tag("hr_features_OR"),
          panel_prose_tag("stenosis_ascertainment")),
  sprintf("Generated %s", Sys.Date()),
  sprintf("Total n = %d  |  Mut+ = %d  |  Neg = %d  |  Adjusted model n = %d",
          nrow(vs), sum(vs$geno == "Variant-positive"),
          sum(vs$geno == "Panel-negative"), nrow(vs_adj)),
  "",
  "## 2x2 table",
  capture.output(print(tab_2x2)),
  "",
  "## Effect sizes",
  sprintf("  Simple 2x2 + 0.5 continuity : OR = %.2f (95%% CI %.2f\u2013%.2f)",
          or_simple, lo_simple, hi_simple),
  sprintf("  Firth, sample-type adjusted : OR = %.2f (95%% CI %.2f\u2013%.2f),  p = %.3g",
          or_firth, lo_firth, hi_firth, p_firth),
  "",
  "## Per-study ascertainment (Supp Table 3 footnote source)",
  capture.output(print(study_df, n = Inf))
),
  file.path(stats_dir, "fig4_venous_stenosis.txt"))

# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# Stats manifest fragments \u2014 fig4 + edfig14_stenosis_asc
# -----------------------------------------------------------------------------
# Every number the Results \u00a75 (stenosis) and \u00a76 (scRNA localisation) prose
# cites is emitted here into the stats manifest.
# Panel-F row ORs + BH-FDRs are extracted as scalars for inline prose and
# also passed through as a full data frame for downstream ED consumers.
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

# Per-study ascertainment cells (for \u00a75 "all 10 positive calls from BCH" line
# and for ED14 ascertainment bar). study_df is built earlier in this script.
.study_get <- function(study_df, study_name, col) {
  row <- study_df[as.character(study_df$study) == study_name, , drop = FALSE]
  if (nrow(row) != 1L) return(0L)
  as.integer(row[[col]])
}
fig4_fragment <- list(
  # Panel a + b \u2014 stenosis 2x2 + forest
  stenosis_n_total           = nrow(vs),
  stenosis_n_mutpos          = unname(sum(vs$geno == "Variant-positive")),
  stenosis_n_neg             = unname(sum(vs$geno == "Panel-negative")),
  stenosis_n_mutpos_stenosis = unname(a),
  stenosis_n_neg_stenosis    = unname(c),
  or_simple                  = unname(or_simple),
  or_simple_lo               = unname(lo_simple),
  or_simple_hi               = unname(hi_simple),
  or_firth                   = unname(or_firth),
  or_firth_lo                = unname(lo_firth),
  or_firth_hi                = unname(hi_firth),
  or_firth_p                 = unname(p_firth),
  firth_n_adjusted           = nrow(vs_adj),
  # Per-study ascertainment (\u00a75 caveat cites "10/10 calls from BCH")
  stenosis_n_bch             = .study_get(study_df, "BCH",  "total"),
  stenosis_n_chop            = .study_get(study_df, "CHOP", "total"),
  stenosis_n_uab             = .study_get(study_df, "UAB",  "total"),
  stenosis_total_positive    = sum(study_df$stenosis),
  stenosis_total_positive    = sum(study_df$stenosis)
)

write_stats_section(section = "fig4", stats = fig4_fragment)

# ---- edfig14 \u2014 stenosis ascertainment bar --------------------------------
edfig14_fragment <- list(
  n_total         = sum(study_df$total),
  n_positive      = sum(study_df$stenosis),
  n_bch           = .study_get(study_df, "BCH",  "total"),
  n_bch_positive  = .study_get(study_df, "BCH",  "stenosis"),
  n_chop          = .study_get(study_df, "CHOP", "total"),
  n_chop_positive = .study_get(study_df, "CHOP", "stenosis"),
  n_uab           = .study_get(study_df, "UAB",  "total"),
  n_uab_positive  = .study_get(study_df, "UAB",  "stenosis"),
  study_table     = as.data.frame(study_df)
)
write_stats_section(section = "edfig14_stenosis_asc", stats = edfig14_fragment)

cat("\n\u2713 venous-stenosis panels complete.\n")
if (!is.null(waffle_dir))
  cat(sprintf("   stenosis_waffle:        %s\n", file.path(waffle_dir, "stenosis_waffle.png")))
cat(sprintf("   stenosis_ascertainment: %s\n", file.path(efig_asc_dir, "stenosis_ascertainment.png")))
