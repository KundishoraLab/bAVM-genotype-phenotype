# 08_F1_cohort_variants.R — cohort_variants group: Variant landscape +
# VAF distributions (panels: variant_landscape, vaf_by_variant) plus the
# per-series cohort-heterogeneity ED panels (per_series_rate, per_series_vaf).
#
# Input:  data/processed/bAVM_analysis_ready.rds
# Output: variant_landscape, vaf_by_variant (main figure panels) +
#         per_series_rate / per_series_vaf (ed_cohort_heterogeneity panels).
# Manuscript figure / panel letters are resolver-driven — see
# analysis/pipeline/panel_registry.R.
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(patchwork)

source(here("analysis", "helper_scripts", "utils.R"))

# ── 1. Load data ─────────────────────────────────────────────────────────────

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))
genotyped <- df %>% filter(!is.na(mutation_positive))
output_dir <- here("results")
fig1_dir <- file.path(output_dir, "Figure1")
efig_dir <- file.path(output_dir, "ExtendedData")

# save_panel(dir, name, plot, w, h) — canonical helper from
# analysis/pipeline/helpers/save_panel.R. Default device + RDS output;
# deterministic seed per `name` keeps PDF/PNG jitter in sync.
source(here("analysis", "pipeline", "helpers", "save_panel.R"))

cat(sprintf("Genotyped patients: %d\n", nrow(genotyped)))

# ── 2. Use unified palettes from utils.R ─────────────────────────────────────

# Save palettes for downstream scripts (all defined in utils.R)
saveRDS(
  list(geno = PAL_VARIANT, binary = PAL_BINARY, detailed = PAL_DETAILED,
       sm = PAL_SM, study = PAL_STUDY, sample = PAL_SAMPLE),
  file.path(output_dir, "palettes.rds")
)

# ── 3. Fig 1B: Somatic variant distribution (horizontal stacked bar) ────────

landscape_data <- df %>%
  mutate(
    variant_label = case_when(
      mut_KRAS_G12D == 1 ~ "KRAS G12D",
      mut_KRAS_G12V == 1 ~ "KRAS G12V",
      mut_KRAS_G12C == 1 ~ "KRAS G12C",
      mut_KRAS_G12A == 1 ~ "KRAS G12A",
      mut_KRAS_Q61H == 1 ~ "KRAS Q61H",
      mut_KRAS_c191_196dup == 1 ~ "KRAS dup",
      mut_BRAF_V600E == 1 ~ "BRAF V600E",
      mut_BRAF_Q636X == 1 ~ "BRAF Q636X",
      # Canonical display labels (see analysis/helper_scripts/utils.R
      # contract). The Fig 1B variant-landscape y-axis renders these
      # factor levels directly; do not regress to the raw "Negative"
      # or shorthand "mut+" strings.
      mutation_positive == FALSE ~ "Panel-negative",
      is.na(mutation_positive) ~ "No tissue",
      # Catch-all: mutation_positive == TRUE but no specific mut_* dummy
      # set. This should never fire after the explicit enumeration
      # above; if it does, surface as "Unassigned variant-positive" so
      # the figure makes the gap visible rather than silently labelling
      # as "No tissue".
      TRUE ~ "Unassigned variant-positive"
    ),
    variant_label = factor(variant_label, levels = rev(c(
      "KRAS G12D", "KRAS G12V", "KRAS G12C", "KRAS G12A", "KRAS Q61H", "KRAS dup",
      "BRAF V600E", "BRAF Q636X", "Panel-negative",
      "Unassigned variant-positive", "No tissue"
    )))
  )

variant_counts <- landscape_data %>%
  count(variant_label) %>%
  mutate(pct = 100 * n / sum(n))

fig1b <- ggplot(variant_counts, aes(x = n, y = variant_label, color = variant_label)) +
  geom_segment(aes(x = 0, xend = n, yend = variant_label),
               color = "grey70", linewidth = 0.6) +
  # v6.11 (2026-05-20): alpha 0.85 -> 1.0 so the lollipop stem is fully
  # occluded behind each disc (was bleeding through at 0.85).
  geom_point(aes(size = n), alpha = 1.0) +
  scale_color_manual(values = PAL_DETAILED, guide = "none") +
  # v6.11 (2026-05-20): max_size 4.5 -> 3.5 so the largest disc
  # (KRAS G12D, n ~ 233) stops getting horizontally clipped at the
  # right edge of the plot region. The smaller cap also keeps the
  # rare-variant discs (n = 1) visually distinct without compressing
  # the dynamic range.
  scale_size_area(max_size = 3.5,
                  name   = "Patients (n)",
                  breaks = c(1, 10, 50, 100, 200),
                  guide  = guide_legend(override.aes = list(colour = "grey40"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Number of Patients", y = NULL
  ) +
  theme_avm() +
  theme(panel.grid.major.y = element_blank())

save_panel(panel_slot_dir("variant_landscape"), "variant_landscape", fig1b, 8, 5)

cat("\n── Variant counts ──\n")
print(table(landscape_data$variant_label))

# ── 4. Fig 1C: VAF violin + box plots by genotype ───────────────────────────

# VAF is only meaningful for variant-positive patients
vaf_data <- genotyped %>%
  filter(mutation_positive == TRUE & !is.na(vaf_prop)) %>%
  mutate(vaf_pct = vaf_prop * 100)

cat(sprintf("\n── VAF analysis: %d variant-positive patients with VAF ──\n", nrow(vaf_data)))

# Kruskal-Wallis test across genotype variants
kw_vaf <- kruskal.test(vaf_pct ~ geno_variant, data = vaf_data)
cat(sprintf("Kruskal-Wallis: chi-sq=%.2f, df=%d, p=%.4f\n",
  kw_vaf$statistic, kw_vaf$parameter, kw_vaf$p.value))

# Pairwise Wilcoxon with BH correction
pw_vaf <- pairwise.wilcox.test(vaf_data$vaf_pct, vaf_data$geno_variant, p.adjust.method = "BH")
cat("\nPairwise Wilcoxon (BH-corrected):\n")
print(pw_vaf$p.value, digits = 3)

fig1c <- ggplot(vaf_data, aes(x = geno_variant, y = vaf_pct, fill = geno_variant)) +
  # v6.7 (2026-05-20): violin fill alpha 0.5 -> 0.2 mirrors the density-
  # fill / CI-ribbon opacity used in Panels E and G so every "background
  # tint" reads at the same translucency across Fig 1. linewidth 0.5
  # unifies with KM curves (D/F/H), density curves (E), regression
  # lines (G).
  # v6.10 (2026-05-20): violin outer linewidth 0.5 -> 0.4 (20% lighter)
  # so the contour reads as a subtle shape outline rather than a
  # competing line element next to the boxplot.
  geom_violin(alpha = 0.2, scale = "width", trim = FALSE, linewidth = 0.4) +
  # v6.8 (2026-05-20): boxplot linewidth 0.5 -> 0.3 — the box+whisker
  # strokes at 0.5 read heavier than the surrounding violin contour
  # and the line geoms elsewhere in Fig 1. 0.3 lets the box read as a
  # quartile marker rather than a competing line element.
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, linewidth = 0.3) +
  # v6.3 (2026-05-21): geom_jitter dropped entirely. Even at size 0.6 /
  # alpha 0.4 the dense KRAS G12D column drew the eye more than the
  # violin shapes themselves; the distribution is fully described by
  # violin + boxplot, and per-stratum n moves to the figure caption.
  scale_fill_manual(values = PAL_VARIANT, guide = "none") +
  # v6.2: wrap two-word x-tick labels (e.g. "KRAS G12D" -> "KRAS\nG12D")
  # so the labels render horizontally without overlapping in the narrow
  # native cell.
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    x = NULL, y = "VAF (%)"
  ) +
  # v6.7 (2026-05-20): 5-tick y-axis (0/2/4/6/8) with .00 formatter so
  # VAF y-axes in Panel C and Panel G match each other and have the
  # same tick count as the x-axes. limits hard-set so the top tick
  # renders even when the empirical VAF max is < 8 %.
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     minor_breaks = NULL,
                     labels = scales::label_number(accuracy = 0.01),
                     limits = c(0, 8),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_avm() +
  theme(axis.text.x = element_text(hjust = 0.5),
        # Extra bottom margin so the wrapped 2-line tick labels do not
        # get clipped at the panel edge.
        plot.margin = margin(t = 4, r = 4, b = 8, l = 4))

save_panel(panel_slot_dir("vaf_by_variant"), "vaf_by_variant", fig1c, 8, 6)

# ── 5. Extended Data Fig. 2: Mutation spectrum by study ────────────────────────────────

# v6.40 (2026-05-27, Andy): order panels (a)/(b) by descending case count so
# the largest series sits at the top (panel a) / left (panel b). Counts:
# BCH 110, CHOP 105, UAB 72, Nikolaev 72, Gao 56, Priemer 21, Hong 21, Goss 16
# (ties UAB/Nikolaev and Priemer/Hong broken to match the published order).
# Local relevel only — the global study_clean order (01_clean_master.R) and
# PAL_STUDY keys are shared by Fig 1 and panel (c) and stay untouched.
ED_STUDY_ORDER <- c("BCH", "CHOP", "UAB", "Nikolaev", "Gao", "Priemer", "Hong", "Goss")

study_spectrum <- genotyped %>%
  mutate(study_clean = fct_relevel(study_clean, ED_STUDY_ORDER)) %>%
  count(study_clean, geno_variant) %>%
  group_by(study_clean) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# Statistical test: are mutation proportions different across studies?
set.seed(MASTER_SEED)   # audit F13: pin Monte-Carlo chi-square P
chisq_study <- chisq.test(table(genotyped$study_clean, genotyped$geno_variant),
  simulate.p.value = TRUE, B = 10000)
cat(sprintf("\n── Extended Data Fig. 2: Study x genotype chi-sq p=%.4f ──\n", chisq_study$p.value))

# Pairwise: compare each study's mut+ rate to the overall rate.
# Drop rows where study_clean is missing — the §4 prose cites the
# range over labelled contributing series, and a single NA-study row
# would be surfaced as "100%" (n=1) inflating the published range
# upper bound (audit 02).
mut_rate_by_study <- genotyped %>%
  filter(!is.na(study_clean)) %>%
  group_by(study_clean) %>%
  summarise(
    n = n(),
    n_mut = sum(mutation_positive == TRUE),
    mut_rate = 100 * n_mut / n,
    .groups = "drop"
  )
overall_rate <- 100 * sum(genotyped$mutation_positive) / nrow(genotyped)
cat(sprintf("Overall variant-positive rate: %.1f%%\n", overall_rate))
cat("Per-study rates:\n")
print(mut_rate_by_study)

# JAMA-style horizontal stacked bars by study
# relabel_geno_factor() swaps the raw "Negative" level to the canonical
# "Panel-negative" display label so the stacked-bar legend matches the
# rest of the manuscript. PAL_VARIANT is dual-keyed for both strings.
study_jama <- study_spectrum %>%
  rename(group = study_clean, category = geno_variant) %>%
  mutate(category = relabel_geno_factor(category))

n_study <- genotyped %>% count(study_clean) %>% deframe()

efig2 <- jama_stacked_bar(
  study_jama, palette = PAL_VARIANT, n_per_group = n_study,
  # v6.26 (2026-05-20): bar_height 0.3 -> 0.42 so inline % annotations
  # no longer clip at the bar top/bottom edges at the Nature-spec
  # 6.60 x 1.62 in cell.
  legend_title = "Genotype", bar_height = 0.42,
  # v6.21 (2026-05-20): base_size 14 -> 7 for Nature double-col native
  # rendering. Forces inline %-labels and tick text to scale into the
  # 5-7 pt band at the new 6.60 x 2.03 in cell size.
  base_size = 7,
  # v6.25 (2026-05-20): drop the dotted inter-row connection segments.
  # With 8 studies x 5 categories the dotted scaffolding read as visual
  # noise on top of an already busy stacked bar.
  show_guide_segments = FALSE
) +
  labs(title = NULL) +
  # v6.21: theme_nature_panel applied after jama_stacked_bar so the
  # in-helper theme_avm gets overridden. .no_grid blanks the grey93
  # gridlines that theme_nature_panel restores.
  theme_nature_panel() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "top")

# ── 6. Extended Data Fig. 2B: VAF distribution by study ────────────────────────────────

vaf_by_study <- genotyped %>%
  filter(!is.na(vaf_prop)) %>%
  # v6.40: same descending-case-count order as panel (a) (see ED_STUDY_ORDER).
  mutate(study_clean = fct_relevel(study_clean, ED_STUDY_ORDER),
         vaf_pct = vaf_prop * 100)

# Kruskal-Wallis for study effect on VAF
kw_study <- kruskal.test(vaf_pct ~ study_clean, data = vaf_by_study)
cat(sprintf("\n── VAF by study: Kruskal-Wallis p=%.4f ──\n", kw_study$p.value))

efig2b <- ggplot(vaf_by_study, aes(x = study_clean, y = vaf_pct, fill = study_clean)) +
  # v6.21 (2026-05-20): Fig 1 polish 1-to-1. Violin + slim boxplot, no
  # jitter dots (matches Fig 1 Panel C decision: dense cohorts drew the
  # eye more than the violin shape). alpha 0.2 fill + linewidth 0.5
  # outer / 0.3 box matches the Fig 1 background-tint convention.
  geom_violin(alpha = 0.2, scale = "width", trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, linewidth = 0.3) +
  labs(
    x = "Study", y = "VAF (%)"
  ) +
  scale_fill_manual(values = PAL_STUDY, guide = "none") +
  # v6.21: explicit 5-tick scale matching Fig 1 Panels C and G so the
  # VAF axes read identically across the manuscript.
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     minor_breaks = NULL,
                     labels = scales::label_number(accuracy = 0.01),
                     limits = c(0, 8),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_nature_panel() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ── Extended Data Fig. 2 combined: panels A + B ─────────────────────────────────────────
# Per-series mutation rate and VAF panels are panels (a) and (b) of
# ed_cohort_heterogeneity (registry token), so they land directly in the
# canonical dir alongside the meta-analysis panels written by 13_F2_genotype_phenotype.R
# and 09_F1_km_age.R. No legacy ExtendedData02/ stub is created.
ecoh_dir <- file.path(efig_dir, "ed_cohort_heterogeneity")
# v6.21 (2026-05-20): per-panel save dims sized for the 1-col x 3-row
# native ED Fig 4 layout: inner cell = 6.60 x 2.03 in (composite 7.20 x
# 6.69 with 0.3 in outer margins).
save_panel(ecoh_dir, "per_series_rate", efig2,  6.60, 2.03)
save_panel(ecoh_dir, "per_series_vaf",  efig2b, 6.60, 2.03)

# ── 7. Save legacy fig1_stats.rds (existing consumers) ──────────────────────

fig1_stats <- list(
  vaf_kw               = kw_vaf,
  vaf_pairwise         = pw_vaf,
  vaf_by_study_kw      = kw_study,
  chisq_study_genotype = chisq_study,
  mut_rate_by_study    = mut_rate_by_study,
  variant_counts       = table(landscape_data$variant_label),
  vaf_summary = vaf_data %>%
    group_by(geno_variant) %>%
    summarise(
      n          = n(),
      median_vaf = median(vaf_pct),
      q25        = quantile(vaf_pct, 0.25),
      q75        = quantile(vaf_pct, 0.75),
      .groups    = "drop"
    )
)

saveRDS(fig1_stats, file.path(output_dir, "stats", "fig1_stats.rds"))

# ── SupplementaryTables ──────────────────────────────────────────────────────────────────
library(writexl)
source(here("analysis", "helper_scripts", "supp_table_writer.R"))
etable_dir <- file.path(output_dir, "SupplementaryTables")
dir.create(etable_dir, recursive = TRUE, showWarnings = FALSE)
# Supplementary Table 01 — Variant-positive rate per contributing series
#   (accompanies Extended Data Fig. 3a in the narrative).
# Nature-Medicine headers: italic N / N (mutant); plain "Study" and
# "Mutation rate (%)". write_supp_table() hard-errors on NA / empty /
# wrong-type so a silently-dropped study can never ship to reviewers.
write_supp_table(
  data    = mut_rate_by_study,
  path    = file.path(etable_dir, "SuppTable01_mutation_rate_by_study.xlsx"),
  sheet   = "Mutation rate by study",
  columns = list(
    col    ("study_clean", label = "Study"),
    col_int("n",           label = "N",          italic = TRUE),
    col_int("n_mut",       label = "N (mutant)", italic = TRUE),
    col_num("mut_rate",    label = "Mutation rate (%)", digits = 1)
  )
)
# Supplementary Table 03 — Per-variant frequency (codon-level) + per-genotype
#   VAF summary among variant-positive lesions (accompanies Fig. 1b/1c and the
#   §1 rare-variant + VAF paragraphs). Sheet 1 carries the codon-level frequency
#   breakdown that defends the "rare variants < 3%" claim cited alongside
#   @supptab[variant_freq]; sheet 2 is the per-genotype VAF summary.
# Single sheet (2026-06-14): per-variant row count + VAF stats collapsed
# from the prior 2-sheet (variant_frequency + vaf_summary) workbook under
# the standing "one sheet per ST" rule. vaf_data already carries the
# per-sample VAF and the mut_* dummies, so we re-derive variant_label on
# vaf_data to get per-variant (not per-genotype-group) VAF stats — what
# the supplement should report.
vaf_by_variant <- vaf_data %>%
  mutate(
    variant_label = case_when(
      mut_KRAS_G12D        == 1 ~ "KRAS G12D",
      mut_KRAS_G12V        == 1 ~ "KRAS G12V",
      mut_KRAS_G12C        == 1 ~ "KRAS G12C",
      mut_KRAS_G12A        == 1 ~ "KRAS G12A",
      mut_KRAS_Q61H        == 1 ~ "KRAS Q61H",
      mut_KRAS_c191_196dup == 1 ~ "KRAS dup",
      mut_BRAF_V600E       == 1 ~ "BRAF V600E",
      mut_BRAF_Q636X       == 1 ~ "BRAF Q636X",
      TRUE                      ~ NA_character_
    )
  ) %>%
  group_by(variant_label) %>%
  summarise(
    n_vaf      = n(),
    median_vaf = median(vaf_pct),
    q25        = quantile(vaf_pct, 0.25),
    q75        = quantile(vaf_pct, 0.75),
    .groups    = "drop"
  )

# 2026-06-15: ST3 back to 2 sheets (Nature MOESM3 (A)/(B) pattern).
# (A) variant frequency over the full cohort (incl. Panel-negative);
# (B) VAF summary over variant-positive rows that have a measured VAF.
st3_freq_df <- variant_counts %>%
  arrange(desc(n)) %>%
  transmute(
    variant_label = as.character(variant_label),
    n_cohort      = as.integer(n),
    pct_cohort    = round(pct, 1)
  )

st3_vaf_df <- variant_counts %>%
  arrange(desc(n)) %>%
  left_join(vaf_by_variant, by = "variant_label") %>%
  filter(!variant_label %in% c("Panel-negative", "No tissue")) %>%
  mutate(
    variant_label  = as.character(variant_label),
    n_vaf_str = dplyr::case_when(
      is.na(n_vaf) ~ "n/a",
      TRUE         ~ formatC(n_vaf, format = "d")
    ),
    median_vaf_str = dplyr::case_when(
      is.na(median_vaf) ~ "n/a",
      TRUE              ~ sprintf("%.2f", median_vaf)
    ),
    iqr_str = dplyr::case_when(
      is.na(n_vaf) ~ "n/a",
      n_vaf < 3    ~ "n<3",
      TRUE         ~ sprintf("%.2f–%.2f", q25, q75)
    )
  ) %>%
  select(variant_label, n_vaf_str, median_vaf_str, iqr_str)

st3_sheets <- list(
  "(A) Variant frequency" = list(
    data = st3_freq_df,
    columns = list(
      col    ("variant_label", label = "Variant"),
      col_int("n_cohort",      label = "N (cohort)",  italic = TRUE),
      col_num("pct_cohort",    label = "% of cohort", digits = 1)
    ),
    footnote = "Variant labels are codon-level (KRAS G12D/V/C/A/Q61H/c.191_196dup; BRAF V600E/Q636X) plus Panel-negative (no variant detected by multiplex dPCR)."
  ),
  "(B) VAF summary" = list(
    data = st3_vaf_df,
    columns = list(
      col    ("variant_label",  label = "Variant"),
      col    ("n_vaf_str",      label = "N (VAF measured)", italic = TRUE),
      col    ("median_vaf_str", label = "Median VAF (%)"),
      col    ("iqr_str",        label = "IQR (Q25–Q75) (%)")
    ),
    footnote = c(
      "Per-variant VAF stats over variant-positive samples with a measured VAF. The Priemer published series contributes 5 KRAS G12V + 1 KRAS G12C without an extractable VAF; overall 314/320 (98.1%) of variant-positive patients have a measured VAF.",
      "'n/a' = no measured VAF for that variant. 'n<3' = too few samples for a Q25–Q75 interquartile range."
    )
  )
)
write_supp_table_workbook(
  sheets = st3_sheets,
  path   = file.path(etable_dir, "SuppTable03_variant_freq.xlsx")
)
cat(sprintf("  ✓ SuppTable03_variant_freq.xlsx (2 sheets: %d + %d rows)\n",
            nrow(st3_freq_df), nrow(st3_vaf_df)))
# Retained for any downstream consumers that reference st3_df.
st3_df <- bind_rows(st3_freq_df, st3_vaf_df)
cat("── Supplementary Tables 01\u201302 saved to results/SupplementaryTables/ ──\n")

# ═════════════════════════════════════════════════════════════════════════════
# ── 8. EXPANSIVE fig1 manifest fragment ─────────────────────────────────────
# ═════════════════════════════════════════════════════════════════════════════
# Every number that Results Section 1 ("Cohort assembly and somatic variant
# landscape") could plausibly reference. Key naming: <scope>_<quantity>.
# Scopes: cohort_, test_, variant_, vaf_, study_, sample_, rate_.
# All keys are scalars or short vectors — prose pulls them via
# `r stats$fig1$<key>`.

# ---- cohort denominators ---------------------------------------------------
cohort_n_total     <- nrow(df)
cohort_n_tissue    <- nrow(genotyped)
cohort_n_not_tested <- cohort_n_total - cohort_n_tissue
cohort_n_mut_pos   <- sum(genotyped$mutation_positive == TRUE,  na.rm = TRUE)
cohort_n_mut_neg   <- sum(genotyped$mutation_positive == FALSE, na.rm = TRUE)
mut_rate_pct       <- 100 * cohort_n_mut_pos / cohort_n_tissue
tissue_pct         <- 100 * cohort_n_tissue  / cohort_n_total

# ---- per-variant counts and proportions ------------------------------------
# (based on geno_variant: KRAS G12D / G12V / Other KRAS / BRAF / Negative)
gv <- genotyped$geno_variant
gv_tbl <- table(gv, useNA = "no")
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) == 1L && is.na(x))   return(y)
  x
}
get_n <- function(lvl) as.integer(gv_tbl[lvl] %||% 0L)

n_g12d       <- get_n("KRAS G12D")
n_g12v       <- get_n("KRAS G12V")
n_other_kras <- get_n("Other KRAS")
n_braf       <- get_n("BRAF")
n_neg        <- get_n("Negative")

# variant counts from variant-positive subset (fine-grained KRAS/BRAF breakout)
mut_pos_df <- genotyped %>% filter(mutation_positive == TRUE)

count_mut <- function(col) sum(mut_pos_df[[col]] == 1, na.rm = TRUE)
n_kras_g12d_fine <- count_mut("mut_KRAS_G12D")
n_kras_g12v_fine <- count_mut("mut_KRAS_G12V")
n_kras_g12c_fine <- count_mut("mut_KRAS_G12C")
n_kras_g12a_fine <- count_mut("mut_KRAS_G12A")
n_kras_q61h_fine <- count_mut("mut_KRAS_Q61H")
n_kras_dup_fine  <- count_mut("mut_KRAS_c191_196dup")
n_braf_v600e     <- count_mut("mut_BRAF_V600E")
n_braf_q636x     <- count_mut("mut_BRAF_Q636X")

n_kras_any <- n_kras_g12d_fine + n_kras_g12v_fine + n_kras_g12c_fine +
              n_kras_g12a_fine + n_kras_q61h_fine + n_kras_dup_fine
n_braf_any <- n_braf_v600e + n_braf_q636x

# KRAS codon 12: G12D + G12V + G12C + G12A + c.191_196dup (all map to codon 12)
n_kras_codon12         <- n_kras_g12d_fine + n_kras_g12v_fine +
                          n_kras_g12c_fine + n_kras_g12a_fine + n_kras_dup_fine
pct_kras_codon12_of_mutpos <- 100 * n_kras_codon12 / max(cohort_n_mut_pos, 1L)

# ---- proportions among variant-positive -----------------------------------
pct_of_mutpos <- function(k) 100 * k / max(cohort_n_mut_pos, 1L)
pct_g12d_of_mutpos      <- pct_of_mutpos(n_kras_g12d_fine)
pct_g12v_of_mutpos      <- pct_of_mutpos(n_kras_g12v_fine)
pct_otherkras_of_mutpos <- pct_of_mutpos(n_kras_g12c_fine + n_kras_g12a_fine +
                                         n_kras_q61h_fine + n_kras_dup_fine)
pct_kras_of_mutpos      <- pct_of_mutpos(n_kras_any)
pct_braf_of_mutpos      <- pct_of_mutpos(n_braf_any)

# ---- VAF summary stats (overall + per-variant) -----------------------------
vaf_all <- vaf_data$vaf_pct
vaf_overall <- list(
  n      = length(vaf_all),
  median = median(vaf_all),
  q25    = unname(quantile(vaf_all, 0.25)),
  q75    = unname(quantile(vaf_all, 0.75)),
  min    = min(vaf_all),
  max    = max(vaf_all),
  mean   = mean(vaf_all),
  sd     = sd(vaf_all)
)

vaf_per_variant <- vaf_data %>%
  group_by(geno_variant) %>%
  summarise(
    n      = n(),
    median = median(vaf_pct),
    q25    = quantile(vaf_pct, 0.25),
    q75    = quantile(vaf_pct, 0.75),
    min    = min(vaf_pct),
    max    = max(vaf_pct),
    mean   = mean(vaf_pct),
    sd     = sd(vaf_pct),
    .groups = "drop"
  )

# pairwise pvals extracted as named pairs
pw_mat <- pw_vaf$p.value
safe_pval <- function(row, col) {
  v <- tryCatch(pw_mat[row, col], error = function(e) NA_real_)
  if (is.null(v) || length(v) == 0L) NA_real_ else v
}

# ---- per-study mutation rate (sorted) --------------------------------------
mut_rate_by_study_tbl <- mut_rate_by_study %>% arrange(desc(mut_rate))

# ---- per-sample-type breakdown (FFPE / Fresh / Literature) ------------------
# sample_type levels: 1=FFPE, 2=Fresh/Frozen, 3=Literature (based on Fig1A DOT)
sample_type_tbl <- genotyped %>%
  mutate(sample_label = case_when(
    sample_type == 1 ~ "FFPE",
    sample_type == 2 ~ "Fresh/Frozen",
    sample_type == 3 ~ "Literature",
    TRUE             ~ "Unknown"
  )) %>%
  count(sample_label) %>%
  arrange(desc(n))

n_ffpe <- as.integer(sample_type_tbl$n[sample_type_tbl$sample_label == "FFPE"] %||% 0L)
n_fresh <- as.integer(sample_type_tbl$n[sample_type_tbl$sample_label == "Fresh/Frozen"] %||% 0L)
n_literature <- as.integer(sample_type_tbl$n[sample_type_tbl$sample_label == "Literature"] %||% 0L)
n_institutional <- n_ffpe + n_fresh  # = BCH+CHOP+UAB tissue-genotyped

# ---- fragment assembly -----------------------------------------------------

fig1_fragment <- list(

  # ── cohort denominators ───────────────────────────────────────────────
  cohort_n_total                = cohort_n_total,
  cohort_n_tissue_tested        = cohort_n_tissue,
  cohort_n_not_tested           = cohort_n_not_tested,
  cohort_n_mut_pos              = cohort_n_mut_pos,
  cohort_n_mut_neg              = cohort_n_mut_neg,
  cohort_pct_tissue_tested      = tissue_pct,
  cohort_mut_rate_pct           = mut_rate_pct,

  # ── variant counts (fine-grained, variant-positive subset) ───────────
  n_kras_g12d                   = n_kras_g12d_fine,
  n_kras_g12v                   = n_kras_g12v_fine,
  n_kras_g12c                   = n_kras_g12c_fine,
  n_kras_g12a                   = n_kras_g12a_fine,
  n_kras_q61h                   = n_kras_q61h_fine,
  n_kras_dup                    = n_kras_dup_fine,
  n_braf_v600e                  = n_braf_v600e,
  n_braf_q636x                  = n_braf_q636x,
  n_kras_any                    = n_kras_any,
  n_braf_any                    = n_braf_any,
  n_kras_codon12                = n_kras_codon12,
  pct_kras_codon12_of_mutpos    = pct_kras_codon12_of_mutpos,
  # institutional vs. published provenance
  n_institutional_cohort        = cohort_n_total -
                                    (sum(df$study_clean %in%
                                         c("Nikolaev", "Priemer", "Hong",
                                           "Goss", "Gao"))),
  n_published_cohort            = sum(df$study_clean %in%
                                      c("Nikolaev", "Priemer", "Hong",
                                        "Goss", "Gao")),
  # ── variant counts (geno_variant factor, 5 levels) ────────────────────
  n_g12d_grp                    = n_g12d,
  n_g12v_grp                    = n_g12v,
  n_other_kras_grp              = n_other_kras,
  n_braf_grp                    = n_braf,
  n_neg_grp                     = n_neg,

  # ── variant proportions among variant-positive ───────────────────────
  pct_g12d_of_mutpos            = pct_g12d_of_mutpos,
  pct_g12v_of_mutpos            = pct_g12v_of_mutpos,
  pct_otherkras_of_mutpos       = pct_otherkras_of_mutpos,
  pct_kras_of_mutpos            = pct_kras_of_mutpos,
  pct_braf_of_mutpos            = pct_braf_of_mutpos,

  # ── sample-type breakdown ─────────────────────────────────────────────
  n_ffpe                        = n_ffpe,
  n_fresh                       = n_fresh,
  n_literature                  = n_literature,
  n_institutional_tissue        = n_institutional,

  # ── VAF — overall ─────────────────────────────────────────────────────
  vaf_n                         = vaf_overall$n,
  vaf_median                    = vaf_overall$median,
  vaf_q25                       = vaf_overall$q25,
  vaf_q75                       = vaf_overall$q75,
  vaf_min                       = vaf_overall$min,
  vaf_max                       = vaf_overall$max,
  vaf_mean                      = vaf_overall$mean,
  vaf_sd                        = vaf_overall$sd,

  # ── VAF — per-variant medians/IQR (most-cited prose numbers) ──────────
  vaf_median_g12d = vaf_per_variant$median[vaf_per_variant$geno_variant == "KRAS G12D"]       %||% NA_real_,
  vaf_q25_g12d    = vaf_per_variant$q25[vaf_per_variant$geno_variant    == "KRAS G12D"]       %||% NA_real_,
  vaf_q75_g12d    = vaf_per_variant$q75[vaf_per_variant$geno_variant    == "KRAS G12D"]       %||% NA_real_,
  vaf_median_g12v = vaf_per_variant$median[vaf_per_variant$geno_variant == "KRAS G12V"]       %||% NA_real_,
  vaf_q25_g12v    = vaf_per_variant$q25[vaf_per_variant$geno_variant    == "KRAS G12V"]       %||% NA_real_,
  vaf_q75_g12v    = vaf_per_variant$q75[vaf_per_variant$geno_variant    == "KRAS G12V"]       %||% NA_real_,
  vaf_median_otherkras = vaf_per_variant$median[vaf_per_variant$geno_variant == "Other KRAS"] %||% NA_real_,
  vaf_median_braf      = vaf_per_variant$median[vaf_per_variant$geno_variant == "BRAF"]       %||% NA_real_,
  # per-variant VAF-measured n (for figure legends)
  vaf_n_g12d      = vaf_per_variant$n[vaf_per_variant$geno_variant == "KRAS G12D"]   %||% 0L,
  vaf_n_g12v      = vaf_per_variant$n[vaf_per_variant$geno_variant == "KRAS G12V"]   %||% 0L,
  vaf_n_otherkras = vaf_per_variant$n[vaf_per_variant$geno_variant == "Other KRAS"]  %||% 0L,
  vaf_n_braf      = vaf_per_variant$n[vaf_per_variant$geno_variant == "BRAF"]        %||% 0L,

  # ── VAF — omnibus + pairwise tests ────────────────────────────────────
  vaf_kw_stat                   = unname(kw_vaf$statistic),
  vaf_kw_df                     = unname(kw_vaf$parameter),
  vaf_kw_p                      = kw_vaf$p.value,
  vaf_pw_g12d_vs_g12v           = safe_pval("KRAS G12V",  "KRAS G12D"),
  vaf_pw_g12d_vs_otherkras      = safe_pval("Other KRAS", "KRAS G12D"),
  vaf_pw_g12v_vs_otherkras      = safe_pval("Other KRAS", "KRAS G12V"),

  # ── per-study statistics ──────────────────────────────────────────────
  chisq_study_genotype_stat     = unname(chisq_study$statistic),
  chisq_study_genotype_p        = chisq_study$p.value,
  vaf_by_study_kw_stat          = unname(kw_study$statistic),
  vaf_by_study_kw_p             = kw_study$p.value,

  # ── full tables (for SupplementaryTables + downstream consumers) ──────────────────
  mut_rate_by_study             = mut_rate_by_study_tbl,
  vaf_per_variant               = vaf_per_variant,
  sample_type_breakdown         = sample_type_tbl,
  variant_counts_full           = as.data.frame(gv_tbl,
                                                stringsAsFactors = FALSE,
                                                responseName     = "n")
)

# ---- emit manifest fragment ------------------------------------------------
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(section = "fig1", stats = fig1_fragment)

# ---- efig02 fragment (study-level heterogeneity) ---------------------------
# Extended Data Fig. 02A: mut+ rate by study (companion bar chart)
# Extended Data Fig. 02B: VAF by study (companion box/jitter)
# Stats overlap with fig1 (we intentionally duplicate compact keys so
# Section-1 prose that specifically cites eFig02 reads naturally).

hi_row  <- mut_rate_by_study_tbl[1, , drop = FALSE]
lo_row  <- mut_rate_by_study_tbl[nrow(mut_rate_by_study_tbl), , drop = FALSE]

efig02_fragment <- list(
  # scope
  n_studies                      = nrow(mut_rate_by_study_tbl),
  n_tissue_tested                = cohort_n_tissue,
  n_mut_pos                      = cohort_n_mut_pos,
  pooled_mut_rate_pct            = mut_rate_pct,

  # panel A — mut+ rate by study
  highest_rate_study             = as.character(hi_row$study_clean),
  highest_rate_pct               = hi_row$mut_rate,
  highest_rate_n                 = hi_row$n,
  highest_rate_n_mut             = hi_row$n_mut,
  lowest_rate_study              = as.character(lo_row$study_clean),
  lowest_rate_pct                = lo_row$mut_rate,
  lowest_rate_n                  = lo_row$n,
  lowest_rate_n_mut              = lo_row$n_mut,
  chisq_study_genotype_stat      = unname(chisq_study$statistic),
  chisq_study_genotype_p         = chisq_study$p.value,
  mut_rate_by_study              = mut_rate_by_study_tbl,

  # panel B — VAF by study
  vaf_by_study_kw_stat           = unname(kw_study$statistic),
  vaf_by_study_kw_df             = unname(kw_study$parameter),
  vaf_by_study_kw_p              = kw_study$p.value,
  vaf_by_study_summary = vaf_data %>%
    dplyr::group_by(study_clean) %>%
    dplyr::summarise(
      n      = dplyr::n(),
      median = median(vaf_pct),
      q25    = quantile(vaf_pct, 0.25),
      q75    = quantile(vaf_pct, 0.75),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median))
)

write_stats_section(section = "edfig02", stats = efig02_fragment)

# ═════════════════════════════════════════════════════════════════════════════
# ── 9. Fig 1 composite PNG (patchwork) ──────────────────────────────────────
# ═════════════════════════════════════════════════════════════════════════════
# Layout:
#   ┌───────────┬───────────┐
#   │           │   Fig 1B  │   A (CONSORT flow) spans both rows of col 1
#   │  Fig 1A   ├───────────┤   B (variant landscape)     top-right
#   │ flowchart │   Fig 1C  │   C (VAF violins)           bottom-right
#   └───────────┴───────────┘
# A is loaded from the pre-rendered flowchart PNG
# (results/Figure1/Fig1A_consort_flow.png) via png::readPNG + rasterGrob.
# If the flowchart file is missing, falls back to B/C stacked.

# Composite assembly is now centralised in 26_F1_assemble.R
# The legacy inline
# block below is gated behind if(FALSE){...} as historical documentation
# of the layout intent; the canonical Fig 1 composite path stays at
# results/Figure1/Fig1_composite.{pdf,png}.
cat("\n── Fig 1 composite: deferred to 26_F1_assemble.R ──\n")


# ═════════════════════════════════════════════════════════════════════════════
# ── 10. Human-readable stats readout (copy-paste friendly) ──────────────────
# ═════════════════════════════════════════════════════════════════════════════
# A table-formatted dump of every prose-ready stat in fig1_fragment. Written
# both to stdout (visible in run_all.R log) and to
# results/stats/fig1_readout.txt for easy grep / copy into a draft.

readout_path <- file.path(output_dir, "stats", "fig1_readout.txt")
readout_con  <- file(readout_path, open = "w")

emit <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = readout_con)
}

emit("")
emit("════════════════════════════════════════════════════════════════════")
emit(" FIGURE 1 — human-readable stats readout")
emit(" Source: analysis/01_main_analysis/08_F1_cohort_variants.R")
emit(" Manifest fragment: results/stats/_manifest_fragments/fig1.rds")
emit(" Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
emit("════════════════════════════════════════════════════════════════════")
emit("")
emit("── COHORT DENOMINATORS ──────────────────────────────────────────────")
emit("  Total assembled cohort ............ n = %d",   cohort_n_total)
emit("  Tissue-genotyped subset ........... n = %d  (%.1f%% of total)",
     cohort_n_tissue, tissue_pct)
emit("  Not genotyped (tissue pending) .... n = %d", cohort_n_not_tested)
emit("  Variant-positive ................. n = %d  (%.1f%% of tissue-genotyped)",
     cohort_n_mut_pos, mut_rate_pct)
emit("  Mutation-negative ................. n = %d  (%.1f%%)",
     cohort_n_mut_neg, 100 - mut_rate_pct)
emit("")
emit("── SAMPLE TYPE (tissue-genotyped subset) ────────────────────────────")
emit("  FFPE .............................. n = %d", n_ffpe)
emit("  Fresh / frozen .................... n = %d", n_fresh)
emit("  Literature (published cohorts) .... n = %d", n_literature)
emit("  Institutional tissue (FFPE+fresh) . n = %d", n_institutional)
emit("")
emit("── VARIANT LANDSCAPE (variant-positive subset, n = %d) ─────────────",
     cohort_n_mut_pos)
emit("  KRAS G12D ......................... n = %3d (%.1f%% of mut+)",
     n_kras_g12d_fine, pct_g12d_of_mutpos)
emit("  KRAS G12V ......................... n = %3d (%.1f%%)",
     n_kras_g12v_fine, pct_g12v_of_mutpos)
emit("  KRAS G12C ......................... n = %3d", n_kras_g12c_fine)
emit("  KRAS G12A ......................... n = %3d", n_kras_g12a_fine)
emit("  KRAS Q61H ......................... n = %3d", n_kras_q61h_fine)
emit("  KRAS c.191_196dup ................. n = %3d", n_kras_dup_fine)
emit("    Σ KRAS (any allele) ............. n = %3d (%.1f%% of mut+)",
     n_kras_any, pct_kras_of_mutpos)
emit("  BRAF V600E ........................ n = %3d", n_braf_v600e)
emit("  BRAF Q636X ........................ n = %3d", n_braf_q636x)
emit("    Σ BRAF (any allele) ............. n = %3d (%.1f%% of mut+)",
     n_braf_any, pct_braf_of_mutpos)
emit("")
emit("── VAF (variant-positive with VAF measured, n = %d) ────────────────",
     vaf_overall$n)
emit("  Overall VAF (%%) ................... median %.2f  [IQR %.2f–%.2f]",
     vaf_overall$median, vaf_overall$q25, vaf_overall$q75)
emit("                                      range %.2f – %.2f", vaf_overall$min, vaf_overall$max)
emit("                                      mean %.2f (SD %.2f)",
     vaf_overall$mean, vaf_overall$sd)
emit("")
for (i in seq_len(nrow(vaf_per_variant))) {
  r <- vaf_per_variant[i, ]
  emit("  %-15s .............. n = %3d  median %.2f  [IQR %.2f–%.2f]",
       r$geno_variant, r$n, r$median, r$q25, r$q75)
}
emit("")
emit("  Kruskal–Wallis (VAF ~ geno_variant):")
emit("    χ² = %.3f  df = %d  p = %s",
     unname(kw_vaf$statistic), unname(kw_vaf$parameter),
     format.pval(kw_vaf$p.value, digits = 3))
emit("  Pairwise Wilcoxon (BH-adjusted):")
emit("    G12D vs. G12V ............. p = %s",
     format.pval(safe_pval("KRAS G12V",  "KRAS G12D"), digits = 3))
emit("    G12D vs. Other KRAS ....... p = %s",
     format.pval(safe_pval("Other KRAS", "KRAS G12D"), digits = 3))
emit("    G12V vs. Other KRAS ....... p = %s",
     format.pval(safe_pval("Other KRAS", "KRAS G12V"), digits = 3))
emit("")
emit("── PER-STUDY MUTATION RATE ──────────────────────────────────────────")
for (i in seq_len(nrow(mut_rate_by_study_tbl))) {
  r <- mut_rate_by_study_tbl[i, ]
  emit("  %-10s ............ %3d / %3d = %5.1f%% variant-positive",
       r$study_clean, r$n_mut, r$n, r$mut_rate)
}
emit("")
emit("  Overall (tissue-genotyped): %.1f%% variant-positive", mut_rate_pct)
emit("  χ² (study × geno_variant, sim p): χ² = %.2f  p = %s",
     unname(chisq_study$statistic),
     format.pval(chisq_study$p.value, digits = 3))
emit("  K–W (VAF by study):               χ² = %.2f  p = %s",
     unname(kw_study$statistic),
     format.pval(kw_study$p.value, digits = 3))
emit("")
emit("════════════════════════════════════════════════════════════════════")
emit("  Copy any of the above into prose. Every scalar is also available")
emit("  inline as  `r stats$fig1$<key>`  after run_all.R rebuilds the")
emit("  manifest. See analysis/pipeline/helpers/stats_schema.R for the full key list.")
emit("════════════════════════════════════════════════════════════════════")
close(readout_con)
cat(sprintf("  ✓ readout: %s\n", readout_path))

cat("\n══ 08_F1_cohort_variants.R complete ══\n")
