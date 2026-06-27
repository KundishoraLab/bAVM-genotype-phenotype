# 15_ED7_vaf_age_at_rupture.R
# ─────────────────────────────────────────────────────────────────────────────
# Mirror of Fig 2e (VAF × age-at-presentation) for age-at-RUPTURE: per
# Hale comment 2026-04-25, we asked whether the inverse VAF × age slope
# extends from presentation to the rupture endpoint. Restricted to
# Mut+ lesions that ruptured at presentation (Option-2a strict cohort
# matching Fig 2d KM curves), so the regression "age" inside this
# subset is age at rupture. Per-variant slopes (G12D, G12V) are
# reported for parity with Fig 2e even though their per-stratum n is
# small and CIs cross zero.
#
# Outputs:
#   - results/stats/_manifest_fragments/vaf_rupture.rds
#   - results/stats/vaf_age_at_rupture.txt
#   - results/ExtendedData/ed_vaf_age_rupture/vaf_age_rupture_scatter.{pdf,png,rds}
#   - results/SupplementaryTables/SuppTable05_vaf_age_rupture.xlsx
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(here); library(ggplot2); library(openxlsx)
})

# Match Fig 2e (vaf_age_scatter) palette + theme so the rupture-endpoint
# panel reads as a sibling figure when placed in the supplement.
# PAL_KM is sourced from utils.R below.
source(here("analysis", "helper_scripts", "utils.R"))

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  mutate(variant_group = case_when(
    !mutation_positive          ~ "Negative",
    mutation == "KRAS G12D"     ~ "KRAS G12D",
    mutation == "KRAS G12V"     ~ "KRAS G12V",
    mutation_gene == "KRAS"     ~ "Other KRAS",
    mutation_gene == "BRAF"     ~ "BRAF",
    TRUE                        ~ "Negative"
  ))

# Mut+ ruptured cohort: same Option-2a strict definition the rupture KM uses.
rupt_df <- df %>%
  filter(mutation_positive,
         rupture_category == "Ruptured at presentation",
         !is.na(vaf_prop), !is.na(age)) %>%
  mutate(vaf_pct = vaf_prop * 100)

slope_one <- function(d) {
  if (nrow(d) < 5L) {
    return(list(n = nrow(d), slope = NA_real_, ci_lo = NA_real_,
                ci_hi = NA_real_, p = NA_real_,
                rho = NA_real_, rho_p = NA_real_))
  }
  fit <- lm(age ~ vaf_pct, data = d)
  cf  <- summary(fit)$coefficients["vaf_pct", ]
  ci  <- confint(fit, "vaf_pct", level = 0.95)
  sp  <- suppressWarnings(cor.test(d$vaf_pct, d$age, method = "spearman"))
  list(
    n       = nrow(d),
    slope   = unname(cf["Estimate"]),
    ci_lo   = unname(ci[1]),
    ci_hi   = unname(ci[2]),
    p       = unname(cf["Pr(>|t|)"]),
    rho     = unname(sp$estimate),
    rho_p   = unname(sp$p.value)
  )
}

pooled <- slope_one(rupt_df)
g12d   <- slope_one(rupt_df %>% filter(variant_group == "KRAS G12D"))
g12v   <- slope_one(rupt_df %>% filter(variant_group == "KRAS G12V"))

# Side-by-side prevalence-vs-timing tabulation that supports the §2→§3
# transition sentence in Results.
# Canonical cohort labels — the resulting strings surface in
# `vaf_age_at_rupture.txt` (lines 118-121) and downstream SuppTable
# `slope_summary$cohort`, so use the full prose forms.
geno <- df %>%
  filter(!is.na(mutation_positive)) %>%
  mutate(geno = ifelse(mutation_positive, "Variant-positive", "Panel-negative"))

prev <- geno %>%
  group_by(geno) %>%
  summarise(
    n          = n(),
    n_rupt     = sum(rupture_category == "Ruptured at presentation",
                     na.rm = TRUE),
    pct_rupt   = 100 * n_rupt / n,
    .groups = "drop"
  )

age_rupt <- geno %>%
  filter(rupture_category == "Ruptured at presentation", !is.na(age)) %>%
  group_by(geno) %>%
  summarise(median_age_rupt = median(age), .groups = "drop")

age_pres <- geno %>%
  filter(!is.na(age)) %>%
  group_by(geno) %>%
  summarise(median_age_pres = median(age), .groups = "drop")

prev_rupt_mut    <- prev$pct_rupt[prev$geno == "Variant-positive"]
prev_rupt_neg    <- prev$pct_rupt[prev$geno == "Panel-negative"]
median_rupt_mut  <- age_rupt$median_age_rupt[age_rupt$geno == "Variant-positive"]
median_rupt_neg  <- age_rupt$median_age_rupt[age_rupt$geno == "Panel-negative"]
median_pres_mut  <- age_pres$median_age_pres[age_pres$geno == "Variant-positive"]
median_pres_neg  <- age_pres$median_age_pres[age_pres$geno == "Panel-negative"]

# ── readout ────────────────────────────────────────────────────────────────
out_dir <- here("results", "stats")
writeLines(c(
  sprintf("# %s — VAF × age-at-rupture (mirror of %s) + prevalence-vs-timing",
          panel_prose_tag("vaf_age_rupture_scatter"),
          panel_prose_tag("vaf_age_scatter")),
  sprintf("Generated %s", Sys.Date()),
  "",
  "## Q1: VAF × age-at-rupture (Variant-positive ruptured cohort)",
  sprintf("  Pooled (n=%d):  slope = %+.2f y/%%VAF (95%% CI %+.2f, %+.2f); P = %.3g; Spearman rho = %+.2f, P = %.3g",
          pooled$n, pooled$slope, pooled$ci_lo, pooled$ci_hi, pooled$p, pooled$rho, pooled$rho_p),
  sprintf("  KRAS G12D (n=%d):  slope = %+.2f y/%%VAF (95%% CI %+.2f, %+.2f); P = %.3g",
          g12d$n, g12d$slope, g12d$ci_lo, g12d$ci_hi, g12d$p),
  sprintf("  KRAS G12V (n=%d):  slope = %+.2f y/%%VAF (95%% CI %+.2f, %+.2f); P = %.3g",
          g12v$n, g12v$slope, g12v$ci_lo, g12v$ci_hi, g12v$p),
  "",
  "## Q2: prevalence-vs-timing side-by-side",
  sprintf("  Variant-positive: rupture prevalence = %.1f%%, median age at rupture = %.1f y, median age at presentation = %.1f y",
          prev_rupt_mut, median_rupt_mut, median_pres_mut),
  sprintf("  Panel-negative:    rupture prevalence = %.1f%%, median age at rupture = %.1f y, median age at presentation = %.1f y",
          prev_rupt_neg, median_rupt_neg, median_pres_neg),
  "",
  "Interpretation: similar prevalence (the WHO of rupture) coexists",
  "with a ~two-decade earlier median age at rupture (the WHEN). The",
  "rupture-prevalence forest in Fig 3d sums cumulative incidence at",
  "long follow-up; the rupture KM in Fig 2d shows the hazard shift."
), file.path(out_dir, "vaf_age_at_rupture.txt"))

# ── manifest fragment ──────────────────────────────────────────────────────
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(
  section = "vaf_rupture",
  stats = list(
    n_rupt_pooled       = pooled$n,
    slope_pooled        = pooled$slope,
    slope_ci_lo_pooled  = pooled$ci_lo,
    slope_ci_hi_pooled  = pooled$ci_hi,
    slope_p_pooled      = pooled$p,
    rho_pooled          = pooled$rho,
    rho_p_pooled        = pooled$rho_p,
    n_g12d              = g12d$n,
    slope_g12d          = g12d$slope,
    slope_g12d_p        = g12d$p,
    n_g12v              = g12v$n,
    slope_g12v          = g12v$slope,
    slope_g12v_p        = g12v$p,
    prev_rupt_mut_pct   = prev_rupt_mut,
    prev_rupt_neg_pct   = prev_rupt_neg,
    median_age_rupt_mut = median_rupt_mut,
    median_age_rupt_neg = median_rupt_neg,
    median_age_pres_mut = median_pres_mut,
    median_age_pres_neg = median_pres_neg
  )
)

cat(sprintf("✓ VAF × age-at-rupture: pooled slope %+.2f y/%%VAF (P = %.3g); prevalence Variant-positive %.1f%% vs Panel-negative %.1f%%\n",
            pooled$slope, pooled$p, prev_rupt_mut, prev_rupt_neg))

# ── Extended Data figure: VAF × age-at-rupture scatter ─────────────────────
# Mirrors the Fig 2e (vaf_age_scatter) helper line-for-line: same palette,
# same x-axis ticks at 2/4/6, same theme_avm. Restricted to G12D + G12V so
# the per-variant slopes are interpretable (same restriction Fig 2e uses).
plot_df <- rupt_df %>%
  filter(variant_group %in% c("KRAS G12D", "KRAS G12V")) %>%
  mutate(variant_group = factor(variant_group,
                                levels = c("KRAS G12D", "KRAS G12V")))

# v6.30 (2026-05-20): port to the canonical vaf_age_scatter_panel()
# helper (utils.R) so this panel is the inverted-axes sibling of Fig 1
# Panel G, sharing geom grammar bar-for-bar (alpha 0.3 / size 0.6 /
# linewidth 0.5 / smooth alpha 0.2).
p_rupt_scatter <- vaf_age_scatter_panel(
  plot_df,
  x_var     = "vaf_pct",
  y_var     = "age",
  color_var = "variant_group",
  palette   = PAL_KM,
  vaf_axis  = "x",
  x_lab     = "VAF (%)",
  y_lab     = "Age at rupture (years)"
) + theme(legend.position = "none")   # shared legend in the ED07 composer gutter

# 2026-05-19 (Phase 2 / Iteration 3): ed_vaf_age_rupture merged into the
# ed_vaf_deep_dive umbrella group. vaf_age_rupture_scatter writes flat
# at the group root so 31_ED7_assemble.R can read its rds.
ed_dir <- here("results", "ExtendedData", "ed_vaf_deep_dive")
dir.create(ed_dir, recursive = TRUE, showWarnings = FALSE)
# v6.29 (2026-05-20): write into the canonical panel_A/ slot (resolved
# via panel_slot_dir) so sync_panel_prefixes doesn't sweep the group-
# root copies as registry-misplaced. Composite-row-1 footprint at
# Nature double-col spec is 6.60 x 1.20 in.
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa_var <- load_panel_assignments()
ed_dir <- panel_slot_dir("vaf_age_rupture_scatter", .pa_var)
dir.create(ed_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(ed_dir, "vaf_age_rupture_scatter.pdf"), p_rupt_scatter,
       width = 6.60, height = 1.20, device = cairo_pdf,
       family = NM$font_family)
ggsave(file.path(ed_dir, "vaf_age_rupture_scatter.png"), p_rupt_scatter,
       width = 6.60, height = 1.20, dpi = 300, type = "cairo")
saveRDS(p_rupt_scatter, file.path(ed_dir, "vaf_age_rupture_scatter.rds"))

# ── Supplementary Table 05 — VAF × age-at-rupture slope summary ────────────
# Single sheet (2026-06-14): per-stratum slope summary only. The per_lesion
# sheet was dropped — its only column besides study/variant/vaf/age was a
# structurally-NA patient_id placeholder (PHI rule prohibits patient
# identifiers anyway). Reviewers who want to reproduce the regression can
# refit from data/processed/bAVM_analysis_ready.rds with the same
# (mutation_positive & ruptured-at-presentation) filter the producer uses.
#
# Per-variant Spearman ρ is now surfaced (was previously hardcoded NA for
# G12D and G12V even though slope_one() already computes it).
supp_dir <- here("results", "SupplementaryTables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
source(here("analysis", "helper_scripts", "supp_table_writer.R"))

# Nature-style p-value formatter: 3 dp for p ≥ 0.001, scientific (1 dp) below.
.fmt_p <- function(p) {
  ifelse(p < 0.001,
         formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

st5_df <- tibble::tibble(
  cohort        = c("Pooled variant-positive ruptured", "KRAS G12D", "KRAS G12V"),
  n             = as.integer(c(pooled$n, g12d$n, g12v$n)),
  slope_ci_str  = sprintf("%.2f (%.2f to %.2f)",
                          c(pooled$slope, g12d$slope, g12v$slope),
                          c(pooled$ci_lo, g12d$ci_lo, g12v$ci_lo),
                          c(pooled$ci_hi, g12d$ci_hi, g12v$ci_hi)),
  p_str         = .fmt_p(c(pooled$p, g12d$p, g12v$p)),
  rho_str       = sprintf("%.2f", c(pooled$rho, g12d$rho, g12v$rho)),
  rho_p_str     = .fmt_p(c(pooled$rho_p, g12d$rho_p, g12v$rho_p))
)

write_supp_table(
  data    = st5_df,
  path    = file.path(supp_dir, "SuppTable05_vaf_age_rupture.xlsx"),
  sheet   = "VAF × age-at-rupture slope",
  columns = list(
    col    ("cohort",       label = "Cohort"),
    col_int("n",            label = "N",                          italic = TRUE),
    col    ("slope_ci_str", label = "Slope, years per 1% VAF (95% CI)"),
    col    ("p_str",        label = "P (slope)",                  italic = TRUE),
    col    ("rho_str",      label = "Spearman ρ",                 italic = TRUE),
    col    ("rho_p_str",    label = "P (Spearman)",               italic = TRUE)
  ),
  footnote = c(
    "Slope and 95% CI from ordinary least squares fit of age at rupture ~ VAF (%). Restricted to variant-positive lesions with measured VAF that ruptured at presentation.",
    "Per-variant rows (KRAS G12D, KRAS G12V) are reported for parity with Fig. 2e even though per-stratum N is small and the slope CI crosses zero in both.",
    "Spearman ρ is included as a rank-based, distribution-free anchor against the OLS slope."
  )
)

cat(sprintf("✓ ED scatter + Supp Table 16 written\n"))
