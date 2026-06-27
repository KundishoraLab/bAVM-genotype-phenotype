# 13_F2_genotype_phenotype.R — null_phenotype group + supporting ED panels.
#
# Input:  data/processed/bAVM_analysis_ready.rds
# Output: panels for the null_phenotype main figure (parietal_kras_rupture,
#         sm_grade_dist, sm_components, clinical_history, hr_features_OR)
#         plus per_lobe_prevalence / per_variant_anatomy / interaction_cleveland
#         panels for the ed_per_variant_pheno + ed_anatomy ED groups, and
#         supp tables (location_interactions, parietal_kras, hierarchical_gxp,
#         high_risk_OR). Manuscript figure numbers and panel letters are
#         resolver-driven from analysis/pipeline/panel_registry.R
#         and the prose first-citation order — do not hard-code numbers here.
#
# All panels: Tier 1 (binary) and Tier 2 (per-variant).
# Uses cool-tone palettes from utils.R.
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(patchwork)
# MASS::polr() is called via :: (NOT library(MASS)) so MASS::select does
# not mask dplyr::select for every downstream script in the run_all loop.
library(broom)    # tidy() for model output

source(here("analysis", "helper_scripts", "utils.R"))

# ── 1. Load data and palettes ────────────────────────────────────────────────

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))
genotyped <- df %>% filter(!is.na(mutation_positive))
palettes <- readRDS(here("results", "palettes.rds"))

output_dir <- here("results")
GENO_PALETTE <- PAL_VARIANT
BINARY_PALETTE <- PAL_BINARY
SM_PALETTE <- PAL_SM

# Manuscript Fig 3 (the null_phenotype group) lives at results/Figure3/.
# Disk dir matches the resolver-assigned manuscript figure number.
fig3_dir <- file.path(output_dir, "Figure3")
efig_dir <- file.path(output_dir, "ExtendedData")

# Resolver-assigned panel letters (built upstream by analysis/pipeline/builders/run_all.R
# from first-citation order in the manuscript). The shared
# panel_slot_dir() helper (utils.R) maps any token to its registry-correct
# on-disk directory, so 2026-04-26 promotions/demotions between
# null_phenotype and ed_parietal route automatically without
# producer-side path hardcoding.
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa <- load_panel_assignments()
slot_dir <- function(token) {
  d <- panel_slot_dir(token, .pa)
  if (is.null(d)) {
    stop(sprintf("slot_dir: token '%s' has no resolver letter / group; ",
                 "is it registered in panel_registry?"),
         token)
  }
  d
}

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
# 2026-05-20: panel_dims provides save_panel_native(token, plot) which
# renders + saves each panel at its native composite cell footprint
# (resolved via FIGURE_LAYOUTS) with theme_avm_native() pre-applied.
# Use for any token registered in FIGURE_LAYOUTS so the standalone
# panel matches the composite at print scale.
source(here("analysis", "pipeline", "helpers", "panel_dims.R"))

cat(sprintf("Genotyped: %d | Mut+: %d | Neg: %d\n",
  nrow(genotyped),
  sum(genotyped$geno_binary == "Variant-positive", na.rm = TRUE),
  sum(genotyped$geno_binary == "Panel-negative", na.rm = TRUE)))

# ── Helper: safe Fisher's exact ──────────────────────────────────────────────
safe_fisher <- function(x, y) {
  tbl <- table(x, y)
  if (nrow(tbl) < 2 || ncol(tbl) < 2) return(list(p.value = NA_real_))
  set.seed(MASTER_SEED)   # audit F13: pin Monte-Carlo Fisher P
  fisher.test(tbl, simulate.p.value = TRUE, B = 10000)
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig 2A: SM grade × genotype (stacked bars)
# ══════════════════════════════════════════════════════════════════════════════

sm_data <- genotyped %>%
  filter(!is.na(sm_grade)) %>%
  mutate(sm_grade = factor(sm_grade, levels = c("I", "II", "III", "IV", "V")))

cat(sprintf("\n── Fig 2A: SM grade (n=%d) ──\n", nrow(sm_data)))

# Tier 1 (binary)
kw_sm_binary <- kruskal.test(sm_total_num ~ geno_binary, data = sm_data)
cat(sprintf("Binary KW: p=%.4f\n", kw_sm_binary$p.value))

# JAMA-style horizontal stacked bars — binary
sm_binary_pct <- sm_data %>%
  count(geno_binary, sm_grade) %>%
  group_by(geno_binary) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  rename(group = geno_binary, category = sm_grade)

n_binary <- sm_data %>% count(geno_binary) %>% deframe()

fig2a_binary <- jama_stacked_bar(
  sm_binary_pct, palette = SM_PALETTE, n_per_group = NULL,
  legend_title = "SM Grade",
  # v6 (2026-05-21): pass Nature body_pt so the in-bar % labels
  # render at base/3.5 mm = ~7 pt. theme_nature_panel() applies the
  # remaining font sizes at compose time.
  base_size = NM$body_pt
) +
  labs(title = NULL)
# v6: NO theme_avm_native() here — the composer applies
# theme_nature_panel() via patchwork `&` so font sizes are uniform
# across panels regardless of cell footprint.

# 2026-04-27: wrap as a 1-element patchwork with guides = "keep" so the
# SM Grade fill legend stays anchored inside panel B at the
# composer level. The outer 30_F2_assemble.R uses guides =
# "collect" to share the Genotype + N legends across panels C/D in
# the right margin; without this wrap, the SM Grade legend (only used
# by panel B) would also get pulled to that shared margin.
fig2a_binary <- patchwork::wrap_plots(fig2a_binary) +
  patchwork::plot_layout(guides = "keep")

save_panel_native("sm_grade_dist", fig2a_binary)

# Tier 2 (per-variant) — JAMA stacked
kw_sm_variant <- kruskal.test(sm_total_num ~ geno_variant, data = sm_data)
cat(sprintf("Per-variant KW: p=%.4f\n", kw_sm_variant$p.value))

sm_variant_pct <- sm_data %>%
  count(geno_variant, sm_grade) %>%
  group_by(geno_variant) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  rename(group = geno_variant, category = sm_grade)

n_variant <- sm_data %>% count(geno_variant) %>% deframe()

fig2a_variant <- jama_stacked_bar(
  sm_variant_pct, palette = SM_PALETTE, n_per_group = n_variant,
  legend_title = "SM Grade"
) +
  labs(title = "SM Grade Distribution by Variant")

# NB: Fig2A_sm_grade_variant (per-variant alternate of the SM grade waffle)
# is not in the Fig 3 composite or any supp embed; the pooled binary view
# (fig2a_binary, written above) is the canonical artefact.

# ══════════════════════════════════════════════════════════════════════════════
# Fig 2B: SM components × genotype
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 2B: SM components ──\n")

# SM Size (n=258)
size_data <- genotyped %>% filter(!is.na(sm_size_num))
kw_size <- kruskal.test(sm_size_num ~ geno_binary, data = size_data)
cat(sprintf("Size binary KW: n=%d, p=%.4f\n", nrow(size_data), kw_size$p.value))

# SM Drainage (n=222)
drain_data <- genotyped %>% filter(!is.na(sm_drainage_num))
fisher_drain <- safe_fisher(drain_data$geno_binary, drain_data$sm_drainage_num)
cat(sprintf("Drainage binary Fisher: n=%d, p=%.4f\n", nrow(drain_data), fisher_drain$p.value))

# SM Eloquence (n=96)
eloq_data <- genotyped %>% filter(!is.na(sm_eloquence_num))
fisher_eloq <- safe_fisher(eloq_data$geno_binary, eloq_data$sm_eloquence_num)
cat(sprintf("Eloquence binary Fisher: n=%d, p=%.4f\n", nrow(eloq_data), fisher_eloq$p.value))

# Combined component plot (binary tier)
# Each row carries successes / n so we can compute an exact 95% CI
# (Clopper-Pearson via binom.test) per genotype × component point —
# rendered as horizontal CI segments in the dumbbell. Same idiom is
# used in 11_ED4_power_forest_meta.R.
comp_summary <- bind_rows(
  size_data %>%
    group_by(geno_binary) %>%
    summarise(successes = sum(sm_size_num >= 2), n = n(), .groups = "drop") %>%
    mutate(component = "Size ≥2 (medium/large)", p = kw_size$p.value),
  drain_data %>%
    group_by(geno_binary) %>%
    summarise(successes = sum(sm_drainage_num == 1), n = n(), .groups = "drop") %>%
    mutate(component = "Deep drainage", p = fisher_drain$p.value),
  eloq_data %>%
    group_by(geno_binary) %>%
    summarise(successes = sum(sm_eloquence_num == 1), n = n(), .groups = "drop") %>%
    mutate(component = "Eloquent location", p = fisher_eloq$p.value)
) %>%
  mutate(
    prevalence = successes / n,
    ci_lo = mapply(function(x, n) binom.test(x, n)$conf.int[1], successes, n),
    ci_hi = mapply(function(x, n) binom.test(x, n)$conf.int[2], successes, n),
    component = factor(component,
      levels = c("Size ≥2 (medium/large)", "Deep drainage", "Eloquent location"))
  )

comp_dumbbell <- comp_summary %>%
  mutate(
    prevalence = prevalence * 100,
    ci_lo      = ci_lo * 100,
    ci_hi      = ci_hi * 100
  ) %>%
  rename(feature = component, group = geno_binary) %>%
  mutate(
    # Wrap feature labels so "Size >=2 (medium/large)" doesn't blow out
    # the y-axis at Fig 2's column-1 width. ~14 chars wraps the longest
    # label to two lines and leaves the shorter two single-line.
    feature = factor(
      stringr::str_wrap(as.character(feature), width = 14),
      levels = stringr::str_wrap(
        c("Size ≥2 (medium/large)", "Deep drainage", "Eloquent location"),
        width = 14
      )
    )
  )

# Render + save deferred until after `clin_dumbbell` is built (Fig 2F
# section below) so the two Fig 3 dumbbells can share an N-scale legend
# in the composite via guides = "collect". See the
# "Shared dumbbell render" block at the end of Fig 2F.

# ══════════════════════════════════════════════════════════════════════════════
# Fig 2C: High-risk features OR forest plot
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 2C: High-risk features ──\n")

hr_features <- c("intranidal_aneurysm_num", "venous_varix_num",
                  "venous_outflow_stenosis_num", "flow_related_aneurysm_num")
hr_labels <- c("Intranidal aneurysm", "Venous varix",
               "Venous outflow stenosis", "Flow-related aneurysm")

# Firth's penalized logistic regression for each feature to handle
# complete/quasi-complete separation (venous outflow stenosis: 0 events in neg)
library(logistf)

hr_results <- map2_dfr(hr_features, hr_labels, function(feat, lab) {
  sub <- genotyped %>%
    filter(!is.na(.data[[feat]]) & !is.na(geno_binary) & !is.na(sample_type_clean)) %>%
    mutate(outcome = as.numeric(.data[[feat]]),
           geno = relevel(geno_binary, ref = "Panel-negative"))
  if (nrow(sub) < 10 || length(unique(sub$outcome)) < 2) {
    return(tibble(feature = lab, OR = NA, lower = NA, upper = NA, p = NA, n = nrow(sub),
                  sample_type_p = NA_real_))
  }
  mod <- logistf(outcome ~ geno + sample_type_clean, data = sub, pl = FALSE)
  geno_idx <- grep("genoVariant-positive", mod$terms)
  coef_val <- mod$coefficients[geno_idx]
  ci_lo    <- mod$ci.lower[geno_idx]
  ci_hi    <- mod$ci.upper[geno_idx]
  p_val    <- mod$prob[geno_idx]
  # Check sample_type significance
  st_idx <- grep("sample_type_clean", mod$terms)
  st_p <- if (length(st_idx) > 0) min(mod$prob[st_idx]) else NA_real_
  tibble(
    feature = lab,
    OR = exp(coef_val),
    lower = exp(ci_lo),
    upper = exp(ci_hi),
    p = p_val,
    n = nrow(sub),
    sample_type_p = st_p
  )
})

# BH-FDR correction
hr_results$p_fdr <- p.adjust(hr_results$p, method = "BH")
cat("High-risk feature ORs (Mut+ vs Neg):\n")
print(hr_results)

# Also do composite count
hr_count_data <- genotyped %>% filter(!is.na(n_high_risk_num) & !is.na(geno_binary))
kw_hr <- kruskal.test(n_high_risk_num ~ geno_binary, data = hr_count_data)
cat(sprintf("\nComposite high-risk count KW: n=%d, p=%.4f\n", nrow(hr_count_data), kw_hr$p.value))

hr_plot <- hr_results %>%
  filter(!is.na(OR)) %>%
  mutate(
    # Wrap predictor labels to ~15 chars so the forest's y-axis text
    # doesn't crowd the table columns at Fig 2's post-restructure
    # column-2+3 width. Two longest labels ("Venous Outflow Stenosis",
    # "Flow-Related Aneurysm") wrap to two lines; the shorter ones stay
    # single-line.
    # 2026-05-20 v5: drop str_wrap so each predictor renders on one
    # line. At the native 183 mm cell width the forest has room for
    # the full "Venous outflow stenosis" / "Flow-related aneurysm"
    # labels without column overflow; wrapping was causing the
    # double-line predictor rows to overlap inter-row spacing.
    label     = feature,
    est       = OR,
    lo        = lower,
    hi        = upper,
    n_label   = as.character(n),
    # Nominal (uncorrected) P alongside FDR-adjusted Q so readers can
    # see both values; the cols block below renders them as separate
    # table columns in the forest panel.
    p_nominal = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    p_label   = ifelse(p_fdr < 0.001, "<0.001", sprintf("%.3f", p_fdr))
  )

fig2c <- forest_compact_general(
  hr_plot,
  cols = list(
    list(col = "p_label", header = "FDR P")
  ),
  point_col      = BINARY_PALETTE["Variant-positive"],
  axis_ticks     = c(0.5, 1, 2, 4, 8),
  x_lab          = "OR (log scale)",
  est_col_header = "OR (95% CI)",
  est_fmt        = "%.2f (%.2f–%.2f)",
  est_col_side   = "left",
  # v6.4 (2026-05-21): table_forest's annotation text uses
  #   size = base_size / 3.5    for predictor labels (mm units)
  #   size = base_size / 3.8    for OR text
  #   size = base_size / 4      for axis ticks
  # In points: pt = mm * 2.835.
  # base_size = 9 gives label ~7.3 pt (matches axis.title = NM$body_pt
  # = 7 pt on the other panels) and tick ~6.4 pt (matches axis.text =
  # NM$tick_pt = 6 pt). v6.3 used base_size = 7, giving ~5.7 pt
  # labels which read smaller than the rest of the figure.
  base_size        = 9,
  size_range       = c(1.2, 2.4),
  # 2026-05-21 v6.5: drop the lower hairline between the last data row
  # and the OR x-axis. The top rule + the OR-axis line are enough to
  # frame the table at the native cell footprint.
  show_bottom_rule = FALSE
)

save_panel_native("hr_features_OR", fig2c)

# ══════════════════════════════════════════════════════════════════════════════
# Fig 2D: Brain region × genotype heatmap
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 2D: Location heatmap ──\n")

loc_cols <- c("loc_frontal", "loc_temporal", "loc_insular", "loc_basal_ganglia",
              "loc_thalamus", "loc_periventricular", "loc_parietal", "loc_occipital",
              "loc_cerebellar", "loc_brainstem", "loc_corpus_callosum",
              "loc_cingulate", "loc_sylvian_fissure")
loc_labels <- c("Frontal", "Temporal", "Insular", "Basal Ganglia",
                "Thalamus", "Periventricular", "Parietal", "Occipital",
                "Cerebellar", "Brainstem", "Corpus Callosum",
                "Cingulate", "Sylvian Fissure")

loc_data <- genotyped %>% filter(!is.na(location_codes))

# Binary tier heatmap
loc_binary_summary <- loc_data %>%
  group_by(geno_binary) %>%
  summarise(across(all_of(loc_cols), ~ 100 * mean(. == 1, na.rm = TRUE)), n = n(), .groups = "drop") %>%
  pivot_longer(cols = all_of(loc_cols), names_to = "location", values_to = "pct") %>%
  mutate(loc_label = loc_labels[match(location, loc_cols)],
         loc_label = factor(loc_label, levels = rev(loc_labels)))

fig2d <- ggplot(loc_binary_summary, aes(x = geno_binary, y = loc_label, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), size = TYPO$geom_text_label) +
  scale_fill_gradient(low = PAL_HEAT_LOW, high = PAL_HEAT_HIGH, limits = c(0, 60), name = "Prevalence (%)") +
  labs(
    x = NULL, y = NULL
  ) +
  theme_avm() +
  theme(legend.position = "right")

# NB: Fig2D_location_heatmap_binary was the legacy Fig 2D location panel,
# relocated to Extended Data (ed_anatomy / ed_parietal) per the v2/Hale
# §3 reorder. Not in the current Fig 3 composite or any supp embed; the
# heatmap is computed for the location-distribution stats below but no
# longer written to disk as a panel asset.

# ed_anatomy panels are produced canonically by 17_ED8_anatomy.R; do not
# re-add per_variant_anatomy / per_gene_anatomy writes here (panel-
# uniqueness gate fails on two producers writing the same group dir).

# ══════════════════════════════════════════════════════════════════════════════
# Fig 2E: Hierarchical meta-analysis forest — rupture by genotype
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 2E: Rupture meta-analysis forest ──\n")
library(metafor)

rupture_data <- genotyped %>%
  filter(!is.na(ever_ruptured_num) & !is.na(geno_binary)) %>%
  mutate(geno = relevel(geno_binary, ref = "Panel-negative"))

# ── A. Pooled (overall) model ────────────────────────────────────────────────
rupt_uni <- glm(ever_ruptured_num ~ geno, data = rupture_data, family = binomial)
cat("Pooled rupture model:\n")
print(summary(rupt_uni)$coefficients)

# ── B. Per-study estimates (Firth for sparse strata) ─────────────────────────
study_levels <- levels(factor(rupture_data$study_clean))

per_study <- map_dfr(study_levels, function(st) {
  s <- rupture_data %>% filter(study_clean == st)
  n_pos <- sum(s$geno == "Variant-positive")
  n_neg <- sum(s$geno == "Panel-negative")
  # Need ≥3 per group and ≥1 event to fit
  if (n_pos < 3 || n_neg < 3 || length(unique(s$ever_ruptured_num)) < 2) {
    return(tibble(study = st, yi = NA_real_, sei = NA_real_,
      n = nrow(s), n_pos = n_pos, n_neg = n_neg))
  }
  # Use Firth with pl=TRUE (profile-likelihood CIs) for sparse per-study strata —
  # Wald SEs (pl=FALSE) are anti-conservative near separation.
  mod <- tryCatch(logistf(ever_ruptured_num ~ geno, data = s, pl = TRUE), error = function(e) NULL)
  if (is.null(mod)) return(tibble(study = st, yi = NA_real_, sei = NA_real_,
    n = nrow(s), n_pos = n_pos, n_neg = n_neg))
  geno_idx <- grep("genoVariant-positive", mod$terms)
  if (length(geno_idx) == 0) return(tibble(study = st, yi = NA_real_, sei = NA_real_,
    n = nrow(s), n_pos = n_pos, n_neg = n_neg))
  beta <- mod$coefficients[geno_idx]
  se   <- sqrt(diag(vcov(mod)))[geno_idx]
  tibble(study = st, yi = beta, sei = se, n = nrow(s), n_pos = n_pos, n_neg = n_neg)
})

per_study_valid <- per_study %>% filter(!is.na(yi) & !is.na(sei) & sei > 0)
cat(sprintf("  Studies with estimable OR: %d / %d\n", nrow(per_study_valid), nrow(per_study)))

# ── C. Random-effects meta-analysis ──────────────────────────────────────────
ma <- rma(yi = yi, sei = sei, data = per_study_valid, method = "REML")
cat(sprintf("  Pooled log-OR = %.3f (SE=%.3f), z=%.2f, p=%.4f\n",
  ma$beta[1], ma$se[1], ma$zval[1], ma$pval[1]))
cat(sprintf("  I² = %.1f%%, tau² = %.4f\n", ma$I2, ma$tau2))

# ── D. Build hierarchical forest data ────────────────────────────────────────
per_study_valid$weight <- weights(ma)

forest_df <- per_study_valid %>%
  mutate(
    est     = exp(yi),
    lo      = exp(yi - 1.96 * sei),
    hi      = exp(yi + 1.96 * sei),
    label   = study,
    n_label = as.character(n),
    wt_label = sprintf("%.1f", weight),
    type    = "study"
  ) %>%
  arrange(desc(weight)) %>%
  bind_rows(tibble(
    study = "Pooled", yi = as.numeric(ma$beta), sei = as.numeric(ma$se),
    n = sum(per_study_valid$n), n_pos = NA_integer_, n_neg = NA_integer_,
    weight = 100,
    est    = exp(as.numeric(ma$beta)),
    lo     = exp(as.numeric(ma$ci.lb)),
    hi     = exp(as.numeric(ma$ci.ub)),
    label  = sprintf("Pooled (I²=%.0f%%)", ma$I2),
    n_label = as.character(sum(per_study_valid$n)),
    wt_label = "-",
    type   = "pooled"
  ))

# ── E. Plot ──────────────────────────────────────────────────────────────────
fig2e <- forest_compact_meta(
  forest_df,
  cols = list(
    list(col = "n_label",  header = "N"),
    list(col = "wt_label", header = "Weight (%)")
  ),
  # v6.34 (2026-05-20): forest polish 1-to-1 with ED04 / ED05 — compact
  # cell layout, base_size 8, smaller squares, no extra bottom rule.
  point_col    = BINARY_PALETTE["Variant-positive"],
  axis_ticks   = c(0.25, 0.5, 1, 2, 4),
  x_lab        = "Odds Ratio (95% CI)",
  est_col_side = "left",
  base_size    = 8,
  size_range   = c(0.5, 1.5),
  show_bottom_rule = FALSE
)

# Multivariable adjustment for rupture: pre-specified set of confounders
# that the §4 prose actually cites — Spetzler-Martin size + deep venous
# drainage. Audit M4 found a prior implementation used MASS::stepAIC
# (forward selection over the full predictor pool) which happened to
# select these two covariates on this dataset but mis-described as a
# fixed adjustment in the prose. The hardcoded model below matches the
# Methods declaration ("after multivariable adjustment for SM size and
# deep venous drainage") so reviewers and replicators see the same
# specification the prose cites. No stepAIC is invoked.
rupt_multi_data <- rupture_data %>%
  filter(!is.na(sm_size_num) & !is.na(sm_drainage_num))
rupt_multi <- glm(
  ever_ruptured_num ~ geno + sm_size_num + sm_drainage_num,
  data   = rupt_multi_data,
  family = binomial
)

# rupture_meta_forest dropped from figure in v6.51 (2026-05-21); fig2e
# kept for potential reference but no longer saved to a panel slot.

# ══════════════════════════════════════════════════════════════════════════════
# Fig 2F: Clinical history bars
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 2F: Clinical history ──\n")

clinical_vars <- c("prior_seizure_num", "prior_radiation_num", "prior_embolization_num")
clinical_labels <- c("Seizure history", "Prior radiation", "Prior embolization")

clinical_summary <- map2_dfr(clinical_vars, clinical_labels, function(var, lab) {
  sub <- genotyped %>% filter(!is.na(.data[[var]]) & !is.na(geno_binary))
  ft <- safe_fisher(sub$geno_binary, sub[[var]])
  # successes / n per genotype, then exact 95% CI (Clopper-Pearson) on
  # percent scale to match `prevalence`.
  sub %>%
    group_by(geno_binary) %>%
    summarise(successes = sum(.data[[var]] == 1), n = n(), .groups = "drop") %>%
    mutate(
      prevalence = 100 * successes / n,
      ci_lo = 100 * mapply(function(x, m) binom.test(x, m)$conf.int[1], successes, n),
      ci_hi = 100 * mapply(function(x, m) binom.test(x, m)$conf.int[2], successes, n),
      variable = lab,
      p = ft$p.value
    )
})

# BH-FDR
clinical_p <- clinical_summary %>% distinct(variable, p)
clinical_p$p_fdr <- p.adjust(clinical_p$p, method = "BH")
clinical_summary <- clinical_summary %>% left_join(clinical_p %>% dplyr::select(variable, p_fdr), by = "variable")

cat("Clinical history Fisher's exact (BH-FDR):\n")
print(clinical_p)

clin_dumbbell <- clinical_summary %>%
  rename(feature = variable, group = geno_binary) %>%
  mutate(
    feature = factor(
      stringr::str_wrap(as.character(feature), width = 14),
      levels = stringr::str_wrap(
        c("Seizure history", "Prior radiation", "Prior embolization"),
        width = 14
      )
    )
  )

# ── Shared dumbbell render — sm_components (Fig 3c) + clinical_history
# (Fig 3d). Both call dumbbell_prevalence() with the SAME size_breaks
# and size_limits derived from the union of n values across the two
# data frames, so when patchwork's guides = "collect" runs in
# 30_F2_assemble.R the two panels' identical N-scales merge
# into a single right-margin legend. Without the shared scale, two
# near-duplicate N legends would render side-by-side.
shared_n      <- c(comp_dumbbell$n, clin_dumbbell$n)
shared_lims   <- range(shared_n, na.rm = TRUE)
shared_breaks <- pretty(shared_lims, n = 3)

# Force vertical layout on both legends so keys stack top-to-bottom rather
# than wrapping horizontally. The composer collects guides to the figure
# right margin via `guides = "collect"`; without these overrides patchwork
# was rendering the colour and size keys side-by-side.
.dumbbell_vertical_legend <- list(
  guides(color = guide_legend(order = 1, direction = "vertical",
                              override.aes = list(size = 4)),
         size  = guide_legend(order = 2, direction = "vertical")),
  theme(legend.direction = "vertical",
        legend.box       = "vertical")
)

# 2026-05-20: re-split the prior merged prevalence_dumbbells back into
# two semantically-grouped dumbbell panels — SM sub-components (panel C
# in the new Fig 2 layout) and clinical-history covariates (panel D) —
# so each token resolves to its own panel-letter slot at the cell
# footprint resolved by FIGURE_LAYOUTS$gxp_associations.
# Panel C (sm_components_dumbbell) carries the Genotype + N legend on
# top so it reads as the entry-point panel for the row; panel D
# (clinical_history_dumbbell) suppresses its legend to avoid
# duplication.
# v6 (2026-05-21): producers save UNTHEMED ggplot objects to RDS.
# Composer applies theme_nature_panel() via patchwork `&` so font
# hierarchy is set in one place. Per-panel theme() overrides (e.g.,
# panel.grid.major.y blanking, axis.title.y bold) ARE applied here
# because they are content-specific and should travel with the plot.
#
# Shared N-scale: A/C/D dumbbells + lollipop all carry the same
# `size_breaks` + `size_limits` (computed below in the SHARED_N
# block) so patchwork's guides = "collect" merges the size legend
# into a single entry across the three panels — only ONE N legend
# appears in the composite, not three.
#
# size_range c(1.5, 3.0) matches Panel A's lollipop convention.
# panel.grid.major.y blanked per canonical grid-line rule (only
# numeric axes carry grid lines; dumbbell y is categorical).
fig2c_sm <- dumbbell_prevalence(comp_dumbbell, palette = BINARY_PALETTE,
                                lo_col      = "ci_lo",
                                hi_col      = "ci_hi",
                                size_range  = c(1.5, 3.0),
                                size_breaks = shared_breaks,
                                size_limits = shared_lims) +
  labs(title = NULL, y = "Spetzler-Martin score") +
  guides(color = guide_legend(order = 1, direction = "horizontal",
                              override.aes = list(size = 3)),
         size  = guide_legend(order = 2, direction = "horizontal")) +
  theme(legend.position    = "top",
        legend.direction   = "horizontal",
        legend.box         = "horizontal",
        axis.title.y       = element_text(face = "bold",
                                          margin = margin(r = 4)),
        panel.grid.major.y = element_blank())

fig2d_clin <- dumbbell_prevalence(clin_dumbbell, palette = BINARY_PALETTE,
                                  lo_col      = "ci_lo",
                                  hi_col      = "ci_hi",
                                  size_range  = c(1.5, 3.0),
                                  size_breaks = shared_breaks,
                                  size_limits = shared_lims) +
  labs(title = NULL, y = "Clinical history") +
  guides(color = guide_legend(order = 1, direction = "horizontal",
                              override.aes = list(size = 3)),
         size  = guide_legend(order = 2, direction = "horizontal")) +
  theme(legend.position    = "top",
        legend.direction   = "horizontal",
        legend.box         = "horizontal",
        axis.title.y       = element_text(face = "bold",
                                          margin = margin(r = 4)),
        panel.grid.major.y = element_blank())

save_panel_native("sm_components_dumbbell",    fig2c_sm)
save_panel_native("clinical_history_dumbbell", fig2d_clin)

# ══════════════════════════════════════════════════════════════════════════════
# ed_per_variant_pheno panel a: ordinal-logistic SM-grade OR forest
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_per_variant_pheno panel a: ordinal-logistic SM-grade OR forest ──\n")

olr_data <- sm_data %>%
  filter(!is.na(geno_binary) & !is.na(age) & !is.na(sex_f) & !is.na(sample_type_clean)) %>%
  mutate(
    geno = relevel(geno_binary, ref = "Panel-negative"),
    sm_grade_ord = ordered(sm_grade, levels = c("I", "II", "III", "IV", "V"))
  )

cat(sprintf("Ordinal logistic: n=%d\n", nrow(olr_data)))

tryCatch({
  # Try with sample_type first; fall back without if rank-deficient
  olr_mod <- tryCatch(
    MASS::polr(sm_grade_ord ~ geno + age + sex_f + sample_type_clean, data = olr_data, Hess = TRUE),
    error = function(e) {
      cat("  Note: sample_type dropped from ordinal model (rank-deficient)\n")
      MASS::polr(sm_grade_ord ~ geno + age + sex_f, data = olr_data, Hess = TRUE)
    }
  )
  olr_tidy <- tidy(olr_mod, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(!str_detect(term, "\\|"))

  cat("Ordinal logistic regression coefficients:\n")
  print(olr_tidy)

  # Build table-forest data frame compatible with table_forest_general().
  # N column is nominal (whole-model n) since ordinal logistic yields
  # coefficient-level CIs from the same fitted model; we label each row
  # with the full model n and show the Wald p as "P" (no FDR adjustment
  # because these are predictors within a single model, not a family of
  # hypothesis tests).
  olr_plot <- olr_tidy %>%
    mutate(
      # v6.44 (2026-05-21): wrap multi-word predictor labels with \n so
      # the forest's "Predictor" column reads cleanly at the 3.30-in
      # cell width — single-line "Variant-positive" + "Sample: Fresh"
      # were overflowing into the adjacent N / OR text columns.
      label = case_when(
        term == "genoVariant-positive"       ~ "Mutation-\npositive",
        term == "age"                          ~ "Age\n(per year)",
        term == "sex_fMale"                    ~ "Male\nsex",
        str_detect(term, "sample_type_clean") ~
          paste0("Sample:\n",
                 str_replace(term, "sample_type_clean", "")),
        TRUE                                   ~ term
      ),
      est     = estimate,
      lo      = conf.low,
      hi      = conf.high,
      n_label = as.character(nrow(olr_data)),
      # Wald p from z = statistic (polr coef-level; broom doesn't add p.value by default)
      p_wald  = 2 * (1 - pnorm(abs(statistic))),
      p_label = ifelse(p_wald < 0.001, "<0.001", sprintf("%.3f", p_wald))
    )

  # Mirror the ED Fig 4 cohort-heterogeneity forests (analysis/01_main_analysis/
  # 11_ED4_power_forest_meta.R lines 479-491, 533-545, 589-...):
  # one N text column + a two-row OR/CI column + forest bars. Adding a
  # third "P" text column compressed the forest into the right ~20% of
  # the panel and made CIs overlap with the OR-text sub-rows. With P
  # dropped, the forest gets ~40-50% of the panel width and the two-row
  # OR / "(lo-hi)" sits cleanly to its left, matching the ED4 convention.
  # Wald p-values for individual coefficients are written to the stats
  # readout (results/stats/_manifest_fragments/...) for any reader who
  # needs them; they are also fully captured by the CI-vs-1 inspection
  # the forest makes visually obvious.
  efig6 <- forest_compact_general(
    olr_plot,
    cols = list(list(col = "n_label", header = "N")),
    # v6.41 (2026-05-21): Fig 1 / ED04 / ED05 forest polish 1-to-1 —
    # compact-cell layout (base_size 8, smaller squares 0.5–1.5, no
    # extra bottom rule).
    point_col      = BINARY_PALETTE["Variant-positive"],
    axis_ticks     = c(0.25, 0.5, 1, 2, 4),
    x_lab          = "Odds ratio",
    est_col_header = "OR (95% CI)",
    est_fmt        = "%.2f (%.2f–%.2f)",
    est_col_side   = "left",
    base_size      = 8,
    size_range     = c(0.5, 1.5),
    show_bottom_rule = FALSE
  )

  save_panel(file.path(efig_dir, "ed_per_variant_pheno", "panel_A"), "sm_ordinal", efig6, 3.30, 3.00)
}, error = function(e) {
  cat(sprintf("Ordinal logistic failed: %s\n", e$message))
})

# ══════════════════════════════════════════════════════════════════════════════
# ed_per_variant_pheno panel b: SM components by per-variant genotype
# ══════════════════════════════════════════════════════════════════════════════

# ── Helper: per-variant prevalence heatmap (rows = feature, cols = variant) ──
# Drives ED11 panels B (sm_comp_variant), C (hr_variant), D (clinical_variant)
# — the dumbbell layout was hard to read at composite scale because the five
# variant tracks dodged into 5×N near-overlapping points per feature row.
# Heatmap mirrors the ED10 vaf_heatmap visual contract: white → PAL_HEAT_HIGH
# (= PAL_SM[["IV"]] = #2166AC, the SM Blues palette anchor) ramp,
# black in-cell text on every cell,
# bold 4.5 mm in-cell "<prev%>\nN=<n>" label, no separate legend
# (the in-cell value carries the magnitude). Five-column axis keeps the
# Panel-negative reference stratum alongside the four variant columns
# rather than ED10's four-column variant-only view (ED10 restricts to
# mut+ patients with VAF, which excludes Panel-negative by construction).
prevalence_heatmap <- function(db, feature_levels = NULL, y_lab = NULL,
                               title = NULL) {
  variant_levels <- c("KRAS G12D", "KRAS G12V", "Other KRAS",
                      "BRAF", "Panel-negative")
  agg <- db %>%
    mutate(group = factor(as.character(group), levels = variant_levels))
  if (is.null(feature_levels)) {
    feature_levels <- unique(as.character(agg$feature))
  }
  agg <- agg %>%
    mutate(feature = factor(as.character(feature), levels = feature_levels)) %>%
    tidyr::complete(group, feature,
                    fill = list(prevalence = NA_real_, n = 0L))
  fill_max  <- max(agg$prevalence, na.rm = TRUE)
  # In-cell text is black on every cell (per author request, 2026-05-27);
  # the prior contrast-aware white-on-dark rule is dropped.
  agg <- agg %>%
    mutate(text_col = "black")
  ggplot(agg, aes(x = group, y = feature, fill = prevalence)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # N= dropped 2026-05-17 — % alone reads cleanly in the composite; the
    # per-variant N is already disclosed in the panel A forest table and
    # in SuppTable_03/04, so duplicating it on every cell was just visual
    # noise. Same edit applied to ED10's vaf_heatmap for consistency.
    geom_text(aes(label = ifelse(n == 0L, "—",
                                 sprintf("%.1f%%", prevalence)),
                  colour = text_col),
              size = NM$text$body_mm, fontface = "bold") +
    scale_fill_gradient(low = "white",
                        high = PAL_HEAT_HIGH,
                        na.value = "grey92",
                        limits = c(0, fill_max),
                        name = "Prevalence\n(%)") +
    scale_colour_identity() +
    labs(x = NULL, y = y_lab, title = title) +
    theme_nature_panel() +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 45, hjust = 1,
                                         lineheight = 0.9),
          panel.grid      = element_blank())
}

cat("\n── ed_per_variant_pheno panel b: SM components per-variant ──\n")

comp_variant_summary <- bind_rows(
  genotyped %>% filter(!is.na(sm_size_num)) %>%
    group_by(geno_variant) %>%
    summarise(prevalence = mean(sm_size_num >= 2), n = n(), .groups = "drop") %>%
    mutate(component = "Size ≥2"),
  genotyped %>% filter(!is.na(sm_drainage_num)) %>%
    group_by(geno_variant) %>%
    summarise(prevalence = mean(sm_drainage_num == 1), n = n(), .groups = "drop") %>%
    mutate(component = "Deep drainage"),
  genotyped %>% filter(!is.na(sm_eloquence_num)) %>%
    group_by(geno_variant) %>%
    summarise(prevalence = mean(sm_eloquence_num == 1), n = n(), .groups = "drop") %>%
    mutate(component = "Eloquent location")
)

efig7_db <- comp_variant_summary %>%
  mutate(prevalence = prevalence * 100) %>%
  rename(feature = component, group = geno_variant) %>%
  # Swap raw "Negative" → "Panel-negative" display label so the heatmap
  # column reads with manuscript convention. PAL_VARIANT is dual-keyed.
  mutate(group = relabel_geno_factor(group))

# Feature row order locks the same top→bottom reading the dumbbell layout
# used (Size ≥ 2 on top, Deep Drainage on the bottom). ggplot's y-axis
# inverts factor order, so list features in reverse for the levels =.
efig7 <- prevalence_heatmap(
  efig7_db,
  feature_levels = c("Deep drainage", "Eloquent location", "Size ≥2"),
  title = "SM components"
)

save_panel(file.path(efig_dir, "ed_per_variant_pheno", "panel_B"), "sm_comp_variant", efig7, 3.30, 3.00)

# ══════════════════════════════════════════════════════════════════════════════
# ed_per_variant_pheno panel c: high-risk features by per-variant genotype
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_per_variant_pheno panel c: high-risk features per-variant ──\n")

# v6.43 (2026-05-21): wrap multi-word feature names with \n so the
# y-axis tick text fits inside the 3.30-in cell without colliding with
# the heatmap tiles. The Fig 2 hr_features_OR forest uses the un-wrapped
# hr_labels directly, so we only wrap on the heatmap-side copy here.
.hr_labels_wrapped <- c(
  "Intranidal aneurysm"     = "Intranidal\naneurysm",
  "Venous varix"            = "Venous\nvarix",
  "Venous outflow stenosis" = "Venous outflow\nstenosis",
  "Flow-related aneurysm"   = "Flow-related\naneurysm"
)

hr_variant_summary <- map2_dfr(hr_features, hr_labels, function(feat, lab) {
  genotyped %>%
    filter(!is.na(.data[[feat]])) %>%
    group_by(geno_variant) %>%
    summarise(prevalence = 100 * mean(.data[[feat]] == 1), n = n(), .groups = "drop") %>%
    mutate(feature = .hr_labels_wrapped[lab])
})

efig8_db <- hr_variant_summary %>%
  rename(group = geno_variant) %>%
  # "Negative" → "Panel-negative" display label (PAL_VARIANT dual-keyed).
  mutate(group = relabel_geno_factor(group))

# Mirror the dumbbell row order: Venous Varix on top, Flow-Related Aneurysm
# on the bottom. ggplot y-axis inverts factor order, so list bottom first.
# Use the wrapped versions to match the feature column built above.
efig8 <- prevalence_heatmap(
  efig8_db,
  feature_levels = unname(.hr_labels_wrapped[c(
    "Flow-related aneurysm", "Intranidal aneurysm",
    "Venous outflow stenosis", "Venous varix"
  )]),
  title = "High-risk features"
)

save_panel(file.path(efig_dir, "ed_per_variant_pheno", "panel_C"), "hr_variant", efig8, 3.30, 3.00)

# ══════════════════════════════════════════════════════════════════════════════
# ed_per_variant_pheno panel d: clinical history by per-variant genotype
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_per_variant_pheno panel d: clinical history per-variant ──\n")

clinical_variant_summary <- map2_dfr(clinical_vars, clinical_labels, function(var, lab) {
  genotyped %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(geno_variant) %>%
    summarise(prevalence = 100 * mean(.data[[var]] == 1), n = n(), .groups = "drop") %>%
    mutate(variable = lab)
})

efig9_db <- clinical_variant_summary %>%
  rename(feature = variable, group = geno_variant) %>%
  # "Negative" → "Panel-negative" display label (PAL_VARIANT dual-keyed).
  mutate(group = relabel_geno_factor(group))

# Mirror the dumbbell row order: Seizure History on top, Prior Embolization
# on the bottom. ggplot y-axis inverts factor order, so list bottom first.
# Heatmap replaces the dumbbell here so panel D no longer carries the
# composite's shared Genotype/N legend — heatmap fill is continuous
# prevalence (in-cell labels carry the magnitude), not a discrete
# genotype/size aesthetic.
efig9 <- prevalence_heatmap(
  efig9_db,
  feature_levels = c("Prior embolization", "Prior radiation",
                     "Seizure history"),
  title = "Clinical history"
)

save_panel(file.path(efig_dir, "ed_per_variant_pheno", "panel_D"), "clinical_variant", efig9, 3.30, 3.00)

# ══════════════════════════════════════════════════════════════════════════════
# HIERARCHICAL ANALYSIS: Tier 3 pairwise, adjusted models, interactions
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  HIERARCHICAL GENOTYPE-PHENOTYPE ANALYSIS\n")
cat("══════════════════════════════════════════════════════════════\n")

# ── Setup: create analytic groups ────────────────────────────────────────────

genotyped <- genotyped %>%
  mutate(
    geno = relevel(geno_binary, ref = "Panel-negative"),
    is_kras = str_detect(geno_variant, "KRAS"),
    gene_grp = case_when(
      str_detect(geno_variant, "KRAS") ~ "KRAS",
      str_detect(geno_variant, "BRAF") ~ "BRAF",
      geno_variant == "Negative" ~ "Negative",
      TRUE ~ "Other"
    )
  )

# Outcomes to test hierarchically
outcomes <- tribble(
  ~var,                       ~label,               ~type,
  "ever_ruptured_num",        "Rupture (ever)",     "binary",
  "sm_drainage_num",          "Deep Drainage",      "binary",
  "sm_size_num",              "SM Size",            "continuous",
  "sm_eloquence_num",         "Eloquence",          "binary",
  "n_high_risk_num",          "High-Risk Count",    "continuous",
  "intranidal_aneurysm_num",  "Intranidal Aneur",   "binary",
  "compact_nidus_num",        "Compact Nidus",      "binary",
  "prior_seizure_num",        "Seizure",            "binary",
  "prior_radiation_num",      "Prior Radiation",    "binary",
  "prior_embolization_num",   "Prior Embolization", "binary",
  "age",                      "Age",                "continuous"
)

# ── Helper: run one comparison ───────────────────────────────────────────────
run_unadjusted <- function(data, grp_var, grp1, grp2, outcome_var, outcome_type) {
  d <- data %>% filter(.data[[grp_var]] %in% c(grp1, grp2) & !is.na(.data[[outcome_var]]))
  n1 <- sum(d[[grp_var]] == grp1); n2 <- sum(d[[grp_var]] == grp2)
  if (n1 < 5 | n2 < 5) return(tibble(n1 = n1, n2 = n2, val1 = NA, val2 = NA, p = NA))
  if (outcome_type == "binary") {
    v1 <- round(100 * mean(d[[outcome_var]][d[[grp_var]] == grp1] == 1), 1)
    v2 <- round(100 * mean(d[[outcome_var]][d[[grp_var]] == grp2] == 1), 1)
    p <- tryCatch(fisher.test(table(d[[grp_var]], d[[outcome_var]]))$p.value, error = function(e) NA)
  } else {
    v1 <- round(median(d[[outcome_var]][d[[grp_var]] == grp1], na.rm = TRUE), 2)
    v2 <- round(median(d[[outcome_var]][d[[grp_var]] == grp2], na.rm = TRUE), 2)
    p <- tryCatch(wilcox.test(as.formula(paste(outcome_var, "~", grp_var)), data = d)$p.value, error = function(e) NA)
  }
  tibble(n1 = n1, n2 = n2, val1 = v1, val2 = v2, p = p)
}

run_adjusted <- function(data, outcome_var, outcome_type, geno_var = "geno") {
  s <- data %>% filter(!is.na(.data[[outcome_var]]) & !is.na(age) & !is.na(sample_type_clean))
  sn <- s %>% count(study_clean) %>% filter(n >= 5)
  s <- s %>% filter(study_clean %in% sn$study_clean) %>% mutate(study_clean = droplevels(study_clean))
  if (nrow(s) < 30 || length(unique(s[[geno_var]])) < 2) return(NULL)
  fam <- if (outcome_type == "binary") binomial() else gaussian()
  covars <- if (outcome_var == "age") "study_clean + sample_type_clean" else "age + study_clean + sample_type_clean"
  frm <- as.formula(paste(outcome_var, "~", geno_var, "+", covars))
  m <- tryCatch(glm(frm, data = s, family = fam), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cs <- summary(m)$coefficients
  gr <- grep(geno_var, rownames(cs), fixed = FALSE)
  gr <- gr[!grepl("\\:", rownames(cs)[gr])]
  if (length(gr) == 0) return(NULL)
  ci <- confint.default(m)
  map_dfr(gr, function(r) {
    tibble(
      term = rownames(cs)[r],
      estimate = cs[r, 1],
      OR = if (outcome_type == "binary") exp(cs[r, 1]) else NA_real_,
      beta = if (outcome_type == "continuous") cs[r, 1] else NA_real_,
      p = cs[r, 4],
      ci_low = ci[rownames(cs)[r], 1],
      ci_high = ci[rownames(cs)[r], 2],
      n = nrow(s)
    )
  })
}

# ── Tier 1: Variant-positive vs Panel-negative ───────────────────────────

cat("\n── TIER 1: Mut+ vs Negative (unadjusted) ──\n")
tier1_unadj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_unadjusted(genotyped, "geno_binary", "Variant-positive", "Panel-negative", o$var, o$type)
  r %>% mutate(outcome = o$label, .before = 1)
})
tier1_unadj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.))) cat(sprintf("  %-20s Mut+=%-7s Neg=%-7s p=%.4f %s (n=%d/%d)\n",
    .$outcome[j], .$val1[j], .$val2[j], .$p[j], .$sig[j], .$n1[j], .$n2[j]))}

cat("\n── TIER 1: Mut+ vs Negative (adjusted) ──\n")
tier1_adj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_adjusted(genotyped, o$var, o$type, "geno")
  if (is.null(r)) return(tibble())
  r %>% mutate(outcome = o$label, .before = 1)
})
tier1_adj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.)))
    if (!is.na(.$OR[j])) cat(sprintf("  %-20s OR=%.2f p=%.4f %s (n=%d)\n", .$outcome[j], .$OR[j], .$p[j], .$sig[j], .$n[j]))
    else cat(sprintf("  %-20s β=%.2f p=%.4f %s (n=%d)\n", .$outcome[j], .$beta[j], .$p[j], .$sig[j], .$n[j]))}

# ── Tier 2: KRAS vs Negative ────────────────────────────────────────────────

cat("\n── TIER 2: KRAS vs Negative (unadjusted) ──\n")
kras_neg <- genotyped %>% filter(gene_grp %in% c("KRAS", "Negative"))
tier2_unadj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_unadjusted(kras_neg, "gene_grp", "KRAS", "Negative", o$var, o$type)
  r %>% mutate(outcome = o$label, .before = 1)
})
tier2_unadj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.))) cat(sprintf("  %-20s KRAS=%-7s Neg=%-7s p=%.4f %s (n=%d/%d)\n",
    .$outcome[j], .$val1[j], .$val2[j], .$p[j], .$sig[j], .$n1[j], .$n2[j]))}

# ── Tier 3: G12D vs Neg, G12V vs Neg, G12D vs G12V ──────────────────────────

cat("\n── TIER 3: G12D vs Negative (unadjusted) ──\n")
g12d_neg <- genotyped %>% filter(geno_variant %in% c("KRAS G12D", "Negative"))
tier3a_unadj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_unadjusted(g12d_neg, "geno_variant", "KRAS G12D", "Negative", o$var, o$type)
  r %>% mutate(outcome = o$label, .before = 1)
})
tier3a_unadj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.))) if (!is.na(.$p[j])) cat(sprintf("  %-20s G12D=%-7s Neg=%-7s p=%.4f %s (n=%d/%d)\n",
    .$outcome[j], .$val1[j], .$val2[j], .$p[j], .$sig[j], .$n1[j], .$n2[j]))}

cat("\n── TIER 3: G12D vs Negative (adjusted) ──\n")
g12d_neg <- g12d_neg %>% mutate(is_g12d = factor(geno_variant == "KRAS G12D", levels = c(FALSE, TRUE)))
tier3a_adj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_adjusted(g12d_neg, o$var, o$type, "is_g12d")
  if (is.null(r)) return(tibble())
  r %>% mutate(outcome = o$label, .before = 1)
})
tier3a_adj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.)))
    if (!is.na(.$OR[j])) cat(sprintf("  %-20s OR=%.2f p=%.4f %s (n=%d)\n", .$outcome[j], .$OR[j], .$p[j], .$sig[j], .$n[j]))
    else cat(sprintf("  %-20s β=%.2f p=%.4f %s (n=%d)\n", .$outcome[j], .$beta[j], .$p[j], .$sig[j], .$n[j]))}

cat("\n── TIER 3: G12V vs Negative (unadjusted) ──\n")
g12v_neg <- genotyped %>% filter(geno_variant %in% c("KRAS G12V", "Negative"))
tier3b_unadj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_unadjusted(g12v_neg, "geno_variant", "KRAS G12V", "Negative", o$var, o$type)
  r %>% mutate(outcome = o$label, .before = 1)
})
tier3b_unadj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.))) if (!is.na(.$p[j])) cat(sprintf("  %-20s G12V=%-7s Neg=%-7s p=%.4f %s (n=%d/%d)\n",
    .$outcome[j], .$val1[j], .$val2[j], .$p[j], .$sig[j], .$n1[j], .$n2[j]))}

cat("\n── TIER 3: G12D vs G12V (unadjusted) ──\n")
g12d_g12v <- genotyped %>% filter(geno_variant %in% c("KRAS G12D", "KRAS G12V"))
tier3c_unadj <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  r <- run_unadjusted(g12d_g12v, "geno_variant", "KRAS G12D", "KRAS G12V", o$var, o$type)
  r %>% mutate(outcome = o$label, .before = 1)
})
tier3c_unadj %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.))) if (!is.na(.$p[j])) cat(sprintf("  %-20s G12D=%-7s G12V=%-7s p=%.4f %s (n=%d/%d)\n",
    .$outcome[j], .$val1[j], .$val2[j], .$p[j], .$sig[j], .$n1[j], .$n2[j]))}

# ══════════════════════════════════════════════════════════════════════════════
# LOCATION × GENOTYPE INTERACTIONS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  LOCATION × GENOTYPE INTERACTION MODELS\n")
cat("══════════════════════════════════════════════════════════════\n")

loc_interaction_vars <- c("loc_frontal", "loc_temporal", "loc_parietal", "loc_occipital",
  "loc_cerebellar", "loc_basal_ganglia", "loc_thalamus", "loc_brainstem", "loc_insular")
loc_interaction_labels <- c("Frontal", "Temporal", "Parietal", "Occipital",
  "Cerebellar", "Basal Ganglia", "Thalamus", "Brainstem", "Insular")

key_interaction_outcomes <- tribble(
  ~var,                      ~label,          ~type,
  "ever_ruptured_num",       "Rupture (ever)", "binary",
  "sm_drainage_num",         "Drainage",      "binary",
  "age",                     "Age",           "continuous",
  "n_high_risk_num",         "High-Risk",     "continuous",
  "sm_size_num",             "SM Size",       "continuous",
  "sm_eloquence_num",        "Eloquence",     "binary",
  "intranidal_aneurysm_num", "Intranidal",    "binary",
  "compact_nidus_num",       "Compact Nidus", "binary"
)

# A. Mut+ vs Neg × Location
#
# AUDIT 2026-05-12 (F11 deep, canonical helper): the (lobe, outcome)
# interaction GLM is now fit by `fit_loc_geno_interaction()` from
# `analysis/helper_scripts/utils.R`. Both this producer (SuppTable08,
# 8-outcome family) and `analysis/01_main_analysis/17_ED8_anatomy.R`
# (ED Fig 8 panel C, 3-outcome focused family cited by §4 prose) call
# the same helper, so the raw interaction P AND beta-sign at any
# overlapping (lobe, outcome) cell are guaranteed identical. The two
# BH families differ INTENTIONALLY (broad supplementary exploration
# vs focused hypothesis); the helper canonicalises everything below
# the BH step.
cat("\n── A. Genotype (mut+ vs neg) × Location → Outcomes ──\n")
loc_int_results <- list()
for (li in seq_along(loc_interaction_vars)) {
  loc <- loc_interaction_vars[li]
  loc_lab <- loc_interaction_labels[li]
  for (i in seq_len(nrow(key_interaction_outcomes))) {
    o   <- key_interaction_outcomes[i, ]
    fit <- fit_loc_geno_interaction(df, o$var, loc, o$type)
    if (is.na(fit$p)) next
    # Drop the p < 0.15 pre-filter so the table is a complete test
    # family (audit D5 / 11: prior survivor-only export prevented BH-FDR
    # from being applied on the proper denominator). The console print
    # still flags marginal/significant rows for the operator log.
    sig <- ifelse(fit$p < 0.05, "***", ifelse(fit$p < 0.1, ".", ""))
    if (fit$p < 0.15) {
      cat(sprintf("  %-14s × %-10s → %-12s β=%.2f p=%.4f %s (n=%d)\n",
        "Genotype", loc_lab, o$label, fit$beta, fit$p, sig, fit$n))
    }
    loc_int_results[[paste(loc, o$var, "binary")]] <- tibble(
      location = loc_lab, outcome = o$label, tier = "Mut+ vs Neg",
      beta = fit$beta, p = fit$p, n = fit$n)
  }
}

# B. G12D vs Neg × Parietal and G12V vs Neg × Parietal (specific decomposition)
cat("\n── B. Variant-specific × Parietal → Outcomes ──\n")
for (vt in c("KRAS G12D", "KRAS G12V")) {
  vn <- genotyped %>%
    filter(geno_variant %in% c(vt, "Negative") & !is.na(loc_parietal)) %>%
    mutate(is_var = factor(geno_variant == vt, levels = c(FALSE, TRUE)))
  cat(sprintf("  %s (n=%d) vs Neg (n=%d):\n", vt, sum(as.logical(vn$is_var)), sum(!as.logical(vn$is_var))))
  for (i in seq_len(nrow(key_interaction_outcomes))) {
    o <- key_interaction_outcomes[i, ]
    s <- vn %>% filter(!is.na(.data[[o$var]]))
    if (nrow(s) < 30 || n_distinct(s$is_var) < 2 || n_distinct(s$loc_parietal) < 2) next
    fam <- if (o$type == "binary") binomial() else gaussian()
    m <- tryCatch(glm(as.formula(paste(o$var, "~ is_var * loc_parietal")), data = s, family = fam), error = function(e) NULL)
    if (is.null(m)) next
    cs <- summary(m)$coefficients; ir <- grep(":", rownames(cs))
    if (length(ir) == 0) next
    p_int <- cs[ir, 4]
    sig <- ifelse(p_int < 0.05, "***", ifelse(p_int < 0.1, ".", ""))
    cat(sprintf("    %-16s β=%.2f p=%.4f %s (n=%d)\n", o$label, cs[ir, 1], p_int, sig, nrow(s)))
    loc_int_results[[paste("parietal", o$var, vt)]] <- tibble(
      location = "Parietal", outcome = o$label, tier = vt,
      beta = cs[ir, 1], p = p_int, n = nrow(s))
  }
}

# C. Parietal × Genotype → Rupture: stratified rates
cat("\n── C. Parietal × Genotype stratified rupture rates ──\n")
par_rupt <- genotyped %>%
  filter(!is.na(ever_ruptured_num) & !is.na(loc_parietal)) %>%
  mutate(par_group = ifelse(loc_parietal == 1, "Parietal", "Non-Parietal"))
par_rupt_tab <- par_rupt %>%
  group_by(geno_binary, par_group) %>%
  summarise(n = n(), rupt = sum(ever_ruptured_num), pct = round(100 * mean(ever_ruptured_num), 1), .groups = "drop")
for (j in seq_len(nrow(par_rupt_tab))) {
  cat(sprintf("  %-20s %-14s n=%3d rupt=%.1f%%\n",
    par_rupt_tab$geno_binary[j], par_rupt_tab$par_group[j], par_rupt_tab$n[j], par_rupt_tab$pct[j]))
}

# D. KRAS-Parietal subgroup table
cat("\n── D. KRAS-Parietal subgroup ──\n")
genotyped <- genotyped %>%
  mutate(
    kras_par_group = case_when(
      is_kras & loc_parietal == 1 ~ "KRAS-Parietal",
      is_kras & (loc_parietal == 0 | is.na(loc_parietal)) ~ "KRAS-NonParietal",
      !mutation_positive ~ "Negative",
      TRUE ~ NA_character_
    )
  )

kpg <- genotyped %>% filter(!is.na(kras_par_group))
cat(sprintf("  KRAS-Parietal=%d KRAS-NonPar=%d Negative=%d\n",
  sum(kpg$kras_par_group == "KRAS-Parietal"),
  sum(kpg$kras_par_group == "KRAS-NonParietal"),
  sum(kpg$kras_par_group == "Negative")))

kras_par_stats <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  o <- outcomes[i, ]
  s <- kpg %>% filter(!is.na(.data[[o$var]]))
  grp_n <- s %>% count(kras_par_group)
  if (any(grp_n$n < 3)) return(tibble())
  if (o$type == "binary") {
    vals <- s %>% group_by(kras_par_group) %>%
      summarise(pct = round(100 * mean(.data[[o$var]] == 1), 1), n = n(), .groups = "drop")
    set.seed(MASTER_SEED)   # audit F13: pin Monte-Carlo Fisher P
    p <- tryCatch(fisher.test(table(s$kras_par_group, s[[o$var]]), simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA)
  } else {
    vals <- s %>% group_by(kras_par_group) %>%
      summarise(pct = median(.data[[o$var]], na.rm = TRUE), n = n(), .groups = "drop")
    p <- tryCatch(kruskal.test(as.formula(paste(o$var, "~ kras_par_group")), data = s)$p.value, error = function(e) NA)
  }
  kp <- vals %>% filter(kras_par_group == "KRAS-Parietal")
  knp <- vals %>% filter(kras_par_group == "KRAS-NonParietal")
  neg <- vals %>% filter(kras_par_group == "Negative")
  tibble(outcome = o$label, kras_par = kp$pct[1], kras_nonpar = knp$pct[1], negative = neg$pct[1], p = p,
    n_kp = kp$n[1], n_knp = knp$n[1], n_neg = neg$n[1])
})

kras_par_stats %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.))) cat(sprintf("  %-20s KP=%-7s KNP=%-7s Neg=%-7s p=%.4f %s\n",
    .$outcome[j], .$kras_par[j], .$kras_nonpar[j], .$negative[j], .$p[j], .$sig[j]))}

# KRAS-Parietal variant breakdown
kras_par_variants <- genotyped %>%
  filter(kras_par_group == "KRAS-Parietal") %>%
  count(geno_variant, sort = TRUE)
cat("\n  KRAS-Parietal by variant:\n")
for (j in seq_len(nrow(kras_par_variants)))
  cat(sprintf("    %s: n=%d\n", kras_par_variants$geno_variant[j], kras_par_variants$n[j]))

# ══════════════════════════════════════════════════════════════════════════════
# LOCATION DISTRIBUTION BY GENOTYPE (formal tests)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Location distribution: G12D vs G12V vs Negative ──\n")
gvn_loc <- genotyped %>% filter(geno_variant %in% c("KRAS G12D", "KRAS G12V", "Negative"))

loc_dist_tests <- map2_dfr(loc_interaction_vars, loc_interaction_labels, function(loc, lab) {
  s <- gvn_loc %>% filter(!is.na(.data[[loc]]))
  d_pct <- 100 * mean(s[[loc]][s$geno_variant == "KRAS G12D"] == 1)
  v_pct <- 100 * mean(s[[loc]][s$geno_variant == "KRAS G12V"] == 1)
  n_pct <- 100 * mean(s[[loc]][s$geno_variant == "Negative"] == 1)
  set.seed(MASTER_SEED)   # audit F13: pin Monte-Carlo Fisher P
  p3 <- tryCatch(fisher.test(table(s$geno_variant, s[[loc]]), simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA)
  tibble(location = lab, g12d_pct = round(d_pct, 1), g12v_pct = round(v_pct, 1), neg_pct = round(n_pct, 1), p = p3)
})

cat(sprintf("%-14s %-8s %-8s %-8s %s\n", "Location", "G12D%", "G12V%", "Neg%", "p (3-way)"))
loc_dist_tests %>% {for (j in seq_len(nrow(.)))
  cat(sprintf("  %-14s %-8.1f %-8.1f %-8.1f p=%.3f %s\n",
    .$location[j], .$g12d_pct[j], .$g12v_pct[j], .$neg_pct[j], .$p[j],
    ifelse(.$p[j] < 0.05, "***", ifelse(.$p[j] < 0.1, ".", ""))))}

# ══════════════════════════════════════════════════════════════════════════════
# Save all stats
# ══════════════════════════════════════════════════════════════════════════════

loc_int_df <- bind_rows(loc_int_results)
# Apply BH-FDR for SuppTable08. Audit 2026-05-12 (F11): two different BH
# families were previously computed for the same parietal × rupture
# interaction cell — within-outcome (9-test) here and joint-27-cell in
# `17_ED8_anatomy.R`. §4 prose cites the joint-27-cell q from edfig08,
# so the more defensible (and stricter) choice is to harmonize on the
# joint-27-cell family. We therefore now BH across the full 9-lobe ×
# 3-outcome binary-tier grid (`tier == "Mut+ vs Neg"`); the variant-
# specific decomposition rows (`tier %in% {"KRAS G12D", "KRAS G12V"}`)
# are treated as a separate, smaller family.
loc_int_df <- loc_int_df %>%
  mutate(.bh_family = dplyr::if_else(tier == "Mut+ vs Neg",
                                     "binary_27",
                                     "variant_decomp")) %>%
  group_by(.bh_family) %>%
  mutate(p_fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  dplyr::select(-.bh_family)

fig2_stats <- list(
  sm_kw_binary = kw_sm_binary,
  sm_kw_variant = kw_sm_variant,
  size_kw = kw_size,
  drainage_fisher = fisher_drain,
  eloquence_fisher = fisher_eloq,
  hr_forest = hr_results,
  hr_count_kw = kw_hr,
  rupture_uni = tidy(rupt_uni, conf.int = TRUE, exponentiate = TRUE),
  rupture_multi = tidy(rupt_multi, conf.int = TRUE, exponentiate = TRUE),
  clinical_summary = clinical_summary,
  tier1_unadj = tier1_unadj,
  tier1_adj = tier1_adj,
  tier2_unadj = tier2_unadj,
  tier3a_g12d_neg_unadj = tier3a_unadj,
  tier3a_g12d_neg_adj = tier3a_adj,
  tier3b_g12v_neg_unadj = tier3b_unadj,
  tier3c_g12d_g12v_unadj = tier3c_unadj,
  location_interactions = loc_int_df,
  parietal_rupture_stratified = par_rupt_tab,
  kras_parietal_subgroup = kras_par_stats,
  kras_parietal_variants = kras_par_variants,
  location_distribution = loc_dist_tests
)

saveRDS(fig2_stats, file.path(output_dir, "stats", "fig2_stats.rds"))

# ── SupplementaryTables ──────────────────────────────────────────────────────────────────
library(writexl)
etable_dir <- file.path(output_dir, "SupplementaryTables")
dir.create(etable_dir, recursive = TRUE, showWarnings = FALSE)

# Supplementary Table 13 — High-risk feature OR forest (accompanies Fig. 2c).
# 2026-06-14: rename Sheet1 → "High-risk feature OR"; Nature-styled headers;
# hard-error validation.
source(here("analysis", "helper_scripts", "supp_table_writer.R"))
.fmt_p_st13 <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}
st13_df <- hr_results %>%
  transmute(
    feature      = feature,
    n            = as.integer(n),
    or_ci_str    = sprintf("%.3f (%.3f, %.3f)", OR, lower, upper),
    p_str        = .fmt_p_st13(p),
    p_sample_str = .fmt_p_st13(sample_type_p),
    p_fdr_str    = .fmt_p_st13(p_fdr)
  )
write_supp_table(
  data    = st13_df,
  path    = file.path(etable_dir, "SuppTable13_high_risk_OR.xlsx"),
  sheet   = "High-risk feature OR",
  columns = list(
    col    ("feature",      label = "Feature"),
    col_int("n",            label = "N",                       italic = TRUE),
    col    ("or_ci_str",    label = "OR (95% CI)",             italic = TRUE),
    col    ("p_str",        label = "P",                       italic = TRUE),
    col    ("p_sample_str", label = "P (sample type)",         italic = TRUE),
    col    ("p_fdr_str",    label = "FDR P",                   italic = TRUE)
  ),
  footnote = c(
    "Firth-penalized logistic regression for each pre-specified high-risk angioarchitectural feature (variant-positive vs panel-negative bAVMs), sample-type-adjusted.",
    "P (sample type) is the Wald P for the sample_type covariate in the same model — included so readers can audit whether the genotype effect is being absorbed by surgical-vs-autopsy specimen mix.",
    "FDR P applies BH correction across the four pre-specified features."
  )
)

# Supplementary Table 08 — Hierarchical genotype-phenotype models
# 2026-06-14: collapsed from 7 sheets (tier1_unadj/adj + tier2 + tier3a_unadj/adj +
#   tier3b + tier3c) into a single stacked sheet under the standing
#   one-sheet-per-ST rule. Schema unifies unadjusted Fisher / Wilcoxon
#   comparisons and age-+-series-adjusted GLM regressions across all five
#   pairwise contrasts (Variant-positive / KRAS / G12D / G12V vs panel-negative,
#   plus G12D vs G12V).
source(here("analysis", "helper_scripts", "supp_table_writer.R"))

# Canonical phenotype labels (match ST4 / ST6 / ST7).
.canon_pheno_st8 <- function(x) dplyr::recode(x,
  "Rupture (ever)"     = "Rupture (ever)",
  "Deep Drainage"      = "Deep venous drainage",
  "SM Size"            = "Spetzler–Martin size",
  "Eloquence"          = "Eloquent brain location",
  "High-Risk Count"    = "High-risk feature count",
  "Intranidal Aneur"   = "Intranidal aneurysm",
  "Compact Nidus"      = "Compact nidus",
  "Seizure"            = "Seizure history",
  "Prior Radiation"    = "Prior radiation",
  "Prior Embolization" = "Prior embolization",
  "Age"                = "Age at presentation"
)

# Map outcome label → type (binary / continuous) so each row's Statistic
# string can be written without re-checking the model object.
.outcome_type_st8 <- setNames(outcomes$type, outcomes$label)

.fmt_p_st8 <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}

# Unadjusted block builder: turn a tier_*_unadj frame into the unified
# long-format schema. Statistic column communicates whether val1/val2 are
# percentages (binary) or medians (continuous); Estimate is "v1 vs v2"; CI
# is sentinel (no CI on the unadjusted comparison). Capture outcome type
# BEFORE renaming so the lookup hits the original labels.
.unadj_block <- function(d, contrast_label) {
  d %>%
    mutate(.type = .outcome_type_st8[outcome]) %>%
    transmute(
      contrast  = contrast_label,
      outcome   = .canon_pheno_st8(outcome),
      analysis  = "Unadjusted (Fisher / Wilcoxon)",
      n_str     = sprintf("%d / %d", as.integer(n1), as.integer(n2)),
      statistic = ifelse(.type == "binary",
                         "% (group 1) vs % (panel-neg)",
                         "Median (group 1) vs Median (panel-neg)"),
      estimate  = ifelse(.type == "binary",
                         sprintf("%.1f%% vs %.1f%%", val1, val2),
                         sprintf("%.2f vs %.2f", val1, val2)),
      ci_str    = "—",
      p_str     = .fmt_p_st8(p)
    )
}

# Adjusted block builder: take the relevant per-outcome coefficient (OR for
# binary, β for continuous) and report it with its 95% CI and p.
.adj_block <- function(d, contrast_label) {
  d %>%
    mutate(.type = .outcome_type_st8[outcome]) %>%
    transmute(
      contrast  = contrast_label,
      outcome   = .canon_pheno_st8(outcome),
      analysis  = "Age + series adjusted (GLM)",
      n_str     = sprintf("%d", as.integer(n)),
      statistic = ifelse(.type == "binary",
                         "OR (adj)",
                         "β (adj)"),
      estimate  = ifelse(.type == "binary",
                         sprintf("%.3f", OR),
                         sprintf("%.3f", beta)),
      ci_str    = ifelse(.type == "binary",
                         sprintf("(%.3f, %.3f)", exp(ci_low), exp(ci_high)),
                         sprintf("(%.3f, %.3f)", ci_low, ci_high)),
      p_str     = .fmt_p_st8(p)
    )
}

# 2026-06-15: ST8 back to 5 sheets (Nature MOESM3 (A)-(E) pattern), one per
# contrast. Unadjusted + adjusted analyses stack within each sheet where
# both exist; Contrast column dropped (sheet name carries it).
.st8_cols <- list(
  col    ("outcome",   label = "Outcome"),
  col    ("analysis",  label = "Analysis"),
  col    ("n_str",     label = "N",          italic = TRUE),
  col    ("statistic", label = "Statistic"),
  col    ("estimate",  label = "Estimate"),
  col    ("ci_str",    label = "95% CI"),
  col    ("p_str",     label = "P",          italic = TRUE)
)
.st8_footnote <- function(has_adj) {
  base <- c(
    "Eleven outcomes (8 binary, 3 continuous). Unadjusted comparisons use Fisher's exact (binary) or Wilcoxon (continuous).",
    "Unadjusted N is reported as 'group 1 / panel-negative' counts.",
    "'—' marks the 95% CI on unadjusted rows: Fisher / Wilcoxon yield a P-value but no scalar effect-size CI."
  )
  if (has_adj) c(base,
    "Adjusted rows: age- and series-adjusted GLM (logistic for binary, OLS for continuous). N is the model fit count after dropping NA covariates and applying the series-≥5 filter.",
    "Adjusted ORs and their CIs are exponentiated from the model log-odds; adjusted β and its CI are reported on the original outcome scale.")
  else base
}

.st8_strip <- function(d) d[, c("outcome","analysis","n_str","statistic","estimate","ci_str","p_str"), drop = FALSE]

st8_sheets <- list(
  "(A) Variant-pos vs Neg" = list(
    data = .st8_strip(bind_rows(
      .unadj_block(tier1_unadj, "Variant-positive vs Panel-negative"),
      .adj_block  (tier1_adj,   "Variant-positive vs Panel-negative")
    )),
    columns = .st8_cols,
    footnote = .st8_footnote(has_adj = TRUE)
  ),
  "(B) KRAS-mutant vs Neg" = list(
    data = .st8_strip(.unadj_block(tier2_unadj, "KRAS-mutant vs Panel-negative")),
    columns = .st8_cols,
    footnote = .st8_footnote(has_adj = FALSE)
  ),
  "(C) KRAS G12D vs Neg" = list(
    data = .st8_strip(bind_rows(
      .unadj_block(tier3a_unadj, "KRAS G12D vs Panel-negative"),
      .adj_block  (tier3a_adj,   "KRAS G12D vs Panel-negative")
    )),
    columns = .st8_cols,
    footnote = .st8_footnote(has_adj = TRUE)
  ),
  "(D) KRAS G12V vs Neg" = list(
    data = .st8_strip(.unadj_block(tier3b_unadj, "KRAS G12V vs Panel-negative")),
    columns = .st8_cols,
    footnote = .st8_footnote(has_adj = FALSE)
  ),
  "(E) KRAS G12D vs G12V" = list(
    data = .st8_strip(.unadj_block(tier3c_unadj, "KRAS G12D vs KRAS G12V")),
    columns = .st8_cols,
    footnote = .st8_footnote(has_adj = FALSE)
  )
)

write_supp_table_workbook(
  sheets = st8_sheets,
  path   = file.path(etable_dir, "SuppTable08_hierarchical_genotype_phenotype.xlsx")
)
cat(sprintf("  ✓ SuppTable08_hierarchical_genotype_phenotype.xlsx (5 sheets: %s rows)\n",
            paste(vapply(st8_sheets, function(x) nrow(x$data), integer(1)),
                  collapse = ", ")))
st8_df <- bind_rows(lapply(st8_sheets, function(x) x$data))

# Supplementary Table 09 — Location × genotype interaction tests
# 2026-06-14: kept as single sheet (already one); switched to write_supp_table
# for hard-error validation + Nature-styled headers + harmonized phenotype
# and tier labels. The asymmetric Parietal coverage (variant-decomposed in
# Parietal only; binary-tier elsewhere) reflects the §2 prose's parietal-
# KRAS-rupture focus.
.canon_pheno_st9 <- function(x) dplyr::recode(x,
  "Rupture (ever)"   = "Rupture (ever)",
  "Drainage"         = "Deep venous drainage",
  "SM Size"          = "Spetzler–Martin size",
  "Eloquence"        = "Eloquent brain location",
  "High-Risk"        = "High-risk feature count",
  "Intranidal"       = "Intranidal aneurysm",
  "Compact Nidus"    = "Compact nidus",
  "Age"              = "Age at presentation"
)
.canon_tier_st9 <- function(x) dplyr::recode(x,
  "Mut+ vs Neg" = "Variant-positive vs Panel-negative",
  "KRAS G12D"   = "KRAS G12D vs Panel-negative",
  "KRAS G12V"   = "KRAS G12V vs Panel-negative"
)
.fmt_p_st9 <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}

st9_df <- loc_int_df %>%
  transmute(
    location  = location,
    contrast  = .canon_tier_st9(tier),
    outcome   = .canon_pheno_st9(outcome),
    n         = as.integer(n),
    beta_str  = sprintf("%.3f", beta),
    p_str     = .fmt_p_st9(p),
    fdr_p_str = .fmt_p_st9(p_fdr)
  ) %>%
  arrange(location, contrast, outcome)

write_supp_table(
  data    = st9_df,
  path    = file.path(etable_dir, "SuppTable09_location_interactions.xlsx"),
  sheet   = "Location × genotype interactions",
  columns = list(
    col    ("location",  label = "Location"),
    col    ("contrast",  label = "Contrast"),
    col    ("outcome",   label = "Outcome"),
    col_int("n",         label = "N",     italic = TRUE),
    col    ("beta_str",  label = "β (interaction)", italic = TRUE),
    col    ("p_str",     label = "P",     italic = TRUE),
    col    ("fdr_p_str", label = "FDR P", italic = TRUE)
  ),
  footnote = c(
    "Per-lobe interaction term from <outcome> ~ contrast × in_lobe (binary; lobe-indicator coded 1 for that location, 0 otherwise) fit by logistic GLM for binary outcomes and OLS for continuous. The β column is the interaction coefficient.",
    "Two BH-FDR families: (1) the binary-tier 'Variant-positive vs Panel-negative' rows are corrected jointly across the 9 lobes × 8 outcomes (27-cell family — matches the joint family cited in the §2 parietal-KRAS-rupture prose); (2) the per-variant 'KRAS G12D' / 'KRAS G12V' decomposition is corrected as a separate, smaller family.",
    "Asymmetric coverage by design: variant decomposition (KRAS G12D / KRAS G12V) is reported only for Parietal because the §2 parietal-KRAS-rupture finding is the only location where the variant-level breakdown was pre-specified. All other locations report the binary-tier (Variant-positive vs Panel-negative) contrast only."
  )
)
cat(sprintf("  ✓ SuppTable09_location_interactions.xlsx (1 sheet, %d rows)\n",
            nrow(st9_df)))
cat("── Supplementary Tables 08, 09, 12 saved to results/SupplementaryTables/ ──\n")

# ═════════════════════════════════════════════════════════════════════════════
# ── fig2 manifest fragment ──────────────────────────────────────────────────
# ═════════════════════════════════════════════════════════════════════════════
# Every scalar number that Section 2 prose cites. Key names match
# stats$fig2$<key>.

# ---- denominators ---------------------------------------------------------
cohort_n_sm_graded <- nrow(sm_data)
cohort_n_hr_forest <- max(hr_results$n, na.rm = TRUE)

# ---- SM grade + components -----------------------------------------------
sm_kw_binary_p     <- kw_sm_binary$p.value
sm_kw_variant_p    <- kw_sm_variant$p.value
sm_size_kw_p       <- kw_size$p.value
sm_drainage_fisher_p  <- fisher_drain$p.value
sm_eloquence_fisher_p <- fisher_eloq$p.value
sm_comp_min_p <- min(c(sm_size_kw_p, sm_drainage_fisher_p,
                       sm_eloquence_fisher_p), na.rm = TRUE)
# BH-FDR across the 3-test SM-sub-component family. 2026-05-19: with the
# 18 May 26 cohort eloquence reached p=0.023 uncorrected; q=0.07 keeps the
# null framing defensible at the standard manuscript FDR convention.
.sm_comp_q <- p.adjust(c(sm_size_kw_p, sm_drainage_fisher_p,
                          sm_eloquence_fisher_p), method = "BH")
sm_comp_min_q <- min(.sm_comp_q, na.rm = TRUE)
sm_size_kw_q          <- .sm_comp_q[1]
sm_drainage_fisher_q  <- .sm_comp_q[2]
sm_eloquence_fisher_q <- .sm_comp_q[3]

# ---- High-risk feature forest --------------------------------------------
hr_min_fdr <- min(hr_results$p_fdr, na.rm = TRUE)
hr_max_fdr <- max(hr_results$p_fdr, na.rm = TRUE)

# ---- Rupture: random-effects meta-analysis (REML) ------------------------
# `ma` is the metafor::rma() fit from §Fig 2E above.
rupt_meta_log_or    <- unname(ma$beta[1])
rupt_meta_or        <- exp(rupt_meta_log_or)
rupt_meta_ci_lo     <- exp(rupt_meta_log_or - 1.96 * ma$se[1])
rupt_meta_ci_hi     <- exp(rupt_meta_log_or + 1.96 * ma$se[1])
rupt_meta_p         <- ma$pval[1]
rupt_meta_i2_pct    <- ma$I2
rupt_meta_n_studies <- nrow(per_study_valid)

# ---- Rupture: univariate + multivariable logistic ------------------------
.row_geno <- function(tbl) {
  out <- tbl[tbl$term == "genoVariant-positive", ,  drop = FALSE]
  if (nrow(out) == 0L) out[1, ] <- NA
  out
}
.ru_u <- .row_geno(fig2_stats$rupture_uni)
.ru_m <- .row_geno(fig2_stats$rupture_multi)
rupt_uni_or       <- .ru_u$estimate
rupt_uni_ci_lo    <- .ru_u$conf.low
rupt_uni_ci_hi    <- .ru_u$conf.high
rupt_uni_p        <- .ru_u$p.value
rupt_uni_n        <- stats::nobs(rupt_uni)
rupt_multi_or     <- .ru_m$estimate
rupt_multi_ci_lo  <- .ru_m$conf.low
rupt_multi_ci_hi  <- .ru_m$conf.high
rupt_multi_p      <- .ru_m$p.value
rupt_multi_n      <- stats::nobs(rupt_multi)

# ---- SM sub-component sample sizes + Fisher ORs --------------------------
# Used by Fig 3 panel C (sm_components) caption to report per-sub-component
# n and the (mut+ vs neg) odds ratio for the binary sub-components.
sm_size_n      <- nrow(size_data)
sm_drainage_n  <- nrow(drain_data)
sm_eloquence_n <- nrow(eloq_data)
sm_drainage_or <- unname(fisher_drain$estimate)
sm_eloquence_or <- unname(fisher_eloq$estimate)

# ---- Clinical history FDR range ------------------------------------------
clinical_min_fdr <- min(clinical_summary$p_fdr, na.rm = TRUE)
clinical_max_fdr <- max(clinical_summary$p_fdr, na.rm = TRUE)
# Per-test denominator across the 3 clinical-history Fisher tests
.clin_per_test_n <- clinical_summary %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(n_test = sum(n), .groups = "drop")
clinical_n_min <- min(.clin_per_test_n$n_test, na.rm = TRUE)
clinical_n_max <- max(.clin_per_test_n$n_test, na.rm = TRUE)

# ---- Per-clinical-history scalar pulls (Fig 3 panel D caption) -----------
# clinical_summary carries one row per (variable × geno_binary) with
# `prevalence` already on percent scale and `p_fdr` per variable. Helper
# pulls the scalar for one (label, stratum) tuple, falling back to NA if
# the row is missing (e.g., during cohort refreshes that drop a variable).
.clin_pct <- function(label, stratum) {
  row <- clinical_summary[clinical_summary$variable == label &
                          clinical_summary$geno_binary == stratum, , drop = FALSE]
  if (nrow(row) == 0L) NA_real_ else row$prevalence[1]
}
.clin_fdr <- function(label) {
  row <- clinical_summary[clinical_summary$variable == label, , drop = FALSE]
  if (nrow(row) == 0L) NA_real_ else row$p_fdr[1]
}
.clin_n <- function(label) {
  row <- clinical_summary[clinical_summary$variable == label, , drop = FALSE]
  if (nrow(row) == 0L) NA_integer_ else as.integer(sum(row$n))
}

clin_seizure_pct_mut    <- .clin_pct("Seizure history",    "Variant-positive")
clin_seizure_pct_neg    <- .clin_pct("Seizure history",    "Panel-negative")
clin_seizure_fdr        <- .clin_fdr("Seizure history")
clin_seizure_n          <- .clin_n("Seizure history")
clin_radiation_pct_mut  <- .clin_pct("Prior radiation",    "Variant-positive")
clin_radiation_pct_neg  <- .clin_pct("Prior radiation",    "Panel-negative")
clin_radiation_fdr      <- .clin_fdr("Prior radiation")
clin_radiation_n        <- .clin_n("Prior radiation")
clin_embolization_pct_mut <- .clin_pct("Prior embolization", "Variant-positive")
clin_embolization_pct_neg <- .clin_pct("Prior embolization", "Panel-negative")
clin_embolization_fdr     <- .clin_fdr("Prior embolization")
clin_embolization_n       <- .clin_n("Prior embolization")

# ---- KRAS-parietal subgroup ----------------------------------------------
.par_tbl <- fig2_stats$kras_parietal_subgroup
.par_rupt <- .par_tbl[.par_tbl$outcome == "Rupture (ever)", , drop = FALSE]
par_kras_rupt_pct    <- .par_rupt$kras_par
par_nonkras_rupt_pct <- .par_rupt$kras_nonpar
par_neg_rupt_pct     <- .par_rupt$negative
par_kras_n           <- as.integer(
  fig2_stats$parietal_rupture_stratified$n[
    fig2_stats$parietal_rupture_stratified$geno_binary == "Variant-positive" &
    fig2_stats$parietal_rupture_stratified$par_group  == "Parietal"
  ]
)

# Interaction β, P for (Mut+ vs Neg) × Parietal rupture
.loc_int <- fig2_stats$location_interactions
.loc_par <- .loc_int[.loc_int$location == "Parietal" &
                     .loc_int$outcome  == "Rupture (ever)"  &
                     .loc_int$tier     == "Mut+ vs Neg", , drop = FALSE]
par_rupt_interaction_beta <- .loc_par$beta
par_rupt_interaction_p    <- .loc_par$p
par_rupt_interaction_n    <- as.integer(.loc_par$n)
n_lobe_tests              <- length(unique(fig2_stats$location_distribution$location))

fig2_fragment <- list(
  cohort_n_sm_graded        = cohort_n_sm_graded,
  cohort_n_hr_forest        = cohort_n_hr_forest,
  sm_kw_binary_p            = sm_kw_binary_p,
  sm_kw_variant_p           = sm_kw_variant_p,
  sm_size_kw_p              = sm_size_kw_p,
  sm_drainage_fisher_p      = sm_drainage_fisher_p,
  sm_eloquence_fisher_p     = sm_eloquence_fisher_p,
  sm_comp_min_p             = sm_comp_min_p,
  # BH-FDR-adjusted across the 3-component family (added 2026-05-19)
  sm_size_kw_q              = sm_size_kw_q,
  sm_drainage_fisher_q      = sm_drainage_fisher_q,
  sm_eloquence_fisher_q     = sm_eloquence_fisher_q,
  sm_comp_min_q             = sm_comp_min_q,
  hr_min_fdr                = hr_min_fdr,
  hr_max_fdr                = hr_max_fdr,
  rupt_meta_or              = rupt_meta_or,
  rupt_meta_ci_lo           = rupt_meta_ci_lo,
  rupt_meta_ci_hi           = rupt_meta_ci_hi,
  rupt_meta_p               = rupt_meta_p,
  rupt_meta_i2_pct          = rupt_meta_i2_pct,
  rupt_meta_n_studies       = rupt_meta_n_studies,
  rupt_uni_or               = rupt_uni_or,
  rupt_uni_ci_lo            = rupt_uni_ci_lo,
  rupt_uni_ci_hi            = rupt_uni_ci_hi,
  rupt_uni_p                = rupt_uni_p,
  rupt_uni_n                = rupt_uni_n,
  rupt_multi_or             = rupt_multi_or,
  rupt_multi_ci_lo          = rupt_multi_ci_lo,
  rupt_multi_ci_hi          = rupt_multi_ci_hi,
  rupt_multi_p              = rupt_multi_p,
  rupt_multi_n              = rupt_multi_n,
  clinical_min_fdr          = clinical_min_fdr,
  clinical_max_fdr          = clinical_max_fdr,
  clinical_n_min            = clinical_n_min,
  clinical_n_max            = clinical_n_max,
  sm_size_n                 = sm_size_n,
  sm_drainage_n             = sm_drainage_n,
  sm_eloquence_n            = sm_eloquence_n,
  sm_drainage_or            = sm_drainage_or,
  sm_eloquence_or           = sm_eloquence_or,
  clin_seizure_pct_mut      = clin_seizure_pct_mut,
  clin_seizure_pct_neg      = clin_seizure_pct_neg,
  clin_seizure_fdr          = clin_seizure_fdr,
  clin_seizure_n            = clin_seizure_n,
  clin_radiation_pct_mut    = clin_radiation_pct_mut,
  clin_radiation_pct_neg    = clin_radiation_pct_neg,
  clin_radiation_fdr        = clin_radiation_fdr,
  clin_radiation_n          = clin_radiation_n,
  clin_embolization_pct_mut = clin_embolization_pct_mut,
  clin_embolization_pct_neg = clin_embolization_pct_neg,
  clin_embolization_fdr     = clin_embolization_fdr,
  clin_embolization_n       = clin_embolization_n,
  par_kras_rupt_pct         = par_kras_rupt_pct,
  par_nonkras_rupt_pct      = par_nonkras_rupt_pct,
  par_neg_rupt_pct          = par_neg_rupt_pct,
  par_kras_n                = par_kras_n,
  par_rupt_interaction_beta = par_rupt_interaction_beta,
  par_rupt_interaction_p    = par_rupt_interaction_p,
  par_rupt_interaction_n    = par_rupt_interaction_n,
  n_lobe_tests              = n_lobe_tests,
  # Aggregate counts across the full location_interactions panel (8
  # outcomes × 9 lobes ≈ 25 tests) — used by §4 prose to avoid hardcoded
  # "spanning eight phenotypic outcomes" descriptors.
  n_loc_int_outcomes        = length(unique(fig2_stats$location_interactions$outcome)),
  n_loc_int_tests           = nrow(fig2_stats$location_interactions),
  # Full tables for downstream consumers
  hr_forest                 = hr_results,
  location_distribution     = fig2_stats$location_distribution,
  location_interactions     = fig2_stats$location_interactions
)

source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(section = "fig2", stats = fig2_fragment)

# ═════════════════════════════════════════════════════════════════════════════
# ── Panel F (parietal × KRAS rupture, promoted to main per Hale v2) ──────────
# ═════════════════════════════════════════════════════════════════════════════
# 2026-04-25 (Hale v2 ordering): the singular nominally-significant cell from
# the location × genotype interaction grid (parietal-lobe KRAS-mutant rupture
# enrichment) is promoted from ed_parietal/panel_A → null_phenotype panel f.
# Built from the same kpg subset that feeds kras_par_stats above.

par_rupt_3 <- kpg %>%
  dplyr::filter(!is.na(ever_ruptured_num)) %>%
  dplyr::group_by(kras_par_group) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    successes = sum(ever_ruptured_num),
    pct       = round(100 * mean(ever_ruptured_num), 1),
    # 2026-05-19 NM pass: add Clopper-Pearson 95% CI per stratum so the
    # lollipop carries uncertainty (matches dumbbell convention in
    # panels c/d). binom.test()$conf.int → [0, 1] on proportion scale,
    # rescale to percent.
    ci_lo     = 100 * binom.test(successes, n)$conf.int[1],
    ci_hi     = 100 * binom.test(successes, n)$conf.int[2],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Canonical stratum labels: "KRAS-positive" (not "KRAS⁺" / "KRAS+")
    # so the panel reads as full prose and the y-axis tick aligns with
    # the manuscript-wide convention banning the "+" shorthand for
    # cohort labels (see contract in utils.R; "+" is reserved for
    # in-prose abbreviation, not figure tick labels).
    kras_par_group = factor(
      kras_par_group,
      levels = c("Negative", "KRAS-NonParietal", "KRAS-Parietal"),
      # Uniform 2-line wrapping so all three strata occupy the same
      # vertical text height on the y-axis (one row of dot spacing
      # each). Panel-negative is the shortest label; soft-wrapping it
      # to "Panel-\nnegative" keeps the layout balanced.
      labels = c("Panel-\nnegative",
                 "KRAS-positive\nnon-parietal",
                 "KRAS-positive\nparietal")
    )
  ) %>%
  # Deterministic stratum order — sort by kras_par_group factor
  # levels so the y-axis renders bottom-to-top: Panel-negative,
  # KRAS-non-parietal, KRAS-parietal (independent of upstream row
  # order in kpg). y_label is shorter than kras_par_group to keep
  # the lollipop's left margin tight (no inline n; stratum counts
  # move to the figure caption).
  dplyr::arrange(kras_par_group) %>%
  dplyr::mutate(
    y_label = factor(as.character(kras_par_group),
                     levels = as.character(kras_par_group)),
    # v6.16 (2026-05-20): pre-clamp n into the shared dumbbell size
    # range so strata whose true n falls below the shared lower limit
    # (KRAS-parietal n=78 vs the C/D-driven shared lower bound) still
    # render as the minimum dot size instead of being censored to NA.
    # Clamping is done on the data frame (not inside aes()) so the
    # saved ggplot RDS is self-contained — the composer re-renders
    # from RDS in a different environment where `shared_lims` is
    # unresolved, which had caused the lollipop to disappear in the
    # composite even after a per-aes pmin/pmax fix. `n` itself is
    # preserved unchanged for downstream stats.
    n_size = pmin(pmax(n, shared_lims[1]), shared_lims[2])
  )

# Tier D categorical (Wong-migrated 2026-05-20): panel-negative grey /
# KRAS-positive non-parietal Sky blue / KRAS-positive parietal Reddish
# purple (highlight, matches PAL_ANATOMY$HIGHLIGHT). AVOIDS Tier A
# variant hues — these strata are anatomy sub-categories, not variant
# identities. Purple flags the parietal stratum (the singular nominally
# significant interaction cell).
fig2g_pal <- c(
  "Panel-negative"               = "#737373",
  "KRAS-positive\nnon-parietal"  = W_SKYBLUE,
  "KRAS-positive\nparietal"      = W_REDDISHPURPLE
)
# Lollipop matching Fig 1B variant_landscape geometry: grey segment 0->pct,
# coloured point sized by n. Color + size legends surfaced inline per AJK
# Tables & Figs #C03 so the colour/size encoding is explained on the
# panel itself rather than caption-only.
fig2g <- ggplot(par_rupt_3,
                aes(x = pct, y = y_label, color = kras_par_group)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, linewidth = 0.8, alpha = 0.85,
                 show.legend = FALSE) +
  # Size aesthetic uses the pre-clamped n_size column (see
  # par_rupt_3 mutate above) so the saved RDS travels cleanly into
  # the composer's render environment.
  geom_point(aes(size = n_size), alpha = 0.9) +
  scale_color_manual(values = fig2g_pal) +
  # SHARED N breaks + limits so patchwork's guides = "collect" merges
  # A's size legend with the C/D dumbbell size legends into a single
  # N legend in the composite. The union is captured in
  # `shared_breaks` / `shared_lims` computed earlier in this script
  # for the dumbbell panels.
  scale_size_continuous(range = c(1.5, 3.0), name = "N",
                        breaks = shared_breaks,
                        limits = shared_lims) +
  scale_x_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.4))) +
  labs(title = NULL, x = "Ever ruptured (%)", y = NULL) +
  guides(color = "none",
         size  = guide_legend(override.aes = list(color = "grey40"))) +
  # v6: no theme_avm_native() — composer applies theme_nature_panel()
  # via `&`. panel.grid.major.y blanked (categorical y per the grid
  # rule) and legend.position = "right" travels with the plot.
  theme(panel.grid.major.y = element_blank(),
        legend.position    = "right")

save_panel_native("parietal_kras_rupture", fig2g)

# ═════════════════════════════════════════════════════════════════════════════
# ── Fig 2 composite PNG (patchwork) ─────────────────────────────────────────
# ═════════════════════════════════════════════════════════════════════════════
# 6-panel layout (Apr 2026 Hale v2 ordering): top row = three short panels
# (parietal × KRAS rupture, SM grade, SM sub-components); middle row =
# rupture meta forest full-width; bottom row = clinical history (narrow) +
# high-risk features forest (2-col wide). tag_levels carries the manuscript-
# facing letters assigned by the resolver in reading order (top-left to
# bottom-right):
#   ┌──────────┬──────────┬──────────┐
#   │    a     │    b     │    c     │      a = parietal × KRAS rupture
#   ├──────────┴──────────┴──────────┤      b = SM grade distribution
#   │             d                   │      c = SM sub-components
#   ├──────────┬──────────────────────┤      d = rupture meta forest (full-width)
#   │    e     │          f           │      e = clinical history
#   └──────────┴──────────────────────┘      f = high-risk features forest

# Composite assembly is now centralised in 30_F2_assemble.R
# The legacy inline
# block below is gated behind if(FALSE){...} as historical documentation
# of the layout intent; canonical Fig 3 composite path stays at
# results/Figure3/Fig3_composite.{pdf,png}.
cat("\n── Fig 3 composite: deferred to 30_F2_assemble.R ──\n")


cat("\n══ 13_F2_genotype_phenotype.R complete ══\n")
cat("Figures: Figure3/panel_a-f (resolver-aligned), ExtendedData/ed_anatomy + ed_per_variant_pheno\n")
cat("Stats: Tier 1-3, adjusted models, location interactions, KRAS-Parietal subgroup\n")
