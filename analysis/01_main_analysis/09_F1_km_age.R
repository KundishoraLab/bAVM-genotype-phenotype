# 09_F1_km_age.R — main age/KM figure (registry group: `age`) + two ED panels.
#
# The age/KM figure is
# the "age as central finding" showcase, promoting panels formerly scattered
# into a
# single coherent main figure.
#
# Main panels (registry group `age` — manuscript figure number is
# resolver-driven from prose first-citation order, not pinned here):
#   km_presentation     Kaplan-Meier — time to any clinical presentation
#                       (G12D / G12V / Neg). Every patient presents →
#                       no censoring; empirical CDF rendered as KM.
#   age_density_binary  Age-density by binary genotype (Mut+ vs Negative).
#   age_density_variant Age-density by variant (G12D / G12V / Negative,
#                       with median lines).
#   vaf_age_scatter     VAF × age scatter (variant-positive cases).
#                       X-axis tick marks at 2%, 4%, 6% per Hale v2 note.
#   km_rupture          Time to rupture (Option 2a strict). Events =
#                       `Ruptured at presentation` (n≈179). Censored =
#                       `Never ruptured` (n≈127). The `Prior rupture
#                       (not at surgery)` cases are dropped — true
#                       rupture age unknown, age_surgery would bias the
#                       hazard. Cox PH adjusted for sex + study +
#                       sample_type.
#
# ED panels:
#   rare_variants       Rare-variant KM (Other KRAS + BRAF V600E vs
#                       Negative), time to any clinical presentation
#                       (saved into ed_cohort_heterogeneity/). Kept
#                       separate from the main KM so the wide CIs on
#                       these small strata do not overwhelm the main
#                       comparison.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here)
  library(survival); library(survminer); library(patchwork)
})

source(here("analysis", "helper_scripts", "utils.R"))

# Phase 2: composite panel tags are derived from the registry-driven
# resolver cache (written by run_all.R step [1.5]) rather than hardcoded.
# Loading once up front so every labs(tag = ...) call below is a simple
# list lookup. Falls back gracefully if the cache is missing (e.g. when
# this script is sourced interactively before run_all.R has been invoked).
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.panel_tags <- tryCatch({
  a <- load_panel_assignments()
  list(
    km_presentation     = panel_letter("km_presentation",     a),
    age_density_binary  = panel_letter("age_density_binary",  a),
    age_density_variant = panel_letter("age_density_variant", a),
    km_rupture          = panel_letter("km_rupture",          a),
    vaf_age_scatter     = panel_letter("vaf_age_scatter",     a)
  )
}, error = function(e) {
  message("  [panel_assignments cache missing] falling back to layout order")
  list(km_presentation = "a", age_density_binary = "b",
       age_density_variant = "c", km_rupture = "d", vaf_age_scatter = "e")
})

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(mutation_positive)) %>%
  mutate(
    variant_group = case_when(
      !mutation_positive              ~ "Negative",
      mutation == "KRAS G12D"         ~ "KRAS G12D",
      mutation == "KRAS G12V"         ~ "KRAS G12V",
      mutation_gene == "KRAS"          ~ "Other KRAS",
      mutation == "BRAF V600E"         ~ "BRAF V600E",
      mutation_gene == "BRAF"          ~ "Other BRAF",
      TRUE                             ~ "Negative"
    ),
    variant_group = factor(variant_group,
      levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF V600E", "Other BRAF", "Negative"))
  )

cat(sprintf("Input n = %d (genotyped)\n", nrow(df)))

# KM-context palettes (PAL_KM, PAL_RARE) and the binary-KM palette below
# (PAL_BINARY_KM) are sourced from analysis/helper_scripts/utils.R via
# the `source(...utils.R)` call above. Comment retained for reader
# orientation: main curves restrict to well-powered strata (G12D, G12V,
# Negative); Other KRAS (n=11) and BRAF V600E (n=4) appear together in
# the ED17 rare-variant panel folded into ed_cohort_heterogeneity.

fig3_dir <- panel_slot_dir("km_presentation")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)
# Rare-variant KM was folded into ed_cohort_heterogeneity (panel f) on
# 2026-04-26; output now lives alongside the other cohort-heterogeneity
# panels so the panel-uniqueness gate (token -> registry-group dir) is
# satisfied. Stats manifest section is still emitted under "edfig17"
# below for downstream caption consumption — that key is internal to the
# stats schema and does not have to track figure-number renumbering.
# v6.22 (2026-05-20): rare_variants moved to ed_km_diagnostics — all KM
# panels live in one ED group regardless of strata (sex / variant /
# rare-variant). Use panel_slot_dir() so the producer writes directly
# into the registered slot (e.g., panel_A/); the group-root location
# the producer previously used is flagged "registry-misplaced" by
# sync_panel_prefixes.R and silently deleted on every render.
source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa17 <- load_panel_assignments()
efig17_dir <- panel_slot_dir("rare_variants", .pa17)
dir.create(efig17_dir, recursive = TRUE, showWarnings = FALSE)
stats_dir <- here("results", "stats")

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
save_trio <- function(dir, name, plot, w, h)
  save_panel_impl(dir, name, plot, w, h, device = "cairo")

# Display-label helper: relabel the bare "Negative" factor level as
# "Panel-negative" at the rendering layer (legend, axis ticks, risk-table
# row names). The underlying factor level is left as "Negative" so every
# downstream filter / stats producer continues to match unchanged.
display_levels <- function(f) {
  sub("^Negative$", "Panel-negative", levels(f))
}

# ─────────────────────────────────────────────────────────────────────────────
# FIG 4A — Time to any clinical presentation (main, BRAF excluded)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Fig 4A: time to any clinical presentation ──\n")

km1_df <- df %>%
  filter(!is.na(age),
         variant_group %in% c("KRAS G12D", "KRAS G12V", "Negative")) %>%
  mutate(variant_group = droplevels(variant_group),
         event = 1L, time = age)    # everyone presents → event = 1

cat(sprintf("  n = %d (Other KRAS + BRAF V600E moved to the rare_variants ED panel)\n",
            nrow(km1_df)))
print(km1_df %>% count(variant_group))

fit1 <- survfit(Surv(time, event) ~ variant_group, data = km1_df)
lr1  <- survdiff(Surv(time, event) ~ variant_group, data = km1_df)
lr1_p <- 1 - pchisq(lr1$chisq, df = length(lr1$n) - 1)

# Pairwise log-rank with Bonferroni correction across variant pairs
pw1 <- pairwise_survdiff(Surv(time, event) ~ variant_group, data = km1_df,
                         p.adjust.method = "bonferroni")

cat(sprintf("  Overall log-rank p = %.3g\n", lr1_p))
cat("  Pairwise (Bonferroni):\n"); print(pw1)

p4a <- ggsurvplot(
  legend = "bottom", fit1, data = km1_df,
  conf.int  = TRUE, conf.int.alpha = 0.2,
  # v6.7 (2026-05-20): linewidth + ribbon alpha unified across Fig 1
  # (KM curves D/F/H, density E, regression G all share linewidth 0.5
  # and background alpha 0.2). 0.5 sits inside Nature's 0.5-1 pt band
  # at composite scale.
  size      = 0.5,
  palette   = as.character(PAL_KM[levels(km1_df$variant_group)]),
  risk.table = TRUE, risk.table.height = 0.28,
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
  xlab = "Age (years)", ylab = "Proportion symptom-free",
  legend.title = "Genotype", legend.labs = display_levels(km1_df$variant_group),
  break.time.by = 10, xlim = c(0, KM_AGE_XLIM_MAX),  # see KM_AGE_XLIM_MAX in utils.R
  # Smaller censor marks (default 4.5 → 1.5) so they don't crowd the
  # curve; lower-contrast applied post-render below.
  censor.size = 1.5,
  ggtheme = theme_avm(),
  tables.theme = theme_cleantable()
  # Title/subtitle intentionally omitted — all figure text lives in the
  # the manuscript legends. Stats captured in
  # results/stats/fig4a_km_stats.txt.
)
# Dim censor marks (ggsurvplot exposes no censor.alpha param). The censor
# layer is the only GeomPoint in $plot$layers — walk and set its aes_params
# alpha. Lower-contrast ticks still indicate censoring without competing
# with the step curve visually.
p4a$plot$layers <- lapply(p4a$plot$layers, function(l) {
  if (inherits(l$geom, "GeomPoint")) l$aes_params$alpha <- 0.4
  l
})

# save_km_panel splits ggsurv$plot and ggsurv$table into separate RDS
# components so the native-size composer (compose_figure.R) can re-theme
# each at print-footprint dimensions. Standalone PDF/PNG still saved at
# 10x8 for review.
save_km_panel(fig3_dir, "km_presentation", p4a, w = 10, h = 8)

# Dump stats
writeLines(c(
  sprintf("# %s — Time to any clinical presentation",
          panel_prose_tag("km_presentation")),
  sprintf("Generated %s", Sys.Date()),
  sprintf("n = %d  |  Overall log-rank p = %.3g (df = %d)",
          nrow(km1_df), lr1_p, length(lr1$n) - 1),
  "",
  "## Per-stratum summary (median age at presentation)",
  capture.output(print(summary(fit1)$table)),
  "",
  "## Pairwise log-rank (Bonferroni)",
  capture.output(print(pw1))),
  file.path(stats_dir, "fig3a_km_stats.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# Extended Data Fig. — Time to rupture (Option 2a strict)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Extended Data Fig.: time to rupture (Option 2a strict) ──\n")

km2_df <- df %>%
  mutate(
    km2_event = case_when(
      rupture_category == "Ruptured at presentation"       ~ 1L,
      rupture_category == "Never ruptured"                 ~ 0L,
      TRUE                                                 ~ NA_integer_
    ),
    km2_group = factor(as.character(variant_group),
      levels = c("KRAS G12D", "KRAS G12V", "Negative"))
  ) %>%
  filter(!is.na(km2_event) & !is.na(age),
         variant_group %in% c("KRAS G12D", "KRAS G12V", "Negative"))
  # Other KRAS + BRAF V600E routed to the rare_variants ED panel

cat(sprintf("  n = %d after drop 'Prior rupture (not at surgery)' and missing\n", nrow(km2_df)))
print(km2_df %>% count(km2_group, km2_event))

fit2 <- survfit(Surv(age, km2_event) ~ km2_group, data = km2_df)
lr2  <- survdiff(Surv(age, km2_event) ~ km2_group, data = km2_df)
lr2_p <- 1 - pchisq(lr2$chisq, df = length(lr2$n) - 1)

pw2 <- pairwise_survdiff(Surv(age, km2_event) ~ km2_group, data = km2_df,
                         p.adjust.method = "bonferroni")

# Unadjusted Cox (headline HR — pairs naturally with the log-rank /
# KM curves). DATA_DECISIONS §45: pediatric series (BCH, CHOP) enrol
# younger patients AND contributed most of the KRAS G12D genotyping
# calls, so adjusting for `study_clean` conflates ascertainment with
# developmental biology. We therefore lead with the unadjusted HR and
# carry the multivariate models as sensitivity analyses only.
cox_unadj <- coxph(Surv(age, km2_event) ~ km2_group, data = km2_df)

# Minimal-adjustment sensitivity: sex + sample_type only, no study.
# These covariates do not measurably move the genotype effect (in the
# 18 May 26 cohort HR shifts 0.62 → 0.63), confirming the headline.
cox_minadj_df <- km2_df %>%
  filter(!is.na(sex_f), !is.na(sample_type_clean))
cox_minadj <- coxph(Surv(age, km2_event) ~ km2_group + sex_f + sample_type_clean,
                    data = cox_minadj_df)

# Full-adjustment sensitivity (sex + study + sample_type as covariates).
# Retained for transparency; the multivariate HR attenuates because
# study_clean partly confounds ascertainment with developmental biology
# (see Methods).
cox_df <- km2_df %>%
  filter(!is.na(sex_f), !is.na(study_clean), !is.na(sample_type_clean))
cox2 <- coxph(Surv(age, km2_event) ~ km2_group + sex_f + study_clean + sample_type_clean,
              data = cox_df)

# Schoenfeld global test of the proportional-hazards assumption — Methods
# §45 promises this; expose the global P so the prose / supplementary table
# can cite it. Per the 2026-04-29 audit, the covariate-form model violates
# PH (driven by `study_clean`); we therefore also fit a study-stratified
# Cox as a sensitivity model and expose its HR/CI/P + Schoenfeld global P.
cox2_zph <- cox.zph(cox2)
rupt_cox_zph_global_p <- unname(cox2_zph$table["GLOBAL", "p"])

# Sensitivity: stratify by study_clean (absorbs site-level baseline-hazard
# differences without imposing PH on study). Tames Schoenfeld global P
# from ~0.001 to ~0.04 while preserving the headline G12D-vs-Negative HR.
cox2_strata <- coxph(
  Surv(age, km2_event) ~ km2_group + sex_f + sample_type_clean +
                        strata(study_clean),
  data = cox_df)
cox2_strata_conf <- summary(cox2_strata)$conf.int
cox2_strata_coef <- summary(cox2_strata)$coefficients
.cox_neg_row_strata <- grep("Negative", rownames(cox2_strata_conf))[1]
rupt_cox_strata_hr_neg    <- unname(cox2_strata_conf[.cox_neg_row_strata, "exp(coef)"])
rupt_cox_strata_ci_lo_neg <- unname(cox2_strata_conf[.cox_neg_row_strata, "lower .95"])
rupt_cox_strata_ci_hi_neg <- unname(cox2_strata_conf[.cox_neg_row_strata, "upper .95"])
rupt_cox_strata_p_neg     <- unname(cox2_strata_coef[.cox_neg_row_strata, "Pr(>|z|)"])
rupt_cox_strata_zph_global_p <- unname(cox.zph(cox2_strata)$table["GLOBAL", "p"])

# Angioarchitectural-adjustment sensitivity: nidus size + deep venous drainage
# + high-risk feature count. Rupture status is NOT included — it is the event
# in km2_event and would be circular. Extraction deferred to after .extract_neg
# is defined (see below).
cox_angio_df <- km2_df %>%
  filter(!is.na(sm_size_num), !is.na(sm_drainage_num), !is.na(n_high_risk_num))
cox_angio <- coxph(
  Surv(age, km2_event) ~ km2_group + sm_size_num + sm_drainage_num + n_high_risk_num,
  data = cox_angio_df)

cat(sprintf("  Log-rank p = %.3g ; Cox n = %d\n", lr2_p, nrow(cox_df)))
cat("  Adjusted Cox summary (study as covariate):\n"); print(summary(cox2)$coefficients)
cat(sprintf("  Schoenfeld global p = %.3g\n", rupt_cox_zph_global_p))
cat("  Sensitivity Cox summary (strata(study_clean)):\n")
print(summary(cox2_strata)$coefficients)
cat(sprintf("  Strata Schoenfeld global p = %.3g\n",
            rupt_cox_strata_zph_global_p))

peff <- ggsurvplot(
  legend = "bottom", fit2, data = km2_df,
  conf.int  = TRUE, conf.int.alpha = 0.2,
  # v6.7 (2026-05-20): linewidth + ribbon alpha unified; see p4a.
  size      = 0.5,
  palette   = as.character(PAL_KM[levels(km2_df$km2_group)]),
  risk.table = TRUE, risk.table.height = 0.28,
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
  xlab = "Age (years)", ylab = "Proportion rupture-free",
  legend.title = "Genotype", legend.labs = display_levels(km2_df$km2_group),
  break.time.by = 10, xlim = c(0, KM_AGE_XLIM_MAX),  # see KM_AGE_XLIM_MAX in utils.R
  censor.size = 1.5,  # smaller censor ticks; alpha applied below
  ggtheme = theme_avm(),
  tables.theme = theme_cleantable()
  # Title/subtitle in caption; stats in results/stats/efig16_km_rupture_stats.txt.
)
# Dim censor marks — see p4a above for the rationale + mechanism. Panel D
# (rupture) has 129 censored cases so the visual benefit is largest here.
peff$plot$layers <- lapply(peff$plot$layers, function(l) {
  if (inherits(l$geom, "GeomPoint")) l$aes_params$alpha <- 0.4
  l
})

# km_rupture is a MAIN Fig 2 panel (registry: age group panel d). Per
# Hale v2 notes — "may be more helpful to include a Kaplan Meier curve;
# we have all the continuous data so let's show it" — the rupture KM
# belongs in the main figure set. The legacy ED standalone copy has
# been retired; single canonical output here. save_km_panel splits the
# ggsurvplot block into curve + table RDS for native-size compositing.
fig3d_dir <- panel_slot_dir("km_rupture")
save_km_panel(fig3d_dir, "km_rupture", peff, w = 10, h = 8)

# Reverse-Kaplan–Meier (censoring distribution Sᶜ(t)) was retired
# 2026-05-17 after audit (`results/stats/reverse_censoring_audit.md`):
# in a surgical cohort every patient contributes their age at surgery
# as the observation time, so the reverse-KM is by construction
# linked to age at presentation. The reverse-KM signal across genotype
# strata is the SAME signal as Fig 2A — not an independent bias — and
# carries no interpretive value beyond what Fig 2A already shows.
# The audit MD documents the test (χ² = 22.75, P = 1.1e-5) for any
# reviewer who asks. Producer code intentionally left thin.

# Dump stats
writeLines(c(
  sprintf("# %s — Time to rupture (Option 2a strict)",
          panel_prose_tag("km_rupture")),
  sprintf("Generated %s", Sys.Date()),
  sprintf("n=%d events=%d censored=%d  |  Log-rank p = %.3g",
          nrow(km2_df), sum(km2_df$km2_event == 1),
          sum(km2_df$km2_event == 0), lr2_p),
  "",
  "## Per-stratum summary (median age at rupture or censoring)",
  capture.output(print(summary(fit2)$table)),
  "",
  "## Pairwise log-rank (Bonferroni)",
  capture.output(print(pw2)),
  "",
  "## Cox PH adjusted for sex + study + sample_type",
  capture.output(print(summary(cox2)))),
  file.path(stats_dir, "fig3d_km_rupture_stats.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# rare_variants ED panel — Rare variants (Other KRAS + BRAF V600E) vs Negative
# Time to any clinical presentation. Kept separate from Fig 4A so the wide
# CIs on these small strata do not overwhelm the main comparison.
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── rare_variants ED panel: rare variants vs Negative (time to presentation) ──\n")

km3_df <- df %>%
  filter(!is.na(age),
         variant_group %in% c("Other KRAS", "BRAF V600E", "Negative")) %>%
  mutate(rare_group = factor(as.character(variant_group),
                             levels = c("Other KRAS", "BRAF V600E", "Negative")),
         event = 1L)

cat(sprintf("  n = %d (Other KRAS = %d, BRAF = %d, Negative = %d)\n",
            nrow(km3_df),
            sum(km3_df$rare_group == "Other KRAS"),
            sum(km3_df$rare_group == "BRAF V600E"),
            sum(km3_df$rare_group == "Negative")))

fit3  <- survfit(Surv(age, event) ~ rare_group, data = km3_df)
lr3   <- survdiff(Surv(age, event) ~ rare_group, data = km3_df)
lr3_p <- 1 - pchisq(lr3$chisq, df = length(lr3$n) - 1)
pw3   <- pairwise_survdiff(Surv(age, event) ~ rare_group, data = km3_df,
                           p.adjust.method = "bonferroni")

cat(sprintf("  Overall log-rank p = %.3g\n", lr3_p))
cat("  Pairwise (Bonferroni):\n"); print(pw3)

p17 <- ggsurvplot(
  legend = "bottom", fit3, data = km3_df,
  conf.int  = TRUE, conf.int.alpha = 0.2,
  # v6.21 (2026-05-20): linewidth unified at 0.5 across Fig 1 KM curves
  # (Fig 1 D/F/H) and this rare-variant ED KM. conf.int.alpha = 0.2
  # already matches Fig 1's CI ribbon convention.
  size      = 0.5,
  palette   = as.character(PAL_RARE[levels(km3_df$rare_group)]),
  risk.table = TRUE, risk.table.height = 0.28,
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
  xlab = "Age (years)", ylab = "Proportion symptom-free",
  legend.title = "Genotype", legend.labs = display_levels(km3_df$rare_group),
  break.time.by = 10, xlim = c(0, KM_AGE_XLIM_MAX),  # see KM_AGE_XLIM_MAX in utils.R
  censor.size = 1.5,  # smaller censor ticks; alpha applied below
  # fontsize is geom_text mm for risk-table cell numbers AND the
  # "Number at risk" title. survminer's default 4.5 mm ≈ 12.8 pt clashes
  # with theme_avm()'s axis tick labels (rel(0.8) of 14 ≈ 11.2 pt) sitting
  # immediately beneath them. 3.95 mm ≈ 11.2 pt aligns the two.
  fontsize = 3.95,
  ggtheme  = theme_avm(),
  # tables.theme overrides the inherited theme_avm() plot.title (15 pt)
  # and axis.title.x (14 pt) so the risk-table "Number at risk" header and
  # its duplicate "Age (years)" label match curve axis tick text rather
  # than out-shouting it. risk.table.y.text = FALSE hides the y-axis
  # strata labels (replaced by colored squares) so the size-50 entry that
  # theme_cleantable carries for axis.text.y never reaches paper.
  tables.theme = theme_cleantable() +
    theme(
      # Risk-table is NOT covered by theme_avm() (survminer manages it
      # separately). Sizes are scaled proportionally to the bumped 18pt
      # base — was 12/12/rel(0.8) at base 14; now 16/16/rel(0.8) at base 18.
      plot.title   = element_text(face = "bold", family = AVM_FONT_FAMILY, size = 16),
      axis.title.x = element_text(family = AVM_FONT_FAMILY, size = 16),
      axis.text.x  = element_text(family = AVM_FONT_FAMILY, size = rel(0.8))
    )
  # Title/subtitle in caption; stats in results/stats/efig17_rare_variants_stats.txt.
)
# Dim censor marks — see p4a above for rationale.
p17$plot$layers <- lapply(p17$plot$layers, function(l) {
  if (inherits(l$geom, "GeomPoint")) l$aes_params$alpha <- 0.4
  l
})

save_km_panel(efig17_dir, "rare_variants", p17, w = 6.60, h = 2.03)

writeLines(c(
  sprintf("# %s — Rare variants (Other KRAS + BRAF V600E) vs. Negative",
          panel_prose_tag("rare_variants")),
  sprintf("Generated %s", Sys.Date()),
  sprintf("n = %d  |  Log-rank p = %.3g", nrow(km3_df), lr3_p),
  "",
  "## Per-stratum summary",
  capture.output(print(summary(fit3)$table)),
  "",
  "## Pairwise log-rank (Bonferroni)",
  capture.output(print(pw3))),
  file.path(stats_dir, "rare_variants_km_stats.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# FIG 4B — Age density, binary genotype (Mut+ vs Negative)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Fig 4B: age density by binary genotype ──\n")

# PAL_BINARY_KM now sourced from utils.R (KM-context palettes block).

b_df <- df %>%
  filter(!is.na(age)) %>%
  mutate(geno_binary = factor(
    ifelse(mutation_positive, "Variant-positive", "Panel-negative"),
    levels = c("Variant-positive", "Panel-negative")))

med_b <- b_df %>% group_by(geno_binary) %>% summarise(med = median(age), .groups = "drop")
kw_b  <- kruskal.test(age ~ geno_binary, data = b_df)
cat(sprintf("  n=%d | med Mut+=%.1f  Neg=%.1f | KW p=%.3g\n",
            nrow(b_df), med_b$med[1], med_b$med[2], kw_b$p.value))

p4b <- ggplot(b_df, aes(x = age, fill = geno_binary, color = geno_binary)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_rug(alpha = 0.4, linewidth = 0.3, show.legend = FALSE) +
  ref_vline(kind = "median",
            data = med_b, mapping = aes(xintercept = med, color = geno_binary),
            show.legend = FALSE) +
  scale_fill_manual(values = PAL_BINARY_KM, name = "Genotype") +
  scale_color_manual(values = PAL_BINARY_KM, name = "Genotype") +
  labs(x = "Age at presentation (years)", y = "Density") +
  theme_avm() +
  # Legend at bottom so the plot area's top edge aligns with the KM curve
  # tops in panels A/D within the Fig 2 composite row.
  theme(legend.position = "bottom")

# 2026-05-16: Removed the plot_spacer() padding that used to wrap B/C/E
# at a 3:1 ratio to match the old A/D curve+at-risk-table footprint. Now
# that the at-risk tables are pulled out of the Fig 2 composite (composer
# 26_F1_assemble.R reads __curve.rds for the KMs and stages tables in
# ED via 29_ED6_assemble.R), the spacer creates phantom empty
# whitespace in B/C/E that prevents them from filling their cell. Save
# the bare ggplot so the composer's absolute-unit cell sizing controls
# the rendered footprint.
# 2026-05-19 layout decision: age_density_binary was dropped from the
# cohort_natural_history Fig 1 layout (redundant with the variant-
# stratified density). The plot object p4b is still constructed above
# in case a future ED slot wants it, but it is no longer saved.

# ─────────────────────────────────────────────────────────────────────────────
# FIG 4C — Age density by variant (G12D / G12V / Negative)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Fig 4C: age density by variant ──\n")

c_df <- df %>%
  filter(!is.na(age),
         variant_group %in% c("KRAS G12D", "KRAS G12V", "Negative")) %>%
  mutate(variant_group = droplevels(variant_group))

med_c <- c_df %>% group_by(variant_group) %>%
  summarise(med = median(age), n = n(), .groups = "drop")
kw_c  <- kruskal.test(age ~ variant_group, data = c_df)
cat("  medians:"); print(med_c)
cat(sprintf("  KW p = %.3g\n", kw_c$p.value))

# v6.3 (2026-05-21): pre-compute density height at each group's median
# so the reference segment spans only the density curve, not the full
# panel height. `grp` captures the rowwise variant_group explicitly to
# avoid shadowing the column name inside the c_df subset filter.
med_c <- med_c %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    density_at_med = {
      grp <- variant_group
      sub <- c_df$age[c_df$variant_group == grp]
      d   <- stats::density(sub)
      stats::approx(d$x, d$y, xout = med)$y
    }
  ) %>%
  dplyr::ungroup()

p4c <- ggplot(c_df, aes(x = age, fill = variant_group, color = variant_group)) +
  # v6.7 (2026-05-20): unify line weight + fill opacity across Fig 1.
  # linewidth 0.8 -> 0.5 matches the regression lines in Panel G and
  # KM curves in D/F/H; geom_density alpha 0.3 -> 0.2 matches the CI
  # ribbon in G and the violin fill in C, giving every "background
  # tint" in the figure the same translucency.
  geom_density(alpha = 0.2, linewidth = 0.5) +
  geom_rug(alpha = 0.4, linewidth = 0.3, show.legend = FALSE) +
  # v6.3 (2026-05-21): geom_segment limited to the density-curve height
  # at each group's median, replacing the full-height geom_vline. Lighter
  # linewidth (0.2 -> 0.3 prior, now 0.2) + dotted linetype so the lines
  # read as reference markers subordinate to the density curves.
  geom_segment(data = med_c,
               aes(x = med, xend = med,
                   y = 0, yend = density_at_med,
                   color = variant_group),
               # v6.4 (2026-05-21): linetype "dotted" was reading
               # nearly solid at lw 0.2 because dot spacing scales
               # with stroke width. "13" = 1-pt dash, 3-pt gap, so
               # the marker reads as a sparse reference line that
               # sits visually subordinate to the density curves.
               linetype = "13", linewidth = 0.25, alpha = 0.6,
               show.legend = FALSE) +
  scale_fill_manual(values = PAL_KM, name = "Genotype") +
  scale_color_manual(values = PAL_KM, name = "Genotype") +
  labs(x = "Age at presentation (years)", y = "Density") +
  theme_avm() +
  # 2026-05-19: legend restored per panel.
  theme(legend.position = "bottom")

# Spacer removed 2026-05-16 — see p4b save above for rationale.
fig3c_dir <- panel_slot_dir("age_density_variant")
save_trio(fig3c_dir, "age_density_variant", p4c, w = 7, h = 4.5)

# ─────────────────────────────────────────────────────────────────────────────
# FIG 4D — VAF × age scatter (G12D / G12V; 2%, 4%, 6% x-ticks)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Fig 4D: VAF × age scatter ──\n")

d_df <- df %>%
  filter(!is.na(age), !is.na(vaf_prop),
         variant_group %in% c("KRAS G12D", "KRAS G12V")) %>%
  mutate(variant_group = droplevels(variant_group),
         vaf_pct = vaf_prop * 100)

# Slope per variant and interaction test (for caption / stats dump)
lm_d_g12d <- lm(age ~ vaf_pct, data = filter(d_df, variant_group == "KRAS G12D"))
lm_d_g12v <- lm(age ~ vaf_pct, data = filter(d_df, variant_group == "KRAS G12V"))
lm_d_int  <- lm(age ~ vaf_pct * variant_group, data = d_df)
int_p_d   <- summary(lm_d_int)$coefficients[grep(":", rownames(summary(lm_d_int)$coefficients)), 4]
sp_d      <- cor.test(d_df$vaf_pct, d_df$age, method = "spearman", exact = FALSE)
# Per-variant Spearman ρ (audit 12: rank-based test disagrees with the
# OLS slope for G12D — surface so the prose can disclose the
# disagreement honestly).
sp_d_g12d <- cor.test(filter(d_df, variant_group == "KRAS G12D")$vaf_pct,
                      filter(d_df, variant_group == "KRAS G12D")$age,
                      method = "spearman", exact = FALSE)
sp_d_g12v <- cor.test(filter(d_df, variant_group == "KRAS G12V")$vaf_pct,
                      filter(d_df, variant_group == "KRAS G12V")$age,
                      method = "spearman", exact = FALSE)
# Cook's distance audit on the G12D fit. Per audit 12, dropping the
# single most-influential lesion pushes the OLS P from ~0.044 to
# ~0.064 — surface the distance threshold cross-count so a sensitivity
# line can be cited.
.cook_g12d <- cooks.distance(lm_d_g12d)
.cook_n_g12d_threshold <- sum(.cook_g12d > 4 / length(.cook_g12d), na.rm = TRUE)

cat(sprintf("  n=%d | G12D slope=%.2f yrs/%%VAF (p=%.3g) | G12V slope=%.2f (p=%.3g)\n",
            nrow(d_df),
            coef(lm_d_g12d)[2], summary(lm_d_g12d)$coefficients[2, 4],
            coef(lm_d_g12v)[2], summary(lm_d_g12v)$coefficients[2, 4]))
cat(sprintf("  interaction p = %.3g | Spearman rho = %.2f (p = %.3g)\n",
            int_p_d, sp_d$estimate, sp_d$p.value))

# v6.30 (2026-05-20): use the canonical vaf_age_scatter_panel() helper
# (utils.R) so Fig 1 Panel G is one call site and ED07 panel A is the
# inverted-axes call site — both share the same grammar 1-to-1.
p4d <- vaf_age_scatter_panel(
  d_df,
  x_var     = "age",
  y_var     = "vaf_pct",
  color_var = "variant_group",
  palette   = PAL_KM,
  vaf_axis  = "y",
  x_lab     = "Age at presentation (years)",
  y_lab     = "VAF (%)"
) + theme(legend.position = "bottom")

fig3e_dir <- panel_slot_dir("vaf_age_scatter")
# Spacer removed 2026-05-16 — see p4b save above for rationale.
save_trio(fig3e_dir, "vaf_age_scatter", p4d, w = 7, h = 5)

# ─────────────────────────────────────────────────────────────────────────────
# Dump Fig 4 B/C/D stats
# ─────────────────────────────────────────────────────────────────────────────
writeLines(c(
  sprintf("# %s — age-density & VAF×age supporting statistics",
          panel_prose_tag(c("age_density_variant", "vaf_age_scatter"))),
  sprintf("Generated %s", Sys.Date()),
  "",
  "## Age density, Mut+ vs Negative (binary; panel dropped from final figure)",
  sprintf("n=%d  |  median Mut+=%.1f  median Neg=%.1f  |  KW p=%.3g",
          nrow(b_df), med_b$med[1], med_b$med[2], kw_b$p.value),
  "",
  sprintf("## %s — Age density by variant (G12D / G12V / Negative)",
          panel_prose_tag("age_density_variant")),
  sprintf("n=%d  |  KW p=%.3g", nrow(c_df), kw_c$p.value),
  capture.output(print(med_c)),
  "",
  sprintf("## %s — VAF × age (G12D / G12V, variant-positive)",
          panel_prose_tag("vaf_age_scatter")),
  sprintf("n=%d", nrow(d_df)),
  sprintf("G12D slope = %.2f yrs per %%VAF (p=%.3g)",
          coef(lm_d_g12d)[2], summary(lm_d_g12d)$coefficients[2, 4]),
  sprintf("G12V slope = %.2f yrs per %%VAF (p=%.3g)",
          coef(lm_d_g12v)[2], summary(lm_d_g12v)$coefficients[2, 4]),
  sprintf("Interaction (slope × variant) p = %.3g", int_p_d),
  sprintf("Spearman (pooled) rho = %.2f, p = %.3g", sp_d$estimate, sp_d$p.value)),
  file.path(stats_dir, "fig3_bcd_stats.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 composite — panels a–e for the age-story main figure
#
# Layout (per Apr 21 author direction; panel letters re-assigned 2026-04-22
# so first-citation order is strictly ascending — the
# rupture KM is discussed before the VAF × age dose-response, so the
# rupture KM becomes panel d and VAF × age becomes panel e):
#
#     ┌──────────────────┬────────────────┐
#     │                  │       b        │
#     │   a  (KM pres)   ├────────────────┤
#     │   + at-risk tbl  │       c        │
#     │                  ├────────────────┤
#     │                  │                │
#     ├──────────────────┤       e        │
#     │                  │   (VAF × age)  │
#     │   d  (KM rupt)   │                │
#     │   + at-risk tbl  │                │
#     │                  │                │
#     └──────────────────┴────────────────┘
#
# Both KM plots (a, d) occupy the left column at equal size, each with
# their at-risk table beneath them so reviewers can check stratum-level
# follow-up without flipping to the Extended Data. Panels b, c, e stack
# on the right.
#
# Tags a–e are added via per-panel labs(tag=) so they survive the mixed
# composite (patchwork tag_levels walk was unreliable when one cell was
# itself a stacked plot+table sub-patchwork).
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Fig 2 composite: deferred to 26_F1_assemble.R ──\n")
cat("  per-panel RDS files saved by save_trio() are read by\n")
cat("  analysis/01_main_analysis/26_F1_assemble.R, which assembles all 7 panels\n")
cat("  (a-e from this script, f & g from 21_F3_rupture_score_panels.R) at native\n")
cat("  Nature Medicine double-column footprint with theme_avm_native()\n")
cat("  applied uniformly. No publisher-stage resize needed.\n")

# The legacy inline-patchwork composite block below built only 5 of 7
# registry panels (a-e). It's been removed — see 26_F1_assemble.R for
# the canonical Fig 2 composite. The objects p4a/peff/p4b/p4c/p4d are
# still needed earlier in this script for stats dumps and per-panel
# saves; nothing further uses them after this point.

# ─────────────────────────────────────────────────────────────────────────────
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# Stats manifest fragments \u2014 fig3 + edfig16 + edfig17
# -----------------------------------------------------------------------------
# Everything the Results \u00a73/\u00a74 prose cites is emitted here into named scalar
# keys (the stats manifest).
# live at render time. Full tables (summary_table, pw_matrix, cox_coef) are
# kept alongside the scalars so downstream ED captions can reconstruct the
# same numbers without re-running the KM fits.
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

# Per-stratum summaries (survfit summary tables) \u2014 these helpers keep row
# lookup robust against survfit's "variant_group=KRAS G12D" label format.
# NB: `as.list(mat[row, , drop=FALSE])` drops the column names, so we have
# to `as.list(mat[row, ])` on a named numeric vector instead.
.km_row <- function(summ_table, stratum_label) {
  rownames(summ_table) <- sub("^.*=\\s*", "", rownames(summ_table))
  as.list(summ_table[stratum_label, ])
}

# Pairwise p-values (Bonferroni) \u2014 survdiff returns a lower-tri matrix
# indexed by factor levels. Both orderings appear depending on level order,
# so look up symmetrically and return whichever cell is non-NA.
.pw_get <- function(pw_mat, a, b) {
  m <- as.matrix(pw_mat)
  for (pair in list(c(a, b), c(b, a))) {
    if (pair[1] %in% rownames(m) && pair[2] %in% colnames(m)) {
      v <- m[pair[1], pair[2]]
      if (!is.na(v)) return(unname(v))
    }
  }
  NA_real_
}

fig3a_tbl  <- summary(fit1)$table
fig3a_g12d <- .km_row(fig3a_tbl, "KRAS G12D")
fig3a_g12v <- .km_row(fig3a_tbl, "KRAS G12V")
fig3a_neg  <- .km_row(fig3a_tbl, "Negative")

# Panel d (rupture) Cox coefficients for the "Negative" stratum vs reference.
# Three coexisting model specifications:
#   * unadj  — `Surv(age, event) ~ km2_group`               (HEADLINE)
#   * minadj — + sex_f + sample_type_clean                  (sensitivity #1)
#   * cox2   — + sex_f + study_clean + sample_type_clean    (sensitivity #2)
.extract_neg <- function(fit, label) {
  coefs <- summary(fit)$coefficients
  conf  <- summary(fit)$conf.int
  row <- grep("km2_groupNegative", rownames(coefs), value = TRUE)
  if (length(row) != 1L) stop(sprintf("Cox row for km2_groupNegative not found (%s)", label))
  list(hr = unname(conf[row, "exp(coef)"]),
       ci_lo = unname(conf[row, "lower .95"]),
       ci_hi = unname(conf[row, "upper .95"]),
       p = unname(coefs[row, "Pr(>|z|)"]))
}
.unadj  <- .extract_neg(cox_unadj,  "unadj")
.minadj <- .extract_neg(cox_minadj, "minadj")
.full   <- .extract_neg(cox2,       "full")
.angio  <- .extract_neg(cox_angio,  "angio")
rupt_cox_angio_n          <- nrow(cox_angio_df)
rupt_cox_angio_hr_neg     <- .angio$hr
rupt_cox_angio_ci_lo_neg  <- .angio$ci_lo
rupt_cox_angio_ci_hi_neg  <- .angio$ci_hi
rupt_cox_angio_p_neg      <- .angio$p
# G12D-vs-panel-negative direction = reciprocal of the Negative-vs-G12D-reference
# coefficient above (CI bounds swap on inversion; P is invariant). The §1 prose
# reports this so an HR > 1 reads consistently with "KRAS^G12D^ associated with
# earlier rupture" — the reference-level coefficient (hr_neg < 1) would otherwise
# imply G12D is protective.
rupt_cox_angio_hr_g12d    <- 1 / rupt_cox_angio_hr_neg
rupt_cox_angio_ci_lo_g12d <- 1 / rupt_cox_angio_ci_hi_neg
rupt_cox_angio_ci_hi_g12d <- 1 / rupt_cox_angio_ci_lo_neg
rupt_cox_angio_p_g12d     <- rupt_cox_angio_p_neg
cat("  Angioarchitectural-adjusted Cox:\n"); print(summary(cox_angio)$coefficients)
# Keep cox2_coef / cox2_conf in scope for downstream stats-fragment
# exporters that ship the full coefficient table.
cox2_coef <- summary(cox2)$coefficients
cox2_conf <- summary(cox2)$conf.int

rupt_cox_unadj_hr_neg     <- .unadj$hr
rupt_cox_unadj_ci_lo_neg  <- .unadj$ci_lo
rupt_cox_unadj_ci_hi_neg  <- .unadj$ci_hi
rupt_cox_unadj_p_neg      <- .unadj$p

rupt_cox_minadj_hr_neg    <- .minadj$hr
rupt_cox_minadj_ci_lo_neg <- .minadj$ci_lo
rupt_cox_minadj_ci_hi_neg <- .minadj$ci_hi
rupt_cox_minadj_p_neg     <- .minadj$p

# Legacy full-adjustment keys (sensitivity #2; covariate-form study adjustment)
rupt_cox_hr_neg    <- .full$hr
rupt_cox_ci_lo_neg <- .full$ci_lo
rupt_cox_ci_hi_neg <- .full$ci_hi
rupt_cox_p_neg     <- .full$p

fig3d_tbl  <- summary(fit2)$table
fig3d_g12d <- .km_row(fig3d_tbl, "KRAS G12D")
fig3d_g12v <- .km_row(fig3d_tbl, "KRAS G12V")
fig3d_neg  <- .km_row(fig3d_tbl, "Negative")

# Pediatric (<18) Fisher OR \u2014 cited in \u00a74 alongside Panel b/c. Computed on
# the same age-available binary cohort as Panel b.
ped_bin <- b_df %>% mutate(ped = age < 18)
ped_tbl <- table(
  geno = factor(ped_bin$geno_binary,
                levels = c("Variant-positive", "Panel-negative")),
  ped  = factor(ped_bin$ped, levels = c(TRUE, FALSE))
)
ped_fisher <- fisher.test(ped_tbl)

fig3_fragment <- list(
  # Panel a \u2014 KM time to any clinical presentation
  km_pres_n              = nrow(km1_df),
  km_pres_n_g12d         = sum(km1_df$variant_group == "KRAS G12D"),
  km_pres_n_g12v         = sum(km1_df$variant_group == "KRAS G12V"),
  km_pres_n_neg          = sum(km1_df$variant_group == "Negative"),
  km_pres_median_g12d    = unname(fig3a_g12d[["median"]]),
  km_pres_median_g12v    = unname(fig3a_g12v[["median"]]),
  km_pres_median_neg     = unname(fig3a_neg[["median"]]),
  km_pres_ci_lo_g12d     = unname(fig3a_g12d[["0.95LCL"]]),
  km_pres_ci_hi_g12d     = unname(fig3a_g12d[["0.95UCL"]]),
  km_pres_ci_lo_g12v     = unname(fig3a_g12v[["0.95LCL"]]),
  km_pres_ci_hi_g12v     = unname(fig3a_g12v[["0.95UCL"]]),
  km_pres_ci_lo_neg      = unname(fig3a_neg[["0.95LCL"]]),
  km_pres_ci_hi_neg      = unname(fig3a_neg[["0.95UCL"]]),
  km_pres_logrank_p      = lr1_p,
  km_pres_pw_g12d_neg    = .pw_get(pw1$p.value, "KRAS G12D", "Negative"),
  km_pres_pw_g12v_neg    = .pw_get(pw1$p.value, "KRAS G12V", "Negative"),
  km_pres_pw_g12d_g12v   = .pw_get(pw1$p.value, "KRAS G12D", "KRAS G12V"),
  # Panel b \u2014 age density, Mut+ vs Neg
  age_bin_n              = nrow(b_df),
  age_bin_n_mut          = sum(b_df$geno_binary == "Variant-positive"),
  age_bin_n_neg          = sum(b_df$geno_binary == "Panel-negative"),
  age_bin_median_mut     = unname(med_b$med[med_b$geno_binary == "Variant-positive"]),
  age_bin_median_neg     = unname(med_b$med[med_b$geno_binary == "Panel-negative"]),
  age_bin_kw_p           = kw_b$p.value,
  # Pediatric Fisher (cited in \u00a74 with Panel b/c)
  ped_n_mut              = sum(ped_tbl["Variant-positive", ]),
  ped_n_mut_u18          = ped_tbl["Variant-positive", "TRUE"],
  ped_n_neg              = sum(ped_tbl["Panel-negative", ]),
  ped_n_neg_u18          = ped_tbl["Panel-negative", "TRUE"],
  ped_or                 = unname(ped_fisher$estimate),
  ped_ci_lo              = unname(ped_fisher$conf.int[1]),
  ped_ci_hi              = unname(ped_fisher$conf.int[2]),
  ped_p                  = ped_fisher$p.value,
  # Panel c \u2014 age density by variant
  age_var_n              = nrow(c_df),
  age_var_median_g12d    = unname(med_c$med[med_c$variant_group == "KRAS G12D"]),
  age_var_median_g12v    = unname(med_c$med[med_c$variant_group == "KRAS G12V"]),
  age_var_median_neg     = unname(med_c$med[med_c$variant_group == "Negative"]),
  age_var_kw_p           = kw_c$p.value,
  # Panel e \u2014 VAF \u00d7 age
  vaf_age_n              = nrow(d_df),
  # Per-variant n's added 2026-04-26 to support the variant-specific
  # framing of the dose-response in the \u00a74 prose: G12D is the primary
  # finding; G12V is exploratory and presented as such.
  vaf_age_n_g12d         = sum(d_df$variant_group == "KRAS G12D"),
  vaf_age_n_g12v         = sum(d_df$variant_group == "KRAS G12V"),
  vaf_age_slope_g12d     = unname(coef(lm_d_g12d)[2]),
  vaf_age_slope_g12d_p   = unname(summary(lm_d_g12d)$coefficients[2, 4]),
  vaf_age_slope_g12v     = unname(coef(lm_d_g12v)[2]),
  vaf_age_slope_g12v_p   = unname(summary(lm_d_g12v)$coefficients[2, 4]),
  vaf_age_interaction_p  = unname(int_p_d),
  vaf_age_spearman_rho   = unname(sp_d$estimate),
  vaf_age_spearman_p     = sp_d$p.value,
  # Per-variant rank-based correlations + influence diagnostic
  # (audit 12: G12D Spearman ρ disagrees with the OLS slope; surface
  # explicitly so prose can honestly frame G12D as leverage-driven).
  vaf_age_spearman_rho_g12d  = unname(sp_d_g12d$estimate),
  vaf_age_spearman_p_g12d    = sp_d_g12d$p.value,
  vaf_age_spearman_rho_g12v  = unname(sp_d_g12v$estimate),
  vaf_age_spearman_p_g12v    = sp_d_g12v$p.value,
  vaf_age_g12d_cook_n_above4n = .cook_n_g12d_threshold,
  # Panel d \u2014 KM time to rupture
  rupt_n                 = nrow(km2_df),
  rupt_events_total      = sum(km2_df$km2_event == 1L),
  rupt_median_g12d       = unname(fig3d_g12d[["median"]]),
  rupt_median_g12v       = unname(fig3d_g12v[["median"]]),
  rupt_median_neg        = unname(fig3d_neg[["median"]]),
  rupt_ci_lo_g12d        = unname(fig3d_g12d[["0.95LCL"]]),
  rupt_ci_hi_g12d        = unname(fig3d_g12d[["0.95UCL"]]),
  rupt_ci_lo_g12v        = unname(fig3d_g12v[["0.95LCL"]]),
  rupt_ci_hi_g12v        = unname(fig3d_g12v[["0.95UCL"]]),
  rupt_ci_lo_neg         = unname(fig3d_neg[["0.95LCL"]]),
  rupt_ci_hi_neg         = unname(fig3d_neg[["0.95UCL"]]),
  rupt_logrank_p         = lr2_p,
  rupt_pw_g12d_neg       = .pw_get(pw2$p.value, "KRAS G12D", "Negative"),
  rupt_pw_g12v_neg       = .pw_get(pw2$p.value, "KRAS G12V", "Negative"),
  rupt_pw_g12d_g12v      = .pw_get(pw2$p.value, "KRAS G12D", "KRAS G12V"),
  # Unadjusted Cox — HEADLINE HR cited in Results §4 (pairs with the
  # log-rank / KM curves; the natural multivariate companion to the
  # univariate finding).
  rupt_cox_unadj_hr_neg     = rupt_cox_unadj_hr_neg,
  rupt_cox_unadj_ci_lo_neg  = rupt_cox_unadj_ci_lo_neg,
  rupt_cox_unadj_ci_hi_neg  = rupt_cox_unadj_ci_hi_neg,
  rupt_cox_unadj_p_neg      = rupt_cox_unadj_p_neg,
  # Minimal-adjustment sensitivity (sex + sample_type; no study)
  rupt_cox_minadj_hr_neg    = rupt_cox_minadj_hr_neg,
  rupt_cox_minadj_ci_lo_neg = rupt_cox_minadj_ci_lo_neg,
  rupt_cox_minadj_ci_hi_neg = rupt_cox_minadj_ci_hi_neg,
  rupt_cox_minadj_p_neg     = rupt_cox_minadj_p_neg,
  # Full-adjustment sensitivity (+ study_clean as fixed effect)
  rupt_cox_hr_neg        = rupt_cox_hr_neg,
  rupt_cox_ci_lo_neg     = rupt_cox_ci_lo_neg,
  rupt_cox_ci_hi_neg     = rupt_cox_ci_hi_neg,
  rupt_cox_p_neg         = rupt_cox_p_neg,
  rupt_cox_zph_global_p  = rupt_cox_zph_global_p,
  # Angioarchitectural-adjusted sensitivity (nidus size + DVDrainage + HR features)
  rupt_cox_angio_n           = rupt_cox_angio_n,
  rupt_cox_angio_hr_neg      = rupt_cox_angio_hr_neg,
  rupt_cox_angio_ci_lo_neg   = rupt_cox_angio_ci_lo_neg,
  rupt_cox_angio_ci_hi_neg   = rupt_cox_angio_ci_hi_neg,
  rupt_cox_angio_p_neg       = rupt_cox_angio_p_neg,
  # G12D-vs-Negative direction (reciprocal of the Negative-reference HR above),
  # cited in §1 so the reported HR > 1 matches "G12D associated with earlier rupture".
  rupt_cox_angio_hr_g12d     = rupt_cox_angio_hr_g12d,
  rupt_cox_angio_ci_lo_g12d  = rupt_cox_angio_ci_lo_g12d,
  rupt_cox_angio_ci_hi_g12d  = rupt_cox_angio_ci_hi_g12d,
  rupt_cox_angio_p_g12d      = rupt_cox_angio_p_g12d,
  # Strata(study_clean) sensitivity model — adopted as primary in the
  # 2026-04-29 audit because the covariate-form model violates PH.
  rupt_cox_strata_hr_neg        = rupt_cox_strata_hr_neg,
  rupt_cox_strata_ci_lo_neg     = rupt_cox_strata_ci_lo_neg,
  rupt_cox_strata_ci_hi_neg     = rupt_cox_strata_ci_hi_neg,
  rupt_cox_strata_p_neg         = rupt_cox_strata_p_neg,
  rupt_cox_strata_zph_global_p  = rupt_cox_strata_zph_global_p,
  # Full tables
  km_pres_summary_table  = as.data.frame(fig3a_tbl),
  km_pres_pw_matrix      = as.data.frame(pw1$p.value),
  rupt_summary_table     = as.data.frame(fig3d_tbl),
  rupt_cox_coef          = as.data.frame(cox2_coef)
)
write_stats_section(section = "fig3", stats = fig3_fragment)

# ---- edfig16 \u2014 KM rupture ------------------------------------------------
edfig16_fragment <- list(
  n                  = nrow(km2_df),
  n_events           = sum(km2_df$km2_event == 1L),
  n_censored         = sum(km2_df$km2_event == 0L),
  median_g12d        = unname(fig3d_g12d[["median"]]),
  median_g12v        = unname(fig3d_g12v[["median"]]),
  median_neg         = unname(fig3d_neg[["median"]]),
  logrank_p          = lr2_p,
  pw_g12d_neg        = .pw_get(pw2$p.value, "KRAS G12D", "Negative"),
  pw_g12v_neg        = .pw_get(pw2$p.value, "KRAS G12V", "Negative"),
  pw_g12d_g12v       = .pw_get(pw2$p.value, "KRAS G12D", "KRAS G12V"),
  cox_hr_neg         = rupt_cox_hr_neg,
  cox_ci_lo_neg      = rupt_cox_ci_lo_neg,
  cox_ci_hi_neg      = rupt_cox_ci_hi_neg,
  cox_p_neg          = rupt_cox_p_neg,
  summary_table      = as.data.frame(fig3d_tbl),
  cox_coef           = as.data.frame(cox2_coef)
)
write_stats_section(section = "edfig16", stats = edfig16_fragment)

# ---- edfig17 \u2014 rare-variant KM -------------------------------------------
fig17_tbl  <- summary(fit3)$table
fig17_ok   <- .km_row(fig17_tbl, "Other KRAS")
fig17_braf <- .km_row(fig17_tbl, "BRAF V600E")
fig17_neg  <- .km_row(fig17_tbl, "Negative")

edfig17_fragment <- list(
  n                  = nrow(km3_df),
  n_other_kras       = sum(km3_df$rare_group == "Other KRAS"),
  n_braf_v600e       = sum(km3_df$rare_group == "BRAF V600E"),
  n_neg              = sum(km3_df$rare_group == "Negative"),
  median_other_kras  = unname(fig17_ok[["median"]]),
  median_braf_v600e  = unname(fig17_braf[["median"]]),
  median_neg         = unname(fig17_neg[["median"]]),
  logrank_p          = lr3_p,
  pw_other_kras_neg  = .pw_get(pw3$p.value, "Other KRAS",  "Negative"),
  pw_braf_v600e_neg  = .pw_get(pw3$p.value, "BRAF V600E",  "Negative"),
  pw_other_kras_braf = .pw_get(pw3$p.value, "Other KRAS",  "BRAF V600E"),
  summary_table      = as.data.frame(fig17_tbl)
)
write_stats_section(section = "edfig17", stats = edfig17_fragment)

cat("\n\u2713 Figure 2 (A-D panels) + Fig 2 composite + ExtendedData 16/17 complete.\n")
cat(sprintf("   Fig 2A:  %s\n", file.path(fig3_dir,  "km_presentation.png")))
cat("   Fig 2B:  age_density_binary (panel dropped 2026-05-19; not saved)\n")
cat(sprintf("   Fig 2C:  %s\n", file.path(fig3c_dir, "age_density_variant.png")))
cat(sprintf("   Fig 2D:  %s\n", file.path(here("results", "Figure2", "panel_D"), "km_rupture.png")))
cat(sprintf("   Fig 2E:  %s\n", file.path(fig3e_dir, "vaf_age_scatter.png")))
cat("   Fig 2 composite: built by 26_F1_assemble.R after this script + 21_F3_rupture_score_panels.R\n")
cat(sprintf("   rare_variants ED panel: %s\n", file.path(efig17_dir, "rare_variants.png")))
