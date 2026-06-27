# 20_ST_supp_tables_8_9_10.R — Generate SuppTable 08, 09, 10 spreadsheets.
#
# These tables back the final Results section (Sections 2 + 4) of the
# manuscript and are part of the manuscript.
#
# Inputs (already produced by upstream scripts):
#   * results/stats/fig3_stats.rds        — 14_ED7_vaf_phenotype.R  (VAF × phenotype)
#   * results/stats/fig4_stats.rds        — (age × phenotype)
#   * data/processed/bAVM_analysis_ready.rds — analysis-ready dataset
#
# Outputs:
#   * results/SupplementaryTables/SuppTable06_vaf_phenotype_correlations.xlsx
#   * results/SupplementaryTables/SuppTable04_age_phenotype_interactions.xlsx
#   * results/SupplementaryTables/SuppTable07_vaf_age_sensitivity.xlsx
#
# Conventions:
#   * Six pooled phenotypes:
#       SM total grade, SM size, drainage, eloquence, rupture, high-risk count.
#   * Spearman ρ + BH-FDR on continuous phenotypes; Wilcoxon / logistic β on
#     binary phenotypes. Per-variant subsheets restrict to the same panel so
#     readers can audit per-phenotype × per-variant cells.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(here)
  library(writexl); library(openxlsx); library(broom); library(logistf)
})

source(here("analysis", "helper_scripts", "utils.R"))

# ── Paths ────────────────────────────────────────────────────────────────────
etable_dir <- here("results", "SupplementaryTables")
dir.create(etable_dir, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(mutation_positive))

# Variant classification aligned with 09_F1_km_age.R
df <- df %>%
  mutate(
    variant_group = case_when(
      !mutation_positive              ~ "Negative",
      mutation == "KRAS G12D"         ~ "KRAS G12D",
      mutation == "KRAS G12V"         ~ "KRAS G12V",
      mutation_gene == "KRAS"          ~ "Other KRAS",
      mutation_gene == "BRAF"          ~ "BRAF",
      TRUE                             ~ "Negative"
    ),
    variant_group = factor(variant_group,
      levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF", "Negative"))
  )

# ═════════════════════════════════════════════════════════════════════════════
# Supplementary Table 5 — VAF × phenotype correlations
# ═════════════════════════════════════════════════════════════════════════════
# Six pooled phenotypes (match Extended Data Fig. 5 panels a–f) plus per-variant subsheets.
# Continuous outcomes report Spearman ρ and BH-FDR-adjusted P; binary outcomes
# use Wilcoxon W with logistic β from GLM for effect direction.

cat("\n══ Supp Table 8 — VAF × phenotype correlations ══\n")

vaf_phenos <- tribble(
  ~var,                ~label,                           ~type,
  "sm_total_num",      "Spetzler\u2013Martin total grade",   "continuous",
  "sm_size_num",       "Spetzler\u2013Martin size",           "continuous",
  "sm_drainage_num",   "Deep venous drainage",            "binary",
  "sm_eloquence_num",  "Eloquence",                       "binary",
  "ever_ruptured_num", "Rupture status (ever)",            "binary",
  "n_high_risk_num",   "High-risk feature count",         "continuous"
)

vaf_phen_fit <- function(sub, phenos = vaf_phenos) {
  map_dfr(seq_len(nrow(phenos)), function(i) {
    o <- phenos[i, ]
    s <- sub %>% filter(!is.na(.data[[o$var]]), !is.na(vaf_prop))
    if (nrow(s) < 10) {
      return(tibble(phenotype = o$label, n = nrow(s),
                    statistic = NA_character_, estimate = NA_real_,
                    p_value = NA_real_))
    }
    if (o$type == "continuous") {
      sp <- suppressWarnings(
        cor.test(s$vaf_prop, s[[o$var]], method = "spearman", exact = FALSE))
      tibble(phenotype = o$label, n = nrow(s),
             statistic = "Spearman \u03c1",
             estimate = unname(sp$estimate), p_value = sp$p.value)
    } else {
      m <- suppressWarnings(
        glm(as.formula(paste(o$var, "~ vaf_prop")),
            data = s, family = binomial))
      cs <- summary(m)$coefficients
      vr <- grep("vaf_prop", rownames(cs))
      if (length(vr) == 0) {
        return(tibble(phenotype = o$label, n = nrow(s),
                      statistic = NA_character_, estimate = NA_real_,
                      p_value = NA_real_))
      }
      tibble(phenotype = o$label, n = nrow(s),
             statistic = "Logistic \u03b2 (VAF, proportion)",
             estimate = cs[vr, 1], p_value = cs[vr, 4])
    }
  }) %>%
    mutate(p_fdr = p.adjust(p_value, method = "BH"))
}

vaf_data <- df %>% filter(mutation_positive, !is.na(vaf_prop))

st8_pooled <- vaf_phen_fit(vaf_data)
st8_g12d   <- vaf_phen_fit(vaf_data %>% filter(variant_group == "KRAS G12D"))
st8_g12v   <- vaf_phen_fit(vaf_data %>% filter(variant_group == "KRAS G12V"))
st8_other  <- vaf_phen_fit(vaf_data %>%
                             filter(variant_group %in% c("Other KRAS", "BRAF")))

cat("  Pooled (n=", nrow(vaf_data), ")\n", sep = "")
print(st8_pooled)

# 2026-06-14: ST6 is now a single-sheet mega-table written by
# 18_ED10_age_adj_vaf_pheno.R, stacking these per-cohort correlations on top
# of the age-adjusted regression block. Persist the correlation rows to a
# private intermediate RDS that 15 reads \u2014 separates compute (here) from
# rendering (there) and avoids the old write_supp_table_sheets cross-script
# accumulation that left a stale duplicate "Pooled mutation-positive"
# sheet for weeks.
st6_corr_path <- here("results", "stats", "_intermediates",
                      "st6_correlations.rds")
dir.create(dirname(st6_corr_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    pooled_variant_positive = st8_pooled,
    kras_g12d               = st8_g12d,
    kras_g12v               = st8_g12v,
    other_kras_braf         = st8_other
  ),
  st6_corr_path
)
cat(sprintf("  \u2192 correlation block cached: %s\n",
            sub(paste0(here::here(), "/"), "", st6_corr_path)))

# ═════════════════════════════════════════════════════════════════════════════
# Supplementary Table 9 — Age × phenotype interaction tests
# ═════════════════════════════════════════════════════════════════════════════
# Tests whether the mut+/neg phenotype association is modified by age group
# (pediatric vs. adult).
#   phenotype ~ geno_binary * age_group
# (numeric phenotypes use OLS; binary phenotypes use logistic GLM). This is
# the specification cited in §4 prose ("all interaction *P* ≥ 0.079").
#
# A second sheet flips the framing (age ~ geno * phenotype) for reviewers
# who wish to audit the alternative "does the age shift differ by phenotype"
# question directly.

cat("\n══ Supp Table 9 — Age × phenotype interaction tests ══\n")

age_phenos <- tribble(
  ~var,                ~label,                           ~type,
  "sm_size_num",       "Spetzler\u2013Martin size",           "continuous",
  "sm_drainage_num",   "Deep venous drainage",            "binary",
  "ever_ruptured_num", "Rupture status (ever)",            "binary",
  "n_high_risk_num",   "High-risk feature count",         "continuous"
)

# Sheet 1 — pheno ~ geno * age_group (prose 0.079)
gt <- df %>%
  filter(!is.na(age), !is.na(mutation_positive)) %>%
  mutate(
    geno_binary = factor(
      ifelse(mutation_positive, "Variant-positive", "Panel-negative"),
      levels = c("Panel-negative", "Variant-positive")),
    age_group = factor(ifelse(age < 18, "Pediatric", "Adult"),
                       levels = c("Adult", "Pediatric"))
  )

# Shared extractor used by both interaction sheets. Returns a one-row tibble
# with phenotype, n, interaction_beta/se/p from the coef row whose name
# contains ":" (the interaction term). Returns NA row on small-n, model
# error, or absent interaction term.
#
# 2026-06-14 fix: logistf objects must be handled separately. summary() on
# an lm/glm returns a 4-column matrix with rownames in $coefficients;
# summary() on a logistf prints the model and exposes m$coefficients as a
# NAMED NUMERIC VECTOR (the point estimates). The old code then called
# rownames(cs) on a vector → NULL → grep(":", NULL) → integer(0) → all-NA
# row. This silently zeroed the two binary phenotypes (deep venous
# drainage + rupture ever) in the §4 "all interaction P ≥ 0.079" table,
# even though logistf had fit them cleanly.
.extract_int_row <- function(s, label, m) {
  na_row <- tibble(phenotype = label, n = nrow(s),
                   interaction_beta = NA_real_, interaction_se = NA_real_,
                   interaction_p = NA_real_)
  if (is.null(m)) return(na_row)

  if (inherits(m, "logistf")) {
    # logistf — pull β / SE / p directly from the fit object.
    coefs <- m$coefficients              # named numeric vector
    ses   <- sqrt(diag(m$var))           # SE = sqrt(diag(var-cov))
    ps    <- m$prob                      # Wald p-values, named
    int_idx <- grep(":", names(coefs))
    if (length(int_idx) == 0) return(na_row)
    return(tibble(
      phenotype        = label,
      n                = nrow(s),
      interaction_beta = unname(coefs[int_idx[1]]),
      interaction_se   = unname(ses[int_idx[1]]),
      interaction_p    = unname(ps[int_idx[1]])
    ))
  }

  # lm / glm path — summary()$coefficients is a 4-col matrix.
  cs      <- summary(m)$coefficients
  int_row <- grep(":", rownames(cs))
  if (length(int_row) == 0) return(na_row)
  tibble(
    phenotype        = label,
    n                = nrow(s),
    interaction_beta = cs[int_row[1], 1],
    interaction_se   = cs[int_row[1], 2],
    interaction_p    = cs[int_row[1], 4]
  )
}

# Sheet 1 — pheno ~ geno * age_group
fit_int <- function(var_name, label, type, data) {
  s <- data %>% filter(!is.na(.data[[var_name]]), !is.na(age_group))
  if (nrow(s) < 30) {
    return(tibble(phenotype = label, n = nrow(s),
                  interaction_beta = NA_real_, interaction_se = NA_real_,
                  interaction_p = NA_real_))
  }
  m <- tryCatch({
    if (type == "continuous") {
      lm(as.formula(paste(var_name, "~ geno_binary * age_group")), data = s)
    } else {
      # Firth penalised logistic for binary outcomes — handles near-separation
      # in sparse pediatric/adult strata without silent NULL failures.
      logistf(as.formula(paste(var_name, "~ geno_binary * age_group")),
              data = s, pl = FALSE)
    }
  }, error = function(e) NULL)
  .extract_int_row(s, label, m)
}

# Sheet 2 — age ~ geno * pheno (flipped framing for reviewer audit).
# Always OLS: age is the continuous outcome; age is never NA in gt (filtered
# at construction). No `type` argument — formula direction is fixed.
fit_int_flipped <- function(var_name, label, data) {
  s <- data %>% filter(!is.na(.data[[var_name]]))
  if (nrow(s) < 30) {
    return(tibble(phenotype = label, n = nrow(s),
                  interaction_beta = NA_real_, interaction_se = NA_real_,
                  interaction_p = NA_real_))
  }
  m <- tryCatch(
    lm(as.formula(paste("age ~ geno_binary *", var_name)), data = s),
    error = function(e) NULL
  )
  .extract_int_row(s, label, m)
}

st9_primary <- pmap_dfr(
  list(age_phenos$var, age_phenos$label, age_phenos$type),
  fit_int, data = gt
) %>%
  mutate(interaction_p_fdr = p.adjust(interaction_p, method = "BH"))

print(st9_primary)

st9_flipped <- pmap_dfr(
  list(age_phenos$var, age_phenos$label),
  fit_int_flipped, data = gt
) %>%
  mutate(interaction_p_fdr = p.adjust(interaction_p, method = "BH"))

# Sheet 1 \u2014 angio-adjusted multivariable rupture Cox (the model cited in \u00a71).
# Genotype rows are vs panel-negative; the G12D HR equals
# stats$fig3$rupt_cox_angio_hr_g12d (1.68) by construction (same cohort/model
# as cox_angio in 09_F1_km_age.R). The VIF column (~1.0) documents that the three
# adjustment covariates are non-collinear.
st4_cox <- rupture_angio_cox_table(df)
cat(sprintf("  Rupture Cox (angio-adjusted): n=%d, events=%d\n",
            attr(st4_cox, "model_n"), attr(st4_cox, "model_events")))
print(st4_cox)

# \u2500\u2500 2026-06-15: ST4 back to 3 sheets (Nature MOESM3 (A)/(B)/(C) pattern) \u2500\u2500\u2500\u2500
# Each block lives on its own tab so per-sheet schemas can be sheet-natural
# (Cox sheet carries VIF, no FDR; interaction sheets carry FDR, no VIF).
# The Analysis column is dropped \u2014 the sheet name carries that information.
source(here("analysis", "helper_scripts", "supp_table_writer.R"))

.fmt_p <- function(p) {
  ifelse(p < 0.001,
         formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

cox_n      <- attr(st4_cox, "model_n")
cox_events <- attr(st4_cox, "model_events")

cox_sheet_df <- st4_cox %>%
  transmute(
    term_pheno    = Term,
    n             = as.integer(cox_n),
    estimate      = sprintf("%.2f", HR),
    ci_str        = sprintf("(%.2f, %.2f)", `CI low`, `CI high`),
    p_str         = .fmt_p(P),
    vif_str       = sprintf("%.2f", VIF)
  )

.interaction_sheet <- function(d) {
  d %>%
    transmute(
      term_pheno    = as.character(phenotype),
      n             = as.integer(n),
      estimate      = sprintf("%.3f", interaction_beta),
      ci_str        = sprintf("(%.3f, %.3f)",
                              interaction_beta - 1.96 * interaction_se,
                              interaction_beta + 1.96 * interaction_se),
      p_str         = .fmt_p(interaction_p),
      fdr_p_str     = .fmt_p(interaction_p_fdr)
    )
}

st4_sheets <- list(
  "(A) Rupture Cox (angio-adj)" = list(
    data = cox_sheet_df,
    columns = list(
      col    ("term_pheno",  label = "Term"),
      col_int("n",           label = "N",        italic = TRUE),
      col    ("estimate",    label = "HR",       italic = TRUE),
      col    ("ci_str",      label = "95% CI"),
      col    ("p_str",       label = "P",        italic = TRUE),
      col    ("vif_str",     label = "VIF",      italic = TRUE)
    ),
    footnote = c(
      sprintf("Multivariable Cox proportional-hazards model for rupture, angioarchitecture-adjusted. Model n = %d, events = %d.",
              cox_n, cox_events),
      "Genotype rows are vs. panel-negative; the KRAS G12D HR matches the prose value (1.68) by construction. VIF (~1) documents non-collinearity of the three angioarchitectural covariates."
    )
  ),
  "(B) Pheno x Geno x AgeGroup" = list(
    data = .interaction_sheet(st9_primary),
    columns = list(
      col    ("term_pheno",  label = "Phenotype"),
      col_int("n",           label = "N",          italic = TRUE),
      col    ("estimate",    label = "\u03b2",     italic = TRUE),
      col    ("ci_str",      label = "95% CI"),
      col    ("p_str",       label = "P",          italic = TRUE),
      col    ("fdr_p_str",   label = "FDR P",      italic = TRUE)
    ),
    footnote = c(
      "Interaction term from phenotype ~ genotype \u00d7 age-group. Continuous outcomes (SM size, high-risk feature count) use OLS; binary outcomes (deep venous drainage, rupture status) use Firth-penalized logistic regression (logistf) to avoid near-separation in sparse pediatric / adult strata.",
      "\u03b2 \u00b1 1.96\u00b7SE Wald CIs. N is per-phenotype complete-case count after dropping rows with NA in the phenotype or age. FDR P uses BH correction within this sheet (4 rows)."
    )
  ),
  "(C) Age x Geno x Phenotype" = list(
    data = .interaction_sheet(st9_flipped),
    columns = list(
      col    ("term_pheno",  label = "Phenotype"),
      col_int("n",           label = "N",          italic = TRUE),
      col    ("estimate",    label = "\u03b2",     italic = TRUE),
      col    ("ci_str",      label = "95% CI"),
      col    ("p_str",       label = "P",          italic = TRUE),
      col    ("fdr_p_str",   label = "FDR P",      italic = TRUE)
    ),
    footnote = c(
      "Reviewer-audit framing \u2014 interaction term from age ~ genotype \u00d7 phenotype (always OLS since age is the continuous outcome). Surfaces the alternative 'does the age shift differ by phenotype' question alongside the primary 'does the phenotype association differ by age' framing in sheet (B).",
      "FDR P uses BH correction within this sheet (4 rows)."
    )
  )
)

write_supp_table_workbook(
  sheets = st4_sheets,
  path   = file.path(etable_dir, "SuppTable04_age_phenotype_interactions.xlsx")
)
cat(sprintf("  \u2713 SuppTable04_age_phenotype_interactions.xlsx (3 sheets: %d + %d + %d rows)\n",
            nrow(st4_sheets[["(A) Rupture Cox (angio-adj)"]]$data),
            nrow(st4_sheets[["(B) Pheno x Geno x AgeGroup"]]$data),
            nrow(st4_sheets[["(C) Age x Geno x Phenotype"]]$data)))

# ═════════════════════════════════════════════════════════════════════════════
# Supplementary Table 10 — VAF × age dose-response sensitivity
# ═════════════════════════════════════════════════════════════════════════════
# Three sheets, all run on the KRAS G12D + KRAS G12V union (variant-positive,
# VAF observed, age observed):
#   Sheet 1: Leave-one-series-out (eight contributing cohorts).
#   Sheet 2: VAF outlier trimming at the 90th, 95th, and 99th percentiles.
#   Sheet 3: Cohort-restriction (institutional / published / full).

cat("\n══ Supp Table 10 — VAF × age sensitivity ══\n")

vaf_age_df <- df %>%
  filter(mutation_positive, !is.na(vaf_prop), !is.na(age),
         variant_group %in% c("KRAS G12D", "KRAS G12V")) %>%
  mutate(variant_group = droplevels(variant_group),
         vaf_pct = vaf_prop * 100,
         cohort_type = ifelse(
           study_clean %in% c("BCH", "UAB", "CHOP"),
           "Institutional", "Published"))

fit_slope <- function(sub, scenario_label) {
  if (nrow(sub) < 10 || n_distinct(sub$variant_group) < 1) {
    return(tibble(scenario = scenario_label, n = nrow(sub),
                  slope = NA_real_, se = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_,
                  p_value = NA_real_))
  }
  m <- lm(age ~ vaf_pct, data = sub)
  cf <- summary(m)$coefficients
  ci <- tryCatch(confint(m, "vaf_pct", level = 0.95),
                 error = function(e) matrix(c(NA_real_, NA_real_), ncol = 2))
  tibble(
    scenario = scenario_label,
    n        = nrow(sub),
    slope    = cf["vaf_pct", "Estimate"],
    se       = cf["vaf_pct", "Std. Error"],
    ci_low   = ci[1, 1],
    ci_high  = ci[1, 2],
    p_value  = cf["vaf_pct", "Pr(>|t|)"]
  )
}

# Sheet 1 — Leave-one-series-out ------------------------------------------
series_list <- sort(unique(vaf_age_df$study_clean))
st10_loso <- map_dfr(series_list, function(exc) {
  sub <- vaf_age_df %>% filter(study_clean != exc)
  fit_slope(sub, sprintf("Excluding %s", exc))
})
# Prepend full-cohort baseline for reference
st10_loso <- bind_rows(
  fit_slope(vaf_age_df, "Full cohort (baseline)"),
  st10_loso
)
cat("  LOSO scenarios: ", nrow(st10_loso), "\n", sep = "")

# Sheet 2 — VAF outlier trimming ------------------------------------------
q_levels <- c(0.90, 0.95, 0.99)
st10_trim <- map_dfr(q_levels, function(q) {
  cut <- quantile(vaf_age_df$vaf_pct, q, na.rm = TRUE)
  sub <- vaf_age_df %>% filter(vaf_pct <= cut)
  fit_slope(sub, sprintf("Trim \u2264 %.0fth pctile (VAF %.2f%%)",
                         q * 100, cut))
})
st10_trim <- bind_rows(
  fit_slope(vaf_age_df, "No trim (baseline)"),
  st10_trim
)

# Sheet 3 — Cohort restriction --------------------------------------------
st10_cohort <- bind_rows(
  fit_slope(vaf_age_df,                                      "Full cohort"),
  fit_slope(vaf_age_df %>% filter(cohort_type == "Institutional"),
            "Institutional only (BCH + UAB + CHOP)"),
  fit_slope(vaf_age_df %>% filter(cohort_type == "Published"),
            "Published series only")
)

cat("\n  LOSO table:\n"); print(st10_loso)
cat("\n  Trim table:\n"); print(st10_trim)
cat("\n  Cohort table:\n"); print(st10_cohort)

# \u2500\u2500 2026-06-14: ST7 collapsed from 7 sheets to a single mega-sheet \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# Read 16_ST_vaf_outlier_sensitivity.R's cached frames (outlier sensitivity + dose-
# response models + LRT sensitivity + distribution metadata) and combine
# with the LOSO / VAF-trim / cohort-restriction blocks built above into one
# long-format sheet keyed on (Analysis, Scenario, Outcome).
st7_int_path <- here("results", "stats", "_intermediates",
                     "st7_outlier_dose.rds")
if (!file.exists(st7_int_path))
  stop("ST7 outlier/dose cache missing: ", st7_int_path,
       ". Run analysis/01_main_analysis/16_ST_vaf_outlier_sensitivity.R first.",
       call. = FALSE)
st7_cache <- readRDS(st7_int_path)

.fmt_p_st7 <- function(p) {
  ifelse(p < 0.001,
         formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

# Canonical phenotype labels (match ST4 / ST6 harmonization).
.canon_outcome_st7 <- function(x) dplyr::recode(x,
  "SM total score"      = "Spetzler\u2013Martin total grade",
  "SM size"             = "Spetzler\u2013Martin size",
  "Deep drainage"       = "Deep venous drainage",
  "Eloquence"           = "Eloquent brain location",
  "High-risk count"     = "High-risk feature count",
  "Rupture"             = "Rupture (ever)",
  "Age at presentation" = "Age at presentation",
  "Age"                 = "Age at presentation"
)

# \u2500\u2500 Block 1: outlier sensitivity (Spearman \u03c1, full vs |z|\u22642 trimmed) \u2500\u2500
block_outlier <- bind_rows(
  st7_cache$sens_rows %>%
    transmute(
      analysis  = "Outlier (|z|>2) sensitivity",
      scenario  = "Full cohort",
      outcome   = .canon_outcome_st7(outcome),
      n         = as.integer(n_full),
      statistic = "Spearman \u03c1",
      estimate  = sprintf("%.3f", rho_full),
      ci_str    = sprintf("(%.3f, %.3f)", lo_full, hi_full),
      p_str     = .fmt_p_st7(p_full)
    ),
  st7_cache$sens_rows %>%
    transmute(
      analysis  = "Outlier (|z|>2) sensitivity",
      scenario  = "|z|\u22642 trimmed",
      outcome   = .canon_outcome_st7(outcome),
      n         = as.integer(n_trim),
      statistic = "Spearman \u03c1",
      estimate  = sprintf("%.3f", rho_trim),
      ci_str    = sprintf("(%.3f, %.3f)", lo_trim, hi_trim),
      p_str     = .fmt_p_st7(p_trim)
    )
) %>% arrange(outcome, scenario)

# \u2500\u2500 Block 2: dose-response models (Model 1/2/3 coefficients per outcome) \u2500\u2500
# dose_results already carries a `type` column (continuous / binary) so no
# join with dose_outcomes is needed.
block_dose <- st7_cache$dose_results %>%
  transmute(
    analysis  = "Dose-response model",
    scenario  = paste0("Pooled ", predictor),
    outcome   = .canon_outcome_st7(outcome),
    n         = as.integer(n),
    statistic = ifelse(type == "binary", "Logistic \u03b2", "OLS \u03b2"),
    estimate  = sprintf("%.3f", estimate),
    ci_str    = sprintf("(%.3f, %.3f)", conf.low, conf.high),
    p_str     = .fmt_p_st7(p.value)
  )

# \u2500\u2500 Block 3: LRT sensitivity (Model 2 vs Model 1, full vs trimmed) \u2500\u2500
# Pure LRT here \u2014 Estimate/CI are sentinel; the Model 2 VAF \u03b2 with CI is
# already captured in Block 2 under "VAF/1% (Model 2)".
block_lrt <- bind_rows(
  st7_cache$lrt_compare %>%
    transmute(
      analysis  = "Dose-response LRT sensitivity",
      scenario  = "Full cohort",
      outcome   = .canon_outcome_st7(outcome),
      n         = as.integer(n_full),
      statistic = "LRT P (Model 2 vs Model 1)",
      estimate  = "\u2014",
      ci_str    = "\u2014",
      p_str     = .fmt_p_st7(lrt_p_full)
    ),
  st7_cache$lrt_compare %>%
    transmute(
      analysis  = "Dose-response LRT sensitivity",
      scenario  = "|z|\u22642 trimmed",
      outcome   = .canon_outcome_st7(outcome),
      n         = as.integer(n_trim),
      statistic = "LRT P (Model 2 vs Model 1)",
      estimate  = "\u2014",
      ci_str    = "\u2014",
      p_str     = .fmt_p_st7(lrt_p_trim)
    )
) %>% arrange(outcome, scenario)

# \u2500\u2500 Block 4: slope sensitivity LOSO + VAF trim + Cohort restriction \u2500\u2500
.slope_block <- function(d, analysis_label) {
  d %>%
    transmute(
      analysis  = analysis_label,
      scenario  = scenario,
      outcome   = "Age at presentation",
      n         = as.integer(n),
      statistic = "OLS \u03b2 (VAF, y per 1%)",
      estimate  = sprintf("%.3f", slope),
      ci_str    = sprintf("(%.3f, %.3f)", ci_low, ci_high),
      p_str     = .fmt_p_st7(p_value)
    )
}
block_loso_df   <- .slope_block(st10_loso,   "LOSO")
block_trim_df   <- .slope_block(st10_trim,   "VAF trim")
block_cohort_df <- .slope_block(st10_cohort, "Cohort restriction")

# Distribution metadata \u2192 compact footnote shared on sheet (A).
.dist <- setNames(st7_cache$vaf_distribution$value,
                  st7_cache$vaf_distribution$metric)
.dist_footnote <- sprintf(
  "Underlying VAF distribution (variant-positive lesions with measured VAF, N = %d): mean = %.2f%%, SD = %.2f, median = %.2f (IQR %.2f\u2013%.2f), range %.2f\u2013%.2f. Skewness = %.2f, excess kurtosis = %.2f. Shapiro\u2013Wilk on raw VAF: W = %.3f, P = %s; on log(VAF): W = %.3f, P = %s. |z|>2: %d lesions (%.1f%%); |z|>3: %d lesions (%.1f%%).",
  .dist[["n"]], .dist[["mean_pct"]], .dist[["sd_pct"]],
  .dist[["median_pct"]], .dist[["IQR_low_pct"]], .dist[["IQR_high_pct"]],
  .dist[["min"]], .dist[["max"]],
  .dist[["skewness"]], .dist[["excess_kurtosis"]],
  .dist[["shapiro_W_raw"]], .fmt_p_st7(.dist[["shapiro_p_raw"]]),
  .dist[["shapiro_W_log"]], .fmt_p_st7(.dist[["shapiro_p_log"]]),
  as.integer(.dist[["n_|z|>2"]]), .dist[["pct_|z|>2"]],
  as.integer(.dist[["n_|z|>3"]]), .dist[["pct_|z|>3"]]
)

st7_path <- file.path(etable_dir, "SuppTable07_vaf_age_sensitivity.xlsx")
if (file.exists(st7_path)) unlink(st7_path)

# Per-analysis sheet columns: each sheet has its own natural schema.
.outlier_cols <- list(
  col    ("scenario",  label = "Scenario"),
  col    ("outcome",   label = "Outcome"),
  col_int("n",         label = "N",        italic = TRUE),
  col    ("estimate",  label = "\u03c1",   italic = TRUE),
  col    ("ci_str",    label = "95% CI"),
  col    ("p_str",     label = "P",        italic = TRUE)
)
.dose_cols <- list(
  col    ("scenario",  label = "Predictor"),
  col    ("outcome",   label = "Outcome"),
  col_int("n",         label = "N",        italic = TRUE),
  col    ("statistic", label = "Model"),
  col    ("estimate",  label = "\u03b2",   italic = TRUE),
  col    ("ci_str",    label = "95% CI"),
  col    ("p_str",     label = "P",        italic = TRUE)
)
.lrt_cols <- list(
  col    ("scenario",  label = "Scenario"),
  col    ("outcome",   label = "Outcome"),
  col_int("n",         label = "N",        italic = TRUE),
  col    ("p_str",     label = "LRT P (Model 2 vs Model 1)", italic = TRUE)
)
.slope_cols <- list(
  col    ("scenario",  label = "Scenario"),
  col_int("n",         label = "N",        italic = TRUE),
  col    ("estimate",  label = "\u03b2 (y per 1% VAF)", italic = TRUE),
  col    ("ci_str",    label = "95% CI"),
  col    ("p_str",     label = "P",        italic = TRUE)
)

# Drop columns no longer needed under per-sheet schemas.
.strip <- function(d, keep) d[, intersect(keep, names(d)), drop = FALSE]

st7_sheets <- list(
  "(A) Outlier sensitivity" = list(
    data = .strip(block_outlier, c("scenario","outcome","n","estimate","ci_str","p_str")),
    columns = .outlier_cols,
    footnote = c(
      "Per-outcome Spearman \u03c1 for VAF \u00d7 phenotype, computed on the full variant-positive cohort and on the |z|\u22642 trimmed subset. Bootstrap 95% CIs (1000 resamples).",
      .dist_footnote
    )
  ),
  "(B) Dose-response model" = list(
    data = .strip(block_dose, c("scenario","outcome","n","statistic","estimate","ci_str","p_str")),
    columns = .dose_cols,
    footnote = c(
      "Nested OLS / logistic fits for each outcome. Estimate is the per-coefficient \u03b2 (continuous outcomes) or log-odds \u03b2 (binary).",
      "Model 1 = mut+ only; Model 2 = mut+ + within-mut VAF (centered); Model 3 = VAF only."
    )
  ),
  "(C) LRT sensitivity" = list(
    data = .strip(block_lrt, c("scenario","outcome","n","p_str")),
    columns = .lrt_cols,
    footnote = c(
      "Likelihood-ratio test, Model 2 (mut+ + within-mut VAF) vs Model 1 (mut+ only), per outcome, under the full cohort and the |z|\u22642 trimmed subset.",
      "The Estimate and 95% CI for the Model 2 VAF \u03b2 are reported in sheet (B) under 'VAF/1% (Model 2)'."
    )
  ),
  "(D) Slope - LOSO" = list(
    data = .strip(block_loso_df, c("scenario","n","estimate","ci_str","p_str")),
    columns = .slope_cols,
    footnote = "OLS \u03b2 (years per 1% VAF) for age-at-presentation under leave-one-series-out perturbation of the pooled (KRAS G12D + G12V) variant-positive cohort."
  ),
  "(E) Slope - VAF trim" = list(
    data = .strip(block_trim_df, c("scenario","n","estimate","ci_str","p_str")),
    columns = .slope_cols,
    footnote = "OLS \u03b2 (years per 1% VAF) for age-at-presentation under VAF percentile-trim perturbation of the pooled variant-positive cohort."
  ),
  "(F) Slope - Cohort restrict" = list(
    data = .strip(block_cohort_df, c("scenario","n","estimate","ci_str","p_str")),
    columns = .slope_cols,
    footnote = "OLS \u03b2 (years per 1% VAF) for age-at-presentation under cohort-restriction perturbation (institutional vs published)."
  )
)

write_supp_table_workbook(
  sheets = st7_sheets,
  path   = st7_path
)
cat(sprintf("  \u2713 SuppTable07_vaf_age_sensitivity.xlsx (6 sheets: %s rows)\n",
            paste(vapply(st7_sheets, function(x) nrow(x$data), integer(1)),
                  collapse = ", ")))

cat("\n══ 20_ST_supp_tables_8_9_10.R complete ══\n")
