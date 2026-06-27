# 18_ED10_age_adj_vaf_pheno.R — Theme F2 (Hale v2)
# ─────────────────────────────────────────────────────────────────────────────
# Andy's concern (Hale v2, item F2): *"Let's adjust for age and see how VAF
# predicts rupture, size, etc."* VAF and age are codependent in the Mut+
# cohort (negative slope in Fig 4D), so an unadjusted `y ~ VAF` coefficient
# could be carrying hidden age variance. This script refits each primary
# phenotype on the Mut+ subset with and without age adjustment and reports
# the VAF coefficient + 95 % CI side-by-side.
#
# Design:
#   Restrict to variant-positive cases with non-missing VAF and age
#     (Mut+ anchor is necessary because VAF = 0 for Negatives is structural,
#      not a dose; mixing them would collapse age-adjustment into the
#      Mut+-vs-Neg contrast).
#   VAF on percent scale (1% unit interpretation, matches Fig 4D).
#   For each outcome fit two models:
#     Model U: y ~ vaf_pct                   (unadjusted)
#     Model A: y ~ age + vaf_pct             (age-adjusted)
#   Continuous outcomes → Gaussian `lm`.
#   Binary outcomes → Firth-penalised logistic (`logistf`) to avoid
#     separation in small subsets with rare positive cells.
#   Also compute VIF for Model A as a collinearity diagnostic (Theme F4).
#
# Outputs:
#   - results/ExtendedData/ed_age_adj_vaf_pheno/age_adj_correlation.{png,pdf}
#     (forest: two rows per outcome, unadjusted vs age-adjusted VAF coef)
#   - results/SupplementaryTables/SuppTable06_vaf_phenotype_correlations.xlsx
#     (appended sheets: `Age-adjusted coefficients`, `Age-adjusted VIF`,
#      `Age-adjusted cohort summary` — folded in via write_supp_table_sheets()
#      alongside the per-variant correlation sheets produced by
#      20_ST_supp_tables_8_9_10.R)
#   - results/stats/vaf_age_adjusted.txt (human-readable dump)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  # NB: car not attached because car::recode masks dplyr::recode for every
  # downstream script in the run_all loop. VIF is computed manually below.
  library(ggplot2); library(here); library(logistf)
  library(openxlsx); library(broom)
})
source(here("analysis", "helper_scripts", "utils.R"))

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
save_panel <- function(dir, name, plot, w, h)
  save_panel_impl(dir, name, plot, w, h, device = "cairo", save_rds = FALSE)

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))

# ── Analysis cohort: Mut+ with VAF and age ───────────────────────────────────
mp <- df %>%
  filter(mutation_positive, !is.na(vaf_prop), !is.na(age)) %>%
  mutate(vaf_pct = vaf_prop * 100)

cat(sprintf("── F2: VAF-phenotype models (age-adjusted vs unadjusted) ──\n"))
cat(sprintf("  Analysis cohort: mut+ with VAF & age, n = %d\n", nrow(mp)))
cat(sprintf("  Mean age = %.1f y, Mean VAF = %.2f %%\n",
            mean(mp$age), mean(mp$vaf_pct)))
cat(sprintf("  Cor(age, VAF) = %.3f (Spearman)\n",
            cor(mp$age, mp$vaf_pct, method = "spearman")))

# ── Outcomes to model ───────────────────────────────────────────────────────
outcomes <- tribble(
  ~var,                           ~label,                      ~type,
  "age",                          "Age (self-check)",           "continuous",
  "sm_size_num",                  "SM size",                    "continuous",
  "sm_total_num",                 "SM total score",             "continuous",
  "n_high_risk_num",              "High-risk feature count",    "continuous",
  "ever_ruptured_num",            "Rupture (ever)",             "binary",
  "sm_drainage_num",              "Deep venous drainage",       "binary",
  "sm_eloquence_num",             "Eloquent brain location",    "binary",
  "intranidal_aneurysm_num",      "Intranidal aneurysm",        "binary",
  "venous_varix_num",             "Venous varix",               "binary",
  "venous_outflow_stenosis_num",  "Venous outflow stenosis",    "binary",
  "flow_related_aneurysm_num",    "Flow-related aneurysm",      "binary"
)

# ── Fit helper ───────────────────────────────────────────────────────────────
# For binary outcomes use Firth; for continuous use Gaussian lm. In both
# cases we return β_VAF (per +1 % VAF) with 95 % CI and p, plus the VIF
# for the adjusted model.
fit_pair <- function(y_var, type, data) {
  d <- data %>% filter(!is.na(.data[[y_var]])) %>% as.data.frame()
  n <- nrow(d)
  if (n < 15) return(NULL)

  fml_u <- as.formula(paste(y_var, "~ vaf_pct"))
  fml_a <- as.formula(paste(y_var, "~ age + vaf_pct"))

  get_vaf <- function(mod, term = "vaf_pct") {
    if (inherits(mod, "logistf")) {
      i <- match(term, names(coef(mod)))
      est <- as.numeric(coef(mod)[i])
      lo  <- as.numeric(mod$ci.lower[i])
      hi  <- as.numeric(mod$ci.upper[i])
      p   <- as.numeric(mod$prob[i])
    } else {
      ci <- confint(mod, term, level = 0.95)
      cf <- summary(mod)$coefficients[term, ]
      est <- cf["Estimate"]; lo <- ci[1]; hi <- ci[2]
      p   <- cf["Pr(>|t|)"]
    }
    tibble(beta = est, lo = lo, hi = hi, p = p)
  }

  if (type == "binary") {
    # Firth for both models; use pl = FALSE to keep Wald-style CIs/p
    # comparable across models (Firth profile CIs can be very wide in
    # sparse cells and make the side-by-side unreadable).
    m_u <- logistf(fml_u, data = d, pl = FALSE)
    m_a <- logistf(fml_a, data = d, pl = FALSE)
  } else {
    m_u <- lm(fml_u, data = d)
    m_a <- lm(fml_a, data = d)
  }

  # VIF on the adjusted model — computed manually because it depends only
  # on the predictor design matrix (same for binary/continuous outcomes
  # once we restrict to the outcome's non-missing rows). VIF(vaf_pct) =
  # 1 / (1 - R^2) from the auxiliary regression vaf_pct ~ age.
  vif_val <- tryCatch({
    r2 <- summary(lm(vaf_pct ~ age, data = d))$r.squared
    1 / (1 - r2)
  }, error = function(e) NA_real_)

  bind_rows(
    get_vaf(m_u) %>% mutate(model = "Unadjusted"),
    get_vaf(m_a) %>% mutate(model = "Age-adjusted")
  ) %>%
    mutate(n = n, vif_vaf_in_adjusted = vif_val, type = type)
}

coef_tbl <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  res <- fit_pair(outcomes$var[i], outcomes$type[i], mp)
  if (is.null(res)) return(NULL)
  res %>% mutate(outcome = outcomes$label[i],
                 variable = outcomes$var[i]) %>%
    select(outcome, variable, type, model, n, beta, lo, hi, p,
           vif_vaf_in_adjusted)
})

cat("\n── Coefficients table ──\n")
print(coef_tbl, width = Inf, n = Inf)

# ── Extended Data Fig. 21: side-by-side forest (unadjusted vs age-adjusted β) ──────────
# For continuous outcomes β is in raw-scale units per +1 % VAF. For binary
# outcomes β is log-OR per +1 % VAF. We facet by type because the
# magnitudes are on different scales, and colour-code the two models.

pdat <- coef_tbl %>%
  filter(outcome != "Age (self-check)") %>%   # self-check printed to log only
  mutate(
    # Keep the declared outcome order but reverse (top-to-bottom in plot)
    outcome = factor(outcome,
                     levels = rev(outcomes$label[outcomes$label != "Age (self-check)"])),
    model   = factor(model, levels = c("Unadjusted", "Age-adjusted"))
  )

# v6.45 (2026-05-21): split into two separate panels (binary vs continuous)
# so each gets its own 3.30-in cell in the ED10 2\u00d72 grid. Clean titles \u2014
# just "Binary outcomes" / "Continuous outcomes"; the \u03b2-definition gloss
# moves to the figure caption.
build_age_adj_panel <- function(sub_df, title_text) {
  ggplot(sub_df, aes(x = beta, y = outcome, color = model, shape = model)) +
    ref_vline(0, kind = "null") +
    geom_linerange(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = 0.55), linewidth = 0.7) +
    geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
    scale_color_manual(values = c("Unadjusted"   = "#1A1A1A",
                                  "Age-adjusted" = unname(PAL_SCORE[["1"]]))) +
    scale_shape_manual(values = c("Unadjusted" = 15, "Age-adjusted" = 17)) +
    # v6.48 (2026-05-21): add x-axis label so the units are explicit on
    # each panel (β = log-OR for binary outcomes, β = outcome units for
    # continuous). Per-figure caption carries the full gloss.
    labs(title = title_text,
         x = "β (per +1% VAF)",
         y = NULL,
         color = NULL, shape = NULL) +
    theme_nature_panel() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position  = "top",
          plot.title       = element_text(face = "bold",
                                          size = NM$label_pt,
                                          family = NM$font_family))
}

p_binary     <- build_age_adj_panel(
  pdat %>% filter(type == "binary"),
  title_text = "Binary outcomes"
)
p_continuous <- build_age_adj_panel(
  pdat %>% filter(type == "continuous"),
  title_text = "Continuous outcomes"
)

# Registry-resolved paths \u2014 split into two tokens so each lands in its
# own panel_X/ slot.
out_dir_bin <- panel_slot_dir("age_adj_binary")
out_dir_con <- panel_slot_dir("age_adj_continuous")
dir.create(out_dir_bin, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_con, recursive = TRUE, showWarnings = FALSE)
save_panel(out_dir_bin, "age_adj_binary",     p_binary,     w = 3.30, h = 3.00)
save_panel(out_dir_con, "age_adj_continuous", p_continuous, w = 3.30, h = 3.00)
saveRDS(p_binary,     file.path(out_dir_bin, "age_adj_binary.rds"))
saveRDS(p_continuous, file.path(out_dir_con, "age_adj_continuous.rds"))

# ── Supplementary Table 06 — single-sheet mega-table ───────────────────────
# 2026-06-14: collapsed from 8 sheets (4 per-cohort correlation +
# Age-adjusted coefficients + Age-adjusted VIF + Age-adjusted cohort summary
# + a duplicate stale Pooled sheet) into one stacked sheet. The correlation
# block is computed in 20_ST_supp_tables_8_9_10.R and cached at
# results/stats/_intermediates/st6_correlations.rds — read it here, append
# the age-adjusted block, write the unified file.
etab_dir <- here("results", "SupplementaryTables")
dir.create(etab_dir, recursive = TRUE, showWarnings = FALSE)
source(here("analysis", "helper_scripts", "supp_table_writer.R"))

st6_corr_path <- here("results", "stats", "_intermediates",
                      "st6_correlations.rds")
if (!file.exists(st6_corr_path))
  stop("ST6 correlation cache missing: ", st6_corr_path,
       ". Run analysis/01_main_analysis/20_ST_supp_tables_8_9_10.R first.",
       call. = FALSE)
st6_corr <- readRDS(st6_corr_path)

# Harmonize phenotype labels: correlation rows use slightly different forms
# than age-adjusted rows ("Spetzler–Martin total grade" vs "SM total score",
# "Eloquence" vs "Eloquent brain location", "Rupture status (ever)" vs
# "Rupture (ever)"). Pick one canonical set for the mega-sheet so the
# Phenotype column reads consistently across blocks.
.canon_pheno <- function(x) dplyr::recode(x,
  "Spetzler–Martin total grade" = "Spetzler–Martin total grade",
  "SM total score"                    = "Spetzler–Martin total grade",
  "Spetzler–Martin size"          = "Spetzler–Martin size",
  "SM size"                           = "Spetzler–Martin size",
  "Eloquence"                         = "Eloquent brain location",
  "Eloquent brain location"           = "Eloquent brain location",
  "Rupture status (ever)"             = "Rupture (ever)",
  "Rupture (ever)"                    = "Rupture (ever)",
  "Deep venous drainage"              = "Deep venous drainage",
  "High-risk feature count"           = "High-risk feature count",
  "Age (self-check)"                  = "Age (self-check)",
  "Intranidal aneurysm"               = "Intranidal aneurysm",
  "Venous varix"                      = "Venous varix",
  "Venous outflow stenosis"           = "Venous outflow stenosis",
  "Flow-related aneurysm"             = "Flow-related aneurysm"
)

.fmt_p <- function(p) {
  ifelse(p < 0.001,
         formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

# 2026-06-15: ST6 back to 5 sheets (Nature MOESM3 (A)-(E) pattern).
# (A)-(D) per-cohort correlations (Spearman ρ + logistic β); (E) age-adjusted
# regressions on the pooled variant-positive cohort.
.corr_sheet <- function(d) {
  d %>%
    transmute(
      phenotype  = .canon_pheno(phenotype),
      n          = as.integer(n),
      method     = ifelse(statistic == "Spearman ρ",
                          "Spearman ρ",
                          "Logistic β (VAF, proportion)"),
      estimate   = sprintf("%.3f", estimate),
      p_str      = .fmt_p(p_value),
      fdr_p_str  = .fmt_p(p_fdr)
    )
}
.corr_cols <- list(
  col    ("phenotype", label = "Phenotype"),
  col_int("n",         label = "N",     italic = TRUE),
  col    ("method",    label = "Method"),
  col    ("estimate",  label = "Estimate"),
  col    ("p_str",     label = "P",     italic = TRUE),
  col    ("fdr_p_str", label = "FDR P", italic = TRUE)
)
.corr_footnote <- "Continuous outcomes use Spearman ρ; binary outcomes use logistic regression with VAF (as proportion) as the predictor. FDR P uses BH correction within this sheet."

# Age-adjusted block: 11 outcomes × 2 models (unadjusted + age-adjusted).
adj_sheet_df <- coef_tbl %>%
  transmute(
    phenotype  = .canon_pheno(outcome),
    n          = as.integer(n),
    method     = ifelse(type == "continuous",
                        ifelse(model == "Unadjusted",
                               "OLS β (unadjusted)",
                               "OLS β (age-adjusted)"),
                        ifelse(model == "Unadjusted",
                               "Firth logistic β (unadjusted)",
                               "Firth logistic β (age-adjusted)")),
    estimate   = sprintf("%.3f", beta),
    ci_str     = sprintf("(%.3f, %.3f)", lo, hi),
    p_str      = .fmt_p(p),
    vif_str    = ifelse(model == "Age-adjusted",
                        sprintf("%.2f", vif_vaf_in_adjusted),
                        "—")
  )

.cohort_footnote <- sprintf(
  "Age-adjusted analysis cohort: variant-positive bAVMs with measured VAF AND measured age (N = %d). Mean age = %.1f y, mean VAF = %.2f %%. Spearman r(age, VAF) = %.3f; Pearson r = %.3f.",
  nrow(mp), mean(mp$age), mean(mp$vaf_pct),
  cor(mp$age, mp$vaf_pct, method = "spearman"),
  cor(mp$age, mp$vaf_pct, method = "pearson")
)

st6_path <- file.path(etab_dir, "SuppTable06_vaf_phenotype_correlations.xlsx")
if (file.exists(st6_path)) unlink(st6_path)

st6_sheets <- list(
  "(A) Pooled variant-positive" = list(
    data = .corr_sheet(st6_corr$pooled_variant_positive),
    columns = .corr_cols,
    footnote = .corr_footnote
  ),
  "(B) KRAS G12D" = list(
    data = .corr_sheet(st6_corr$kras_g12d),
    columns = .corr_cols,
    footnote = .corr_footnote
  ),
  "(C) KRAS G12V" = list(
    data = .corr_sheet(st6_corr$kras_g12v),
    columns = .corr_cols,
    footnote = .corr_footnote
  ),
  "(D) Other KRAS + BRAF" = list(
    data = .corr_sheet(st6_corr$other_kras_braf),
    columns = .corr_cols,
    footnote = .corr_footnote
  ),
  "(E) Age-adjusted (pooled)" = list(
    data = adj_sheet_df,
    columns = list(
      col    ("phenotype", label = "Phenotype"),
      col_int("n",         label = "N",         italic = TRUE),
      col    ("method",    label = "Method"),
      col    ("estimate",  label = "β",    italic = TRUE),
      col    ("ci_str",    label = "95% CI"),
      col    ("p_str",     label = "P",         italic = TRUE),
      col    ("vif_str",   label = "VIF (VAF)", italic = TRUE)
    ),
    footnote = c(
      "Paired Unadjusted vs Age-adjusted regression of each outcome on VAF, restricted to variant-positive lesions with both VAF and age measured. OLS β for continuous outcomes; Firth-penalized logistic β for binary.",
      "VIF (VAF) is the variance-inflation factor for VAF in the age-adjusted model (collinearity diagnostic; values >5 would flag concern). '—' = not applicable (only the age-adjusted row carries a VIF).",
      .cohort_footnote
    )
  )
)

write_supp_table_workbook(
  sheets = st6_sheets,
  path   = st6_path
)
cat(sprintf("  ✓ SuppTable06 (5 sheets: %s rows)\n",
            paste(vapply(st6_sheets, function(x) nrow(x$data), integer(1)),
                  collapse = ", ")))
st6_df <- bind_rows(lapply(st6_sheets, function(x) x$data))

# ── Human-readable stats dump ────────────────────────────────────────────────
stats_dir <- here("results", "stats")
writeLines(c(
  sprintf("# %s — age-adjusted VAF-phenotype models",
          panel_prose_tag(c("age_adj_binary", "age_adj_continuous"))),
  sprintf("Generated %s", Sys.Date()),
  sprintf("Cohort: mut+ with VAF & age, n = %d", nrow(mp)),
  sprintf("Spearman r(age, VAF) = %.3f",
          cor(mp$age, mp$vaf_pct, method = "spearman")),
  "",
  "## Coefficients",
  capture.output(print(coef_tbl, width = Inf, n = Inf)),
  "",
  "## VIF (adjusted model \u2014 VAF column from coef_tbl, age-adjusted rows)",
  capture.output(print(
    coef_tbl %>%
      filter(model == "Age-adjusted") %>%
      distinct(outcome, variable, type, n, vif_vaf_in_adjusted),
    width = Inf))),
  file.path(stats_dir, "vaf_age_adjusted.txt"))

cat("\n\u2713 Extended Data Fig. 21 figure + SuppTable06 saved.\n")
cat(sprintf("   %s (single sheet: VAF \u00d7 phenotype, %d rows)\n",
            file.path(etab_dir, "SuppTable06_vaf_phenotype_correlations.xlsx"),
            nrow(st6_df)))
