# 16_ST_vaf_outlier_sensitivity.R — VAF distribution, outlier sensitivity
# (Clogg coefficient comparison), and dose-response models.
#
# Supplementary-table producer ONLY — no figure is rendered. The earlier 2x2 ED
# panel (vaf_outlier_combined) was retired from the manuscript; this script now
# computes the numeric frames behind Supplementary Table 7 and the VAF-outlier
# stats fragment.
#
# Output:
#   results/stats/_intermediates/st7_outlier_dose.rds — frames assembled into
#     Supplementary Table 7 by 20_ST_supp_tables_8_9_10.R
#   stats fragment `edfig_vaf_outlier` — VAF distribution / Shapiro / outlier %
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(openxlsx)
  library(broom)
})

source(here("analysis", "helper_scripts", "utils.R"))

set.seed(MASTER_SEED)   # audit F13: canonical seed from utils.R

# ── Load data ────────────────────────────────────────────────────────────────
df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))

vaf_df <- df %>%
  filter(mutation_positive == TRUE & !is.na(vaf_prop)) %>%
  mutate(vaf_pct = vaf_prop * 100)

v  <- vaf_df$vaf_pct
m  <- mean(v)
s  <- sd(v)
cat(sprintf("VAF cohort: n=%d | mean=%.2f | sd=%.2f | median=%.2f\n",
            length(v), m, s, median(v)))

# Outlier mask
vaf_df$is_outlier <- abs((vaf_df$vaf_pct - m) / s) > 2
cat(sprintf("  |z|>2 outliers: %d (%.1f%%)\n",
            sum(vaf_df$is_outlier), 100 * mean(vaf_df$is_outlier)))


# ═════════════════════════════════════════════════════════════════════════════
# Outlier sensitivity — full vs trimmed, with bootstrap CIs and
# Clogg-style z-test for coefficient difference
# ═════════════════════════════════════════════════════════════════════════════
cat("── Outlier sensitivity ──\n")

# Fisher-z transform a correlation; SE = 1/sqrt(n-3).
# Clogg-style z for difference between two independent correlations is
# approximate here (trimmed sample is a subset of full), so we additionally
# bootstrap each rho to report honest 95% CIs.
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
inv_fisher_z <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)

boot_spearman <- function(x, y, n_boot = 2000) {
  n <- length(x)
  if (n < 8) return(c(NA, NA))
  rs <- replicate(n_boot, {
    i <- sample(seq_len(n), replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  stats::quantile(rs, c(0.025, 0.975), na.rm = TRUE)
}

outcomes <- tribble(
  ~var,               ~label,
  "sm_total_num",     "SM total score",
  "sm_size_num",      "SM size",
  "sm_drainage_num",  "Deep drainage",
  "sm_eloquence_num", "Eloquence",
  "n_high_risk_num",  "High-risk count",
  "ever_ruptured_num","Rupture",
  "age",              "Age at presentation"
)

sens_rows <- map_dfr(seq_len(nrow(outcomes)), function(i) {
  vr <- outcomes$var[i]; lb <- outcomes$label[i]
  ok   <- !is.na(vaf_df$vaf_pct) & !is.na(vaf_df[[vr]])
  ok_t <- ok & !vaf_df$is_outlier

  full <- suppressWarnings(cor.test(vaf_df$vaf_pct[ok],   vaf_df[[vr]][ok],   method = "spearman"))
  trim <- suppressWarnings(cor.test(vaf_df$vaf_pct[ok_t], vaf_df[[vr]][ok_t], method = "spearman"))

  bf <- boot_spearman(vaf_df$vaf_pct[ok],   vaf_df[[vr]][ok])
  bt <- boot_spearman(vaf_df$vaf_pct[ok_t], vaf_df[[vr]][ok_t])

  # Clogg-style z (Fisher-z on each rho); approximate since samples overlap
  z1 <- fisher_z(full$estimate); n1 <- sum(ok)
  z2 <- fisher_z(trim$estimate); n2 <- sum(ok_t)
  clogg_z <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  clogg_p <- 2 * pnorm(-abs(clogg_z))

  tibble(
    outcome = lb,
    rho_full = full$estimate, lo_full = bf[1], hi_full = bf[2],
    p_full   = full$p.value,  n_full  = n1,
    rho_trim = trim$estimate, lo_trim = bt[1], hi_trim = bt[2],
    p_trim   = trim$p.value,  n_trim  = n2,
    clogg_z  = clogg_z,
    clogg_p  = clogg_p
  )
})

print(sens_rows)

# ═════════════════════════════════════════════════════════════════════════════
# Dose-response models (Mut only / Mut+VAF / VAF only)
# ═════════════════════════════════════════════════════════════════════════════
cat("── Dose-response models ──\n")

# Decorrelated parameterisation:
#   - mut01                 : 1 if mutation present, else 0 → effect of having
#     any mutation at the mean within-mutant VAF
#   - vaf_pct_centered_in_mut : (VAF% - mean VAF% among mut+) if mutation
#     present, else 0 → within-mutant dose effect; orthogonal to mut01 (VIF ~ 1).
# This is scale-interpretable: β is "change in outcome per +1% VAF".
# The prior parameterisation (raw 0–1 proportion, uncentered) yielded
# β ≈ −187 for age simply because a 1.0-unit change on the 0–1 scale is
# biologically impossible; percent units give β ≈ −1.87 years per +1% VAF.
mean_vaf_mut_pct <- mean(df$vaf_prop[df$mutation_positive & !is.na(df$vaf_prop)]) * 100
cat(sprintf("Mean VAF among mut+ cases: %.2f%% (used to centre within-mutant VAF)\n",
            mean_vaf_mut_pct))

dr_df <- df %>%
  mutate(mut01    = as.integer(mutation_positive),
         vaf_pct_dose              = ifelse(mutation_positive, vaf_prop * 100, 0),
         vaf_pct_centered_in_mut   = ifelse(mutation_positive,
                                             vaf_prop * 100 - mean_vaf_mut_pct, 0))

dose_outcomes <- tribble(
  ~var,               ~label,             ~type,
  "age",              "Age",              "continuous",
  "sm_size_num",      "SM size",          "continuous",
  "n_high_risk_num",  "High-risk count",  "continuous",
  "ever_ruptured_num","Rupture",          "binary",
  "sm_drainage_num",  "Deep drainage",    "binary",
  "sm_eloquence_num", "Eloquence",        "binary"
)

fit_pair <- function(y_var, type, data) {
  d <- data %>% filter(!is.na(.data[[y_var]]) & !is.na(vaf_pct_dose) & !is.na(mut01))
  if (nrow(d) < 15) return(NULL)
  fml <- function(rhs) as.formula(paste(y_var, "~", rhs))
  if (type == "binary") {
    m1 <- glm(fml("mut01"),                             data = d, family = binomial)
    m2 <- glm(fml("mut01 + vaf_pct_centered_in_mut"),   data = d, family = binomial)
    m3 <- glm(fml("vaf_pct_dose"),                      data = d, family = binomial)
  } else {
    m1 <- glm(fml("mut01"),                             data = d)
    m2 <- glm(fml("mut01 + vaf_pct_centered_in_mut"),   data = d)
    m3 <- glm(fml("vaf_pct_dose"),                      data = d)
  }
  lrt <- anova(m1, m2, test = "LRT")
  list(m1 = m1, m2 = m2, m3 = m3, lrt = lrt, n = nrow(d), type = type)
}

dose_results <- map_dfr(seq_len(nrow(dose_outcomes)), function(i) {
  vr <- dose_outcomes$var[i]; lb <- dose_outcomes$label[i]; tp <- dose_outcomes$type[i]
  res <- fit_pair(vr, tp, dr_df)
  if (is.null(res)) return(NULL)

  get_term <- function(mod, term_name) {
    cf <- tidy(mod, conf.int = TRUE) %>% filter(term == term_name)
    if (nrow(cf) == 0) tibble(estimate = NA, conf.low = NA, conf.high = NA, p.value = NA)
    else cf %>% select(estimate, conf.low, conf.high, p.value)
  }

  bind_rows(
    tibble(outcome = lb, predictor = "Mut+ (Model 1)",      n = res$n) %>% bind_cols(get_term(res$m1, "mut01")),
    tibble(outcome = lb, predictor = "Mut+ (Model 2)",      n = res$n) %>% bind_cols(get_term(res$m2, "mut01")),
    tibble(outcome = lb, predictor = "VAF/1% (Model 2)",    n = res$n) %>% bind_cols(get_term(res$m2, "vaf_pct_centered_in_mut")),
    tibble(outcome = lb, predictor = "VAF/1% (Model 3)",    n = res$n) %>% bind_cols(get_term(res$m3, "vaf_pct_dose"))
  ) %>% mutate(
    lrt_p = res$lrt$`Pr(>Chi)`[2],
    type  = res$type
  )
})

print(dose_results)

# ─────────────────────────────────────────────────────────────────────────────
# LRT sensitivity: refit dose-response on |z| <= 2 trimmed cohort.
#
# We already trim the SPEARMAN-rho analysis (panel C); this mirrors the same
# trim on the LRT framework so readers can see whether "VAF adds information
# beyond binary mutation status" (Model 2 vs Model 1) survives removal of the
# right-tail BCH cases. Trim is applied on the variant-positive VAF
# distribution — mutation-negative rows (all VAF = 0) are always kept because
# they carry the "mut+ vs neg" contrast (Model 1 effect) and shouldn't be
# flagged as outliers on a zero-valued column.
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Panel D sensitivity: LRT on |z| <= 2 trimmed cohort ──\n")

mut_vaf <- df$vaf_prop[df$mutation_positive & !is.na(df$vaf_prop)] * 100
trim_mean <- mean(mut_vaf);  trim_sd <- sd(mut_vaf)
dr_trim <- dr_df %>%
  filter(!mutation_positive |
         (is.finite(vaf_pct_dose) &
          abs((vaf_pct_dose - trim_mean) / trim_sd) <= 2))
cat(sprintf("  Full n (all rows, any VAF) = %d;  trimmed n = %d  (dropped %d mut+ |z|>2)\n",
            nrow(dr_df), nrow(dr_trim), nrow(dr_df) - nrow(dr_trim)))

lrt_compare <- map_dfr(seq_len(nrow(dose_outcomes)), function(i) {
  vr <- dose_outcomes$var[i]; lb <- dose_outcomes$label[i]; tp <- dose_outcomes$type[i]
  full <- fit_pair(vr, tp, dr_df)
  trim <- fit_pair(vr, tp, dr_trim)
  if (is.null(full) || is.null(trim)) return(NULL)

  # Also pull the VAF coefficient (Model 2 within-mutant dose) for side-by-side.
  beta_vaf <- function(mod) {
    cf <- tidy(mod, conf.int = TRUE) %>%
      filter(term == "vaf_pct_centered_in_mut")
    if (nrow(cf) == 0) tibble(beta = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
    else cf %>% transmute(beta = estimate, lo = conf.low, hi = conf.high, p = p.value)
  }

  bf <- beta_vaf(full$m2); bt <- beta_vaf(trim$m2)
  tibble(
    outcome          = lb,
    type             = tp,
    n_full           = full$n,
    n_trim           = trim$n,
    lrt_p_full       = full$lrt$`Pr(>Chi)`[2],
    lrt_p_trim       = trim$lrt$`Pr(>Chi)`[2],
    beta_vaf_full    = bf$beta,
    beta_vaf_trim    = bt$beta,
    beta_vaf_p_full  = bf$p,
    beta_vaf_p_trim  = bt$p
  )
})

# Concise printout so readers can eyeball robustness without opening the xlsx.
lrt_compare_print <- lrt_compare %>%
  mutate(
    lrt_p_full      = sprintf("%.3g", lrt_p_full),
    lrt_p_trim      = sprintf("%.3g", lrt_p_trim),
    beta_vaf_full   = sprintf("%.3f", beta_vaf_full),
    beta_vaf_trim   = sprintf("%.3f", beta_vaf_trim),
    beta_vaf_p_full = sprintf("%.3g", beta_vaf_p_full),
    beta_vaf_p_trim = sprintf("%.3g", beta_vaf_p_trim)
  )
cat("\n  LRT (Model 2 vs Model 1) and Model 2 VAF β — full vs |z|<=2 trimmed:\n")
print(lrt_compare_print, width = Inf)

# ═════════════════════════════════════════════════════════════════════════════
# Stats manifest fragment — edfig_vaf_outlier
# -----------------------------------------------------------------------------
# Surfaces the numeric handles the caption needs (Shapiro-Wilk W + p, n, %
# outliers, full-vs-trimmed cohort sizes) so the caption-time R substitution
# can render the actual values rather than asserting a panel annotation
# that does not exist. Per-figure rule: stats live in the caption, never
# on-panel. Written BEFORE the SuppTable10 append below so an xlsx-side
# failure does not block the caption-side stats from updating.
# ═════════════════════════════════════════════════════════════════════════════
sw_raw  <- shapiro.test(v)
sw_log  <- shapiro.test(log(v[v > 0]))
edfig_vaf_outlier_fragment <- list(
  n_vaf            = length(v),
  vaf_mean_pct     = m,
  vaf_sd_pct       = s,
  shapiro_W_raw    = unname(sw_raw$statistic),
  shapiro_p_raw    = sw_raw$p.value,
  shapiro_W_log    = unname(sw_log$statistic),
  shapiro_p_log    = sw_log$p.value,
  n_outliers       = sum(vaf_df$is_outlier),
  pct_outliers     = 100 * mean(vaf_df$is_outlier),
  n_full_cohort    = sum(!is.na(vaf_df$vaf_pct)),
  n_trim_cohort    = sum(!vaf_df$is_outlier)
)
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(section = "edfig_vaf_outlier",
                    stats   = edfig_vaf_outlier_fragment)

# ═════════════════════════════════════════════════════════════════════════════
# Supplementary Table: numeric support
# ═════════════════════════════════════════════════════════════════════════════
# 2026-06-14: ST7 is now a single mega-sheet written by 20_ST_supp_tables_8_9_10.R.
# Publish this script's four frames (distribution metadata + |z|>2 outlier
# sensitivity + dose-response models + LRT sensitivity) to a private RDS
# that 14 reads, same architecture as ST6.
cat("── Publishing ST7 outlier/dose/LRT frames to RDS for 14_supp_tables ──\n")

source(here("analysis", "helper_scripts", "utils.R"))

st7_int_path <- here("results", "stats", "_intermediates",
                     "st7_outlier_dose.rds")
dir.create(dirname(st7_int_path), recursive = TRUE, showWarnings = FALSE)

vaf_distribution <- tibble(
  metric = c("n", "mean_pct", "sd_pct", "median_pct", "IQR_low_pct", "IQR_high_pct",
             "min", "max", "skewness", "excess_kurtosis",
             "shapiro_W_raw", "shapiro_p_raw",
             "shapiro_W_log", "shapiro_p_log",
             "n_|z|>2", "pct_|z|>2", "n_|z|>3", "pct_|z|>3"),
  value  = c(length(v), m, s, median(v), quantile(v,.25), quantile(v,.75),
             min(v), max(v),
             mean((v-m)^3)/s^3, mean((v-m)^4)/s^4 - 3,
             shapiro.test(v)$statistic, shapiro.test(v)$p.value,
             shapiro.test(log(v[v>0]))$statistic, shapiro.test(log(v[v>0]))$p.value,
             sum(abs((v-m)/s)>2), 100*mean(abs((v-m)/s)>2),
             sum(abs((v-m)/s)>3), 100*mean(abs((v-m)/s)>3))
)

saveRDS(
  list(
    vaf_distribution = vaf_distribution,
    sens_rows        = sens_rows,       # |z|>2 outlier sensitivity (Spearman)
    dose_results     = dose_results,    # Model 1 / 2 / 3 coefficient table
    lrt_compare      = lrt_compare,     # LRT p, full vs |z|<=2 trimmed
    dose_outcomes    = dose_outcomes    # (var, label, type) for binary vs continuous
  ),
  st7_int_path
)
cat(sprintf("   → ST7 frames cached: %s\n",
            sub(paste0(here::here(), "/"), "", st7_int_path)))


cat("\n\u2713 ed_vaf_outlier complete.\n")
cat("   ST7 will be assembled by 20_ST_supp_tables_8_9_10.R from the cached RDS.\n")
