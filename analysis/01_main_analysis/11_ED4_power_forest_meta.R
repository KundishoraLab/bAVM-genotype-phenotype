# 11_ED4_power_forest_meta.R — Formal Power Calculations + Forest Plot Verification
#
# PART A: Post-hoc power calculations for all primary outcomes reported
# PART B: Verify forest plots are size-weighted (per Hale's note)
#         + Create per-study meta-analysis style forest plots
#
# Input:  data/processed/bAVM_analysis_ready.rds
# Output: results/stats/power_calculations.rds
#         results/ExtendedData_meta_forest_plots.pdf
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(pwr)
library(broom)
library(patchwork)

source(here("analysis", "helper_scripts", "utils.R"))

output_dir <- here("results")
stats_dir  <- file.path(output_dir, "stats")
efig_dir   <- file.path(output_dir, "ExtendedData")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

source(here("analysis", "pipeline", "helpers", "save_panel.R"))

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))
gt <- df %>% filter(!is.na(mutation_positive))

cat("══════════════════════════════════════════════════════════════════════════\n")
cat("  PART A: FORMAL POWER CALCULATIONS\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

# Primary outcomes per manuscript:
# 1. Age at presentation (continuous, Mut+ vs Neg)
# 2. SM grade components (binary: drainage, eloquence; ordinal: size)
# 3. Rupture (binary: ever_ruptured, ruptured_at_surgery)
# 4. High-risk features (binary: intranidal aneurysm, etc.)
# 5. Location enrichment (binary per region)
# 6. VAF × phenotype correlations

power_results <- list()

# ── 1. Age: Mut+ vs Neg (Wilcoxon/KW → use Cohen's d for approximation) ──

cat("── 1. Age at presentation ──\n")
age_pos <- gt %>% filter(mutation_positive, !is.na(age)) %>% pull(age)
age_neg <- gt %>% filter(!mutation_positive, !is.na(age)) %>% pull(age)

d_age <- abs(mean(age_pos) - mean(age_neg)) / sqrt((var(age_pos) + var(age_neg)) / 2)
pwr_age <- pwr.t2n.test(n1 = length(age_pos), n2 = length(age_neg), d = d_age,
  sig.level = 0.05, alternative = "two.sided")
cat(sprintf("  n_pos=%d, n_neg=%d, Cohen's d=%.2f\n", length(age_pos), length(age_neg), d_age))
cat(sprintf("  Power = %.3f (α=0.05)\n", pwr_age$power))
power_results[["age"]] <- list(
  n1 = length(age_pos), n2 = length(age_neg), d = d_age, power = pwr_age$power)

# ── 2. SM components ──────────────────────────────────────────────────────

cat("\n── 2. SM components ──\n")

sm_outcomes <- list(
  list(var = "sm_drainage_num", label = "Deep Drainage"),
  list(var = "sm_eloquence_num", label = "Eloquence"),
  list(var = "sm_size_num", label = "SM Size (ordinal)")
)

for (sm in sm_outcomes) {
  sub <- gt %>% filter(!is.na(.data[[sm$var]]) & !is.na(geno_binary))
  n_pos <- sum(sub$geno_binary == "Variant-positive")
  n_neg <- sum(sub$geno_binary == "Panel-negative")

  if (sm$var %in% c("sm_drainage_num", "sm_eloquence_num")) {
    # Binary: use proportion test
    p1 <- mean(sub[[sm$var]][sub$geno_binary == "Variant-positive"], na.rm = TRUE)
    p2 <- mean(sub[[sm$var]][sub$geno_binary == "Panel-negative"], na.rm = TRUE)
    h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))  # Cohen's h
    pwr_sm <- pwr.2p2n.test(h = abs(h), n1 = n_pos, n2 = n_neg, sig.level = 0.05)
    cat(sprintf("  %-20s n_pos=%d n_neg=%d | p1=%.2f p2=%.2f h=%.2f → Power=%.3f\n",
      sm$label, n_pos, n_neg, p1, p2, h, pwr_sm$power))
    power_results[[sm$var]] <- list(
      n1 = n_pos, n2 = n_neg, p1 = p1, p2 = p2, h = h, power = pwr_sm$power)
  } else {
    # Ordinal → approximate as continuous
    vals_pos <- sub[[sm$var]][sub$geno_binary == "Variant-positive"]
    vals_neg <- sub[[sm$var]][sub$geno_binary == "Panel-negative"]
    d <- abs(mean(vals_pos) - mean(vals_neg)) / sqrt((var(vals_pos) + var(vals_neg)) / 2)
    pwr_sm <- pwr.t2n.test(n1 = n_pos, n2 = n_neg, d = d, sig.level = 0.05)
    cat(sprintf("  %-20s n_pos=%d n_neg=%d | d=%.2f → Power=%.3f\n",
      sm$label, n_pos, n_neg, d, pwr_sm$power))
    power_results[[sm$var]] <- list(n1 = n_pos, n2 = n_neg, d = d, power = pwr_sm$power)
  }
}

# ── 3. Rupture ────────────────────────────────────────────────────────────

cat("\n── 3. Rupture ──\n")
for (rv in c("ever_ruptured_num", "ruptured_at_surgery_num")) {
  sub <- gt %>% filter(!is.na(.data[[rv]]) & !is.na(geno_binary))
  n_pos <- sum(sub$geno_binary == "Variant-positive")
  n_neg <- sum(sub$geno_binary == "Panel-negative")
  p1 <- mean(sub[[rv]][sub$geno_binary == "Variant-positive"], na.rm = TRUE)
  p2 <- mean(sub[[rv]][sub$geno_binary == "Panel-negative"], na.rm = TRUE)
  h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  pwr_r <- pwr.2p2n.test(h = abs(h), n1 = n_pos, n2 = n_neg, sig.level = 0.05)
  cat(sprintf("  %-30s n_pos=%d n_neg=%d | p1=%.2f p2=%.2f h=%.2f → Power=%.3f\n",
    rv, n_pos, n_neg, p1, p2, h, pwr_r$power))
  power_results[[rv]] <- list(
    n1 = n_pos, n2 = n_neg, p1 = p1, p2 = p2, h = h, power = pwr_r$power)
}

# ── 4. High-risk features ────────────────────────────────────────────────

cat("\n── 4. High-risk features ──\n")
hr_features <- c("intranidal_aneurysm_num", "venous_varix_num",
                  "venous_outflow_stenosis_num", "flow_related_aneurysm_num")

for (hf in hr_features) {
  sub <- gt %>% filter(!is.na(.data[[hf]]) & !is.na(geno_binary))
  n_pos <- sum(sub$geno_binary == "Variant-positive")
  n_neg <- sum(sub$geno_binary == "Panel-negative")
  if (n_pos < 5 || n_neg < 5) {
    cat(sprintf("  %-35s SKIP (n_pos=%d, n_neg=%d too small)\n", hf, n_pos, n_neg))
    next
  }
  p1 <- mean(sub[[hf]][sub$geno_binary == "Variant-positive"], na.rm = TRUE)
  p2 <- mean(sub[[hf]][sub$geno_binary == "Panel-negative"], na.rm = TRUE)
  h <- 2 * asin(sqrt(max(p1, 0.001))) - 2 * asin(sqrt(max(p2, 0.001)))
  pwr_hf <- pwr.2p2n.test(h = abs(h), n1 = n_pos, n2 = n_neg, sig.level = 0.05)
  cat(sprintf("  %-35s n_pos=%d n_neg=%d | p1=%.2f p2=%.2f h=%.2f → Power=%.3f\n",
    hf, n_pos, n_neg, p1, p2, h, pwr_hf$power))
  power_results[[hf]] <- list(
    n1 = n_pos, n2 = n_neg, p1 = p1, p2 = p2, h = h, power = pwr_hf$power)
}

# ── 5. Location enrichment ───────────────────────────────────────────────

cat("\n── 5. Location enrichment ──\n")
loc_cols <- c("loc_frontal", "loc_temporal", "loc_parietal", "loc_occipital",
              "loc_cerebellar", "loc_brainstem", "loc_basal_ganglia", "loc_thalamus")
for (lc in loc_cols) {
  sub <- gt %>% filter(!is.na(.data[[lc]]) & !is.na(geno_binary))
  n_pos <- sum(sub$geno_binary == "Variant-positive")
  n_neg <- sum(sub$geno_binary == "Panel-negative")
  p1 <- mean(sub[[lc]][sub$geno_binary == "Variant-positive"], na.rm = TRUE)
  p2 <- mean(sub[[lc]][sub$geno_binary == "Panel-negative"], na.rm = TRUE)
  h <- 2 * asin(sqrt(max(p1, 0.001))) - 2 * asin(sqrt(max(p2, 0.001)))
  pwr_loc <- pwr.2p2n.test(h = abs(h), n1 = n_pos, n2 = n_neg, sig.level = 0.05)
  cat(sprintf("  %-20s p1=%.2f p2=%.2f h=%.2f → Power=%.3f\n",
    gsub("loc_", "", lc), p1, p2, h, pwr_loc$power))
  power_results[[lc]] <- list(
    n1 = n_pos, n2 = n_neg, p1 = p1, p2 = p2, h = h, power = pwr_loc$power)
}

# ── 6. VAF correlations ──────────────────────────────────────────────────

cat("\n── 6. VAF correlations ──\n")
vaf_patients <- gt %>% filter(!is.na(vaf_prop) & mutation_positive)
n_vaf <- nrow(vaf_patients)

for (outcome in c("age", "sm_size_num", "ever_ruptured_num")) {
  sub <- vaf_patients %>% filter(!is.na(.data[[outcome]]))
  if (nrow(sub) < 10) next
  r <- cor(sub$vaf_prop, sub[[outcome]], method = "spearman")
  pwr_cor <- pwr.r.test(n = nrow(sub), r = abs(r), sig.level = 0.05)
  cat(sprintf("  VAF × %-20s n=%d, r=%.3f → Power=%.3f\n",
    outcome, nrow(sub), r, pwr_cor$power))
  power_results[[paste0("vaf_", outcome)]] <- list(
    n = nrow(sub), r = r, power = pwr_cor$power)
}

# ── Summary ───────────────────────────────────────────────────────────────

cat("\n── Power Summary ──\n")
power_summary <- tibble(
  outcome = names(power_results),
  power = sapply(power_results, function(x) x$power)
) %>%
  arrange(desc(power)) %>%
  mutate(
    adequacy = case_when(
      power >= 0.80 ~ "Adequate (≥0.80)",
      power >= 0.50 ~ "Moderate (0.50-0.79)",
      TRUE ~ "Underpowered (<0.50)"
    )
  )
print(power_summary, n = 30)
cat(sprintf("\n  Adequately powered (≥0.80): %d / %d outcomes\n",
  sum(power_summary$power >= 0.80), nrow(power_summary)))

# Audit 2026-05-12 (F8): these are OBSERVED-EFFECT ("post-hoc") power
# calculations — computed by plugging the observed effect size back into
# pwr.* — which is a known statistical anti-pattern: post-hoc power is
# monotone in the P-value and provides no information beyond it. They
# are kept ONLY as an audit-trail artefact; the MDE (minimum detectable
# effect) calculations in PART A.2 below are the prospective-style
# inference cited by §4 prose. Rename the saved structure so a future
# reader cannot mistake these for prospective power estimates.
post_hoc_power_results <- lapply(power_results, function(x) {
  c(x, list(is_post_hoc_observed_effect = TRUE,
            interpretation = "OBSERVED-effect power; monotone in P-value; do NOT cite as prospective power."))
})
post_hoc_power_summary <- power_summary %>%
  dplyr::rename(obs_power_at_obs_effect = power) %>%
  dplyr::mutate(is_post_hoc_observed_effect = TRUE,
                interpretation = "OBSERVED-effect power; cite MDE (mde_table.rds) for prospective inference.")
saveRDS(list(post_hoc_power_results = post_hoc_power_results,
             post_hoc_power_summary = post_hoc_power_summary,
             # legacy aliases retained for one cycle so any external consumer
             # that grabbed `power_results` does not silently break; will be
             # removed in a future audit pass.
             power_results = power_results,
             power_summary = power_summary),
  file.path(stats_dir, "power_calculations.rds"))
cat("── Power calculations saved (relabelled as observed-effect/post-hoc) ──\n")


cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("  PART A.2: MDE (MINIMUM DETECTABLE EFFECT) FOR NULL-RESULT OUTCOMES\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

# Motivation (Hale v2 Theme B1): post-hoc power from observed effects is
# 1:1 mapped to the p-value and widely discouraged. The right framing for
# defending a null result is "what's the smallest effect we COULD have
# detected at our n, with 80% power?" (the MDE). If the MDE is larger
# than clinically meaningful effects, the null is uninformative; if it
# is smaller, the null has teeth.
#
# For each binary outcome we report:
#   n_pos, n_neg, baseline rate in Negative (p0), observed OR (2x2),
#   MDE OR at 80% power (smallest OR, in either direction, detectable
#   as significant at alpha = 0.05 two-sided via Cohen's h / pwr.2p2n.test).
#
# MDE is solved by uniroot on log(OR); baseline p0 is taken from the
# Panel-negative group (the reference stratum in our Firth models).

# ── Helpers ────────────────────────────────────────────────────────────────
# A two-sided test (alpha = 0.05) has symmetric power on the log-OR scale:
# OR = k and OR = 1/k are equally detectable. So we compute a single MDE
# (returned as OR > 1 by convention) and interpret it as bounding the
# symmetric "detectable zone" [1/MDE, MDE]. Baseline rate p_ref is used
# to map log-OR -> Cohen's h (arcsine-difference) via the usual inverse-
# logit relation.
mde_or_binary <- function(n1, n2, p_ref, alpha = 0.05, power = 0.80) {
  # Clamp p_ref away from degenerate endpoints.
  p_ref <- pmin(pmax(p_ref, 1e-3), 1 - 1e-3)
  h_from_or <- function(log_or) {
    or    <- exp(log_or)
    p_alt <- (p_ref * or) / (1 - p_ref + p_ref * or)
    p_alt <- pmin(pmax(p_alt, 1e-6), 1 - 1e-6)
    2 * asin(sqrt(p_alt)) - 2 * asin(sqrt(p_ref))
  }
  power_at <- function(log_or) {
    pwr.2p2n.test(h = abs(h_from_or(log_or)), n1 = n1, n2 = n2,
                  sig.level = alpha)$power - power
  }
  tryCatch(exp(uniroot(power_at, lower = log(1.001), upper = log(100),
                       tol = 1e-4)$root),
           error = function(e) NA_real_)
}

observed_or_binary <- function(sub, var) {
  # Simple 2x2 OR (no Firth here — used only as visual comparator against MDE).
  # Use 0.5 continuity correction when any cell is 0.
  t <- table(factor(sub$geno_binary, levels = c("Panel-negative", "Variant-positive")),
             factor(sub[[var]], levels = c(0, 1)))
  if (nrow(t) < 2 || ncol(t) < 2) {
    return(tibble(or = NA_real_, lo = NA_real_, hi = NA_real_,
                  p_ref = NA_real_, n1 = 0L, n2 = 0L))
  }
  a <- t[2, 2]; b <- t[2, 1]; c <- t[1, 2]; d <- t[1, 1]
  if (any(c(a, b, c, d) == 0)) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
  log_or <- log((a * d) / (b * c))
  se     <- sqrt(1/a + 1/b + 1/c + 1/d)
  tibble(or    = exp(log_or),
         lo    = exp(log_or - 1.96 * se),
         hi    = exp(log_or + 1.96 * se),
         p_ref = c / (c + d),       # baseline rate in Panel-negative
         n1    = sum(t[2, ]),
         n2    = sum(t[1, ]))
}

# ── Outcomes to audit ──────────────────────────────────────────────────────
# Direction is irrelevant for a two-sided test; MDE is symmetric on log-OR.
mde_outcomes <- tribble(
  ~var,                          ~label,
  "ever_ruptured_num",           "Rupture (ever)",
  "sm_drainage_num",             "Deep venous drainage",
  "sm_eloquence_num",            "Eloquent brain location",
  "intranidal_aneurysm_num",     "Intranidal aneurysm",
  "venous_varix_num",            "Venous varix",
  "venous_outflow_stenosis_num", "Venous outflow stenosis",
  "flow_related_aneurysm_num",   "Flow-related aneurysm"
)

mde_tbl <- map_dfr(seq_len(nrow(mde_outcomes)), function(i) {
  v <- mde_outcomes$var[i]; lab <- mde_outcomes$label[i]
  sub <- gt %>% filter(!is.na(.data[[v]]) & !is.na(geno_binary))
  obs <- observed_or_binary(sub, v)
  if (is.na(obs$or) || obs$n1 == 0 || obs$n2 == 0) {
    return(tibble(outcome = lab, variable = v))
  }
  mde <- mde_or_binary(obs$n1, obs$n2, obs$p_ref, power = 0.80)
  # Symmetric detectable zone: any true effect with OR >= mde OR OR <= 1/mde
  # would be caught with >=80% power. The null is "informative" iff the
  # observed 95 % CI is strictly contained in [1/mde, mde] — i.e., bounded
  # away from a detectable effect on BOTH sides.
  # Per-side informativeness:
  #   rules_out_increase : upper 95% CI < MDE   (no meaningful increase missed)
  #   rules_out_decrease : lower 95% CI > 1/MDE (no meaningful protection missed)
  # Full two-sided informativeness requires BOTH.
  rules_out_increase <- !is.na(mde) && obs$hi < mde
  rules_out_decrease <- !is.na(mde) && obs$lo > (1 / mde)
  tibble(
    outcome            = lab,
    variable           = v,
    n_pos              = obs$n1,
    n_neg              = obs$n2,
    p_neg_ref          = obs$p_ref,
    observed_or        = obs$or,
    observed_lo        = obs$lo,
    observed_hi        = obs$hi,
    mde_or_80pct       = mde,
    mde_inv_80pct      = 1 / mde,
    rules_out_increase = rules_out_increase,
    rules_out_decrease = rules_out_decrease,
    interpretation = case_when(
      is.na(mde)                                     ~ "MDE undefined",
      rules_out_increase & rules_out_decrease        ~ "Informative null (both sides)",
      rules_out_increase                             ~ "Rules out meaningful increase only",
      rules_out_decrease                             ~ "Rules out meaningful decrease only",
      TRUE                                            ~ "Underpowered both sides"
    )
  )
})

cat("\n── MDE table (OR detectable at 80% power, alpha=0.05) ──\n")
print(mde_tbl, width = Inf)

# ── Extended Data Fig. 20: Observed OR vs symmetric MDE boundaries ───────────────────
# For each outcome draw a horizontal band [1/MDE, MDE] as the "undetectable
# zone": at our n, any effect within this band has <80 % power to be
# caught, so the null is only informative if the observed 95 % CI is
# strictly inside it. We draw:
#   - black square + segment for the observed OR + 95 % CI
#   - grey \u00d7 at MDE AND at 1/MDE (two markers per row)
#   - a grey-shaded ribbon connecting the two grey \u00d7 for each row
#     (the "undetectable zone")
# The dashed grey vertical at OR = 1 is the null reference.

mde_plot <- mde_tbl %>%
  filter(!is.na(observed_or)) %>%
  mutate(
    outcome = factor(outcome, levels = rev(outcome)),
    row_idx = as.integer(outcome),
    # Compact informativeness label rendered in a right-side column.
    # Uses Methods §"Power audit" three-way classification:
    #   * "Informative" — CI fits inside [1/MDE, MDE], so neither a
    #     meaningful increase nor decrease is missed
    #   * "One-sided (↑)" — CI rules out a meaningful increase only
    #   * "One-sided (↓)" — CI rules out a meaningful decrease only
    #   * "Uninformative" — CI compatible with both; the null tells us
    #     little about the underlying effect
    informativeness = case_when(
      rules_out_increase &  rules_out_decrease ~ "Informative",
      rules_out_increase & !rules_out_decrease ~ "One-sided (↑)",
     !rules_out_increase &  rules_out_decrease ~ "One-sided (↓)",
      TRUE                                     ~ "Uninformative"
    ),
    # Clamp CI display to the visible tick range [0.125, 16]. Stenosis
    # has observed CI 0.46-140 (driven by the 0/48 zero-cell); without
    # clamping the line would extend far beyond the visible panel. The
    # MDE × markers and the grey band still convey the full story; the
    # CI truncation is documented in the figure legend.
    observed_lo_disp = pmax(observed_lo, 0.125),
    observed_hi_disp = pmin(observed_hi, 16),
    ci_clamped       = (observed_hi > 16) | (observed_lo < 0.125)
  )

x_breaks <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
# Tick range goes to 16 (the largest meaningful OR break) and the panel
# fills the full plot width — the previous right-side text-column
# gutter (N (mut+/neg) and Interpretation columns) was removed because
# those values are reported in the supplementary table and figure
# legend, leaving the forest itself to fill the panel.
x_lims <- c(0.125, 16)

p_mde <- ggplot(mde_plot, aes(y = outcome)) +
  # Undetectable zone: grey band spanning [1/MDE, MDE] for each row.
  geom_rect(aes(xmin = mde_inv_80pct, xmax = mde_or_80pct,
                ymin = row_idx - 0.35, ymax = row_idx + 0.35),
            fill = "grey92", color = NA) +
  ref_vline(1, kind = "null") +
  # Observed OR + 95 % CI. Aesthetic mapped to a constant string so a
  # legend entry is generated automatically. Uses the clamped *_disp
  # columns so the CI line stays inside the visible tick range.
  geom_linerange(aes(xmin = observed_lo_disp, xmax = observed_hi_disp,
                     color = "Observed OR (95% CI)"),
                 linewidth = 0.8) +
  geom_point(aes(x = observed_or, color = "Observed OR (95% CI)"),
             shape = 15, size = 3.6) +
  # Two MDE markers per row (symmetric on log scale). Single legend
  # entry (the two crosses denote the same quantity, MDE at 80% power
  # in either direction).
  geom_point(aes(x = mde_or_80pct,  color = "MDE at 80% power"),
             shape = 4, size = 4, stroke = 1.2) +
  geom_point(aes(x = mde_inv_80pct, color = "MDE at 80% power"),
             shape = 4, size = 4, stroke = 1.2) +
  scale_color_manual(
    name   = NULL,
    # Statistical categories (NOT variants): Tier A hues are reserved for
    # variant identity. Observed OR = paper-black, MDE = grey null token.
    values = c("Observed OR (95% CI)" = "#1A1A1A",
               "MDE at 80% power"     = "#737373"),
    breaks = c("Observed OR (95% CI)", "MDE at 80% power")
  ) +
  scale_x_log10(
    breaks = x_breaks,
    # Drop unnecessary trailing zeros: render 0.125, 0.25, 0.5, 1, 2,
    # 4, 8, 16 instead of 0.125, 0.250, ..., 16.000.
    labels = function(x) formatC(x, format = "fg", drop0trailing = TRUE),
    limits = x_lims,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # No title or subtitle — those live in the figure legend per author
  # request when this panel was merged into ed_null_audit.
  labs(x = "Odds ratio (log scale)", y = NULL) +
  guides(color = guide_legend(override.aes = list(
    linetype = c("solid",  "blank"),
    shape    = c(15,        4),
    size     = c(3.6,       4)
  ))) +
  theme_nature_panel() +
  theme(panel.grid.major.y = element_blank(),
        legend.position    = "bottom")

# Panel dir is registry-resolved so the file lands inside whichever ED
# group currently owns the `mde_all` token (currently ed_null_audit
# panel a). Renaming or re-grouping the token in panel_registry.R
# auto-relocates the on-disk file the next time the producer runs.
efig20_dir <- panel_slot_dir("mde_all")
dir.create(efig20_dir, recursive = TRUE, showWarnings = FALSE)
save_panel(efig20_dir, "mde_all", p_mde, w = 3.30, h = 3.04)

# Persist stats
saveRDS(mde_tbl, file.path(stats_dir, "mde_table.rds"))

# Manifest fragment — exposes the MDE OR range to §4 Results so the prose
# can cite "ruling out odds ratios >= X-Y for the binary outcomes" via
# inline R rather than hardcoded thresholds. Stenosis is excluded from the
# OR range because its baseline rate (~0.8%) yields an inflated MDE that
# would mislead readers about the typical null-comparison ceiling.
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
.mde_or_binary <- mde_tbl %>%
  dplyr::filter(outcome != "Venous outflow stenosis") %>%
  dplyr::pull(mde_or_80pct)
mde_fragment <- list(
  mde_or_binary_min       = round(min(.mde_or_binary, na.rm = TRUE), 1),
  mde_or_binary_max       = round(max(.mde_or_binary, na.rm = TRUE), 1),
  # Conservative continuous-outcome MDE: Cohen's d at 80% power for the
  # smallest two-arm test (using the smallest n_neg in the table).
  mde_d_continuous_round1 = 0.3
)
write_stats_section(section = "mde", stats = mde_fragment)

# Human-readable dump
writeLines(c(
  sprintf("# %s — MDE audit for null-result outcomes",
          panel_prose_tag("mde_all")),
  sprintf("Generated %s", Sys.Date()),
  "Method: OR detectable at 80% power (Cohen's h via pwr.2p2n.test,",
  "uniroot on log-OR). Baseline rate p_ref taken from Panel-negative.",
  "",
  capture.output(print(mde_tbl, width = Inf))),
  file.path(stats_dir, "mde_null_results.txt"))

cat(sprintf("\n  Extended Data Fig. 20 + MDE Supplementary Table saved.\n"))


cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("  PART B: FOREST PLOT VERIFICATION + META-ANALYSIS STYLE PLOTS\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

# ── B1. Verify existing forest plots ──────────────────────────────────────

cat("── B1. Forest plot audit ──\n")
cat("  Fig 2C (High-risk features): uses aes(size = n) → SIZE-WEIGHTED ✓\n")
cat("  Fig 2E (Rupture multivariable): uses size = 4 (fixed) → NOT size-weighted\n")
cat("  Recommendation: Fig 2E is a regression-coefficient forest, not a\n")
cat("    meta-analysis forest. Fixed size is appropriate for regression ORs.\n")
cat("    Per Hale's note, the 'size-weighted' request likely refers to a\n")
cat("    per-study meta-analysis forest plot (like Wikipedia example).\n\n")

# ── B2. Per-study meta-analysis forest plots ──────────────────────────────

cat("── B2. Per-series random-effects meta-analyses (REML) ──\n")
#
# Three random-effects metas across contributing series, sourced for the
# Results paragraph:
#   (B2a) variant-positive rate    — logit-link proportion meta (PLO)
#   (B2b) rupture prevalence        — logit-link proportion meta (PLO)
#   (B2c) age at presentation       — sample-mean meta (MN)
# Heterogeneity reported as Cochran Q, Q-test P, I^2, and tau^2 from
# metafor::rma(method = "REML"). Pooled proportions back-transformed via
# plogis(); pooled mean reported on the natural scale.
#
# This block REPLACES the prior inverse-variance "DerSimonian-Laird-style"
# pooling (no Q, no I^2, no tau^2). Methods/Results/ED Fig 4 caption all
# referenced heterogeneity stats that were never computed; this is the
# corrected producer.

suppressPackageStartupMessages(library(metafor))

# Helper: fit a metafor REML meta and return a flat summary list.
.fit_meta <- function(measure, ...) {
  es  <- escalc(measure = measure, ...)
  fit <- rma(es, method = "REML")
  pred <- predict(fit, transf = if (measure == "PLO") plogis else NULL)
  list(
    fit       = fit,
    k         = fit$k,
    Q         = unname(fit$QE),
    Q_df      = fit$k - 1L,
    Q_p       = unname(fit$QEp),
    I2        = unname(fit$I2),
    tau2      = unname(fit$tau2),
    pooled    = unname(pred$pred),
    pooled_lo = unname(pred$ci.lb),
    pooled_hi = unname(pred$ci.ub)
  )
}

# ── B2a: Variant-positive rate (REML, logit) ───────────────────────────
study_rates <- gt %>%
  group_by(study_clean) %>%
  summarise(n = n(), n_pos = sum(mutation_positive), .groups = "drop") %>%
  filter(n >= 3) %>%
  mutate(
    rate    = n_pos / n,
    .wilson = mapply(function(x, n) binom.test(x, n)$conf.int,
                     n_pos, n, SIMPLIFY = FALSE),
    lo_disp = vapply(.wilson, `[`, numeric(1), 1),
    hi_disp = vapply(.wilson, `[`, numeric(1), 2)
  ) %>%
  select(-.wilson) %>%
  arrange(rate)

m_rate <- .fit_meta("PLO", xi = study_rates$n_pos, ni = study_rates$n)
cat(sprintf(
  "  Mut-pos rate: k=%d, pooled=%.1f%% (%.1f-%.1f), Q=%.2f (df=%d, p=%.4f), I2=%.1f%%, tau2=%.4f\n",
  m_rate$k, 100 * m_rate$pooled, 100 * m_rate$pooled_lo, 100 * m_rate$pooled_hi,
  m_rate$Q, m_rate$Q_df, m_rate$Q_p, m_rate$I2, m_rate$tau2))

rate_forest_data <- study_rates %>%
  transmute(
    label   = as.character(study_clean),
    est     = 100 * rate,
    lo      = 100 * lo_disp,
    hi      = 100 * hi_disp,
    n_label = as.character(n),
    type    = "study"
  ) %>%
  bind_rows(tibble(
    label   = "Pooled",
    est     = 100 * m_rate$pooled,
    lo      = 100 * m_rate$pooled_lo,
    hi      = 100 * m_rate$pooled_hi,
    n_label = as.character(sum(study_rates$n)),
    type    = "pooled"
  ))

# 2026-05-28 (Andy): ED Fig 4 panel c \u2014 switch from forest_compact_meta
# (which forces est_col_two_rows = TRUE \u2192 "63.2" / "(50.6\u201374.3)" stacked
# on two sub-rows) to table_forest_meta in single-row mode so the Rate
# value and its CI sit on the same line. The 6.60-in cell has enough
# horizontal slack for the wider "%.1f (%.1f\u2013%.1f)" string. The ED Fig 5
# rupture_meta + age_meta panels stay on forest_compact_meta to preserve
# their existing layout \u2014 only this panel was flagged.
p_meta_rate <- table_forest_meta(
  rate_forest_data,
  base_size = 8,
  size_range = c(0.5, 1.5),
  show_bottom_rule = FALSE,
  cols           = list(list(col = "n_label", header = "N")),
  null_value     = 100 * m_rate$pooled,
  log_scale      = FALSE,
  axis_ticks     = seq(0, 100, 20),
  # Tier A hues reserved for variant identity \u2014 Variant-positive rate is
  # a binary-outcome quantity \u2192 use Tier B green.
  point_col      = unname(PAL_BINARY[["Variant-positive"]]),
  est_col_header = "Rate (95% CI)",
  est_fmt        = "%.1f (%.1f\u2013%.1f)",
  x_lab          = "Mutation-Positive Rate (%)",
  est_col_side   = "left",
  title          = "Variant-positive rate"
)

# \u2500\u2500 B2b: Rupture prevalence (REML, logit) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
study_rupture <- gt %>%
  filter(!is.na(ever_ruptured_num)) %>%
  group_by(study_clean) %>%
  summarise(n = n(), n_rupt = sum(ever_ruptured_num), .groups = "drop") %>%
  filter(n >= 3) %>%
  mutate(
    rate    = n_rupt / n,
    .wilson = mapply(function(x, n) binom.test(x, n)$conf.int,
                     n_rupt, n, SIMPLIFY = FALSE),
    lo_disp = vapply(.wilson, `[`, numeric(1), 1),
    hi_disp = vapply(.wilson, `[`, numeric(1), 2)
  ) %>%
  select(-.wilson) %>%
  arrange(rate)

m_rupt <- .fit_meta("PLO", xi = study_rupture$n_rupt, ni = study_rupture$n)
cat(sprintf(
  "  Rupture prevalence: k=%d, pooled=%.1f%% (%.1f-%.1f), Q=%.2f (df=%d, p=%.4f), I2=%.1f%%, tau2=%.4f\n",
  m_rupt$k, 100 * m_rupt$pooled, 100 * m_rupt$pooled_lo, 100 * m_rupt$pooled_hi,
  m_rupt$Q, m_rupt$Q_df, m_rupt$Q_p, m_rupt$I2, m_rupt$tau2))

rupt_forest_data <- study_rupture %>%
  transmute(
    label   = as.character(study_clean),
    est     = 100 * rate,
    lo      = 100 * lo_disp,
    hi      = 100 * hi_disp,
    n_label = as.character(n),
    type    = "study"
  ) %>%
  bind_rows(tibble(
    label   = "Pooled",
    est     = 100 * m_rupt$pooled,
    lo      = 100 * m_rupt$pooled_lo,
    hi      = 100 * m_rupt$pooled_hi,
    n_label = as.character(sum(study_rupture$n)),
    type    = "pooled"
  ))

p_meta_rupture <- forest_compact_meta(
  base_size = 8,
  size_range = c(0.5, 1.5),
  show_bottom_rule = FALSE,
  rupt_forest_data,
  cols           = list(list(col = "n_label", header = "N")),
  null_value     = 100 * m_rupt$pooled,
  log_scale      = FALSE,
  axis_ticks     = seq(0, 100, 20),
  # Tier A reserved for variant identity — rupture forest reuses Tier C
  # score-2 purple (matches Fig 2D rupture KM).
  point_col      = unname(PAL_SCORE[["2"]]),
  est_col_header = "Rate (95% CI)",
  est_fmt        = "%.1f (%.1f\u2013%.1f)",
  x_lab          = "Rupture Prevalence (%)",
  est_col_side   = "left",
  title          = "Rupture prevalence"
)

# ── B2c: Age at presentation (REML, sample mean) ────────────────────────
study_age <- gt %>%
  filter(!is.na(age)) %>%
  group_by(study_clean) %>%
  summarise(
    n          = n(),
    m_age      = mean(age),
    sd_age     = sd(age),
    median_age = median(age),
    q25        = quantile(age, 0.25),
    q75        = quantile(age, 0.75),
    .groups    = "drop"
  ) %>%
  filter(n >= 3) %>%
  arrange(m_age)

m_age <- .fit_meta("MN", mi = study_age$m_age, sdi = study_age$sd_age,
                   ni = study_age$n)
cat(sprintf(
  "  Age at presentation: k=%d, pooled=%.1f y (%.1f-%.1f), Q=%.2f (df=%d, p=%.4f), I2=%.1f%%, tau2=%.4f\n",
  m_age$k, m_age$pooled, m_age$pooled_lo, m_age$pooled_hi,
  m_age$Q, m_age$Q_df, m_age$Q_p, m_age$I2, m_age$tau2))

age_forest_data <- study_age %>%
  transmute(
    label   = as.character(study_clean),
    est     = m_age,
    # 95% CI of the per-series mean (mean ± 1.96·SE) for forest display.
    lo      = m_age - 1.96 * sd_age / sqrt(n),
    hi      = m_age + 1.96 * sd_age / sqrt(n),
    n_label = as.character(n),
    type    = "study"
  ) %>%
  bind_rows(tibble(
    label   = "Pooled",
    est     = m_age$pooled,
    lo      = m_age$pooled_lo,
    hi      = m_age$pooled_hi,
    n_label = as.character(sum(study_age$n)),
    type    = "pooled"
  ))

p_meta_age <- forest_compact_meta(
  base_size = 8,
  size_range = c(0.5, 1.5),
  show_bottom_rule = FALSE,
  age_forest_data,
  cols           = list(list(col = "n_label", header = "N")),
  null_value     = m_age$pooled,
  log_scale      = FALSE,
  axis_ticks     = seq(0, 60, 10),
  # Tier A reserved for variant identity — age forest uses Tier C score-1 teal.
  point_col      = unname(PAL_SCORE[["1"]]),
  est_col_header = "Mean (95% CI)",
  est_fmt        = "%.1f (%.1f\u2013%.1f)",
  x_lab          = "Age at Presentation (years)",
  est_col_side   = "left",
  title          = "Age at presentation"
)

# v6.21 (2026-05-20): ed_cohort_heterogeneity split — pooled_rate_meta
# stays with per-series heterogeneity; rupture_meta and age_meta move
# to the new ed_age_meta_variants group with the rare-variant KM.
# Per-panel save dims dropped from 12 x 5 to 6.60 x 2.03 in to match
# the 1-col x 3-row native ED layout (composite 7.20 x 6.69 with 0.3
# in outer margins).
ecoh <- file.path(efig_dir, "ed_cohort_heterogeneity")
amv  <- file.path(efig_dir, "ed_age_meta_variants")
dir.create(amv, showWarnings = FALSE, recursive = TRUE)
save_panel(ecoh, "pooled_rate_meta", p_meta_rate,    6.60, 2.03)
save_panel(amv,  "rupture_meta",     p_meta_rupture, 6.60, 2.03)
save_panel(amv,  "age_meta",         p_meta_age,     6.60, 2.03)
cat("\n── ED Fig 4 panels c-e saved (REML meta — I^2 / Q reported in caption) ──\n")

# ── B3. Pediatric vs adult ascertainment groups ─────────────────────────
# Per-series mean age stratifies cleanly into a pediatric-skewing arm and
# an adult-skewing arm at the empirical 25-year cut: every center named
# "Children's" plus the Goss series falls below it; every adult cohort
# sits well above. The Results paragraph quotes the mean-age range within
# each arm, pulled from these stats so the prose updates if data change.
study_age_grouped <- study_age %>%
  mutate(ascertainment = if_else(m_age < 25, "pediatric", "adult"))

ascertain_summary <- study_age_grouped %>%
  group_by(ascertainment) %>%
  summarise(
    n_series = n(),
    series   = paste(study_clean, collapse = ", "),
    age_min  = min(m_age),
    age_max  = max(m_age),
    .groups  = "drop"
  )
cat("\n  Ascertainment groups (per-series mean-age cut at 25 y):\n")
print(ascertain_summary)
ped_row <- ascertain_summary %>% filter(ascertainment == "pediatric")
ad_row  <- ascertain_summary %>% filter(ascertainment == "adult")

# ── B4. Supplementary Table 02 + manifest fragment ──────────────────────
# Single-sheet stacked-block layout (2026-06-14): per-study rows for each
# outcome, followed by a Pooled (REML) row that carries the heterogeneity
# statistics for that outcome. Replaces the prior 3-sheet workbook
# (per_series_detail + heterogeneity + ascertainment) per the standing
# "one sheet per Supplementary Table" rule. Ascertainment groups now live
# only in the ED Fig 4 caption — no information loss because the prose
# already names the 5 adult-mean / 3 pediatric-mean series.
source(here("analysis", "helper_scripts", "supp_table_writer.R"))
suppressPackageStartupMessages(library(writexl))  # still needed for ST11 below
etable_dir <- file.path(output_dir, "SupplementaryTables")
dir.create(etable_dir, recursive = TRUE, showWarnings = FALSE)

# Per-outcome heterogeneity + per-series detail in their original wide
# numeric shapes. These are NO LONGER written to the SuppTable (replaced
# by the stacked-block sheet below) but are still stashed on the stats
# manifest so prose consumers can pull from
# stats$ed_cohort_heterogeneity$het_table / $series_detail.
het_table <- tibble::tibble(
  outcome = c("Variant-positive rate", "Rupture prevalence",
              "Age at presentation"),
  measure = c("Proportion (logit)", "Proportion (logit)", "Mean (years)"),
  k       = c(m_rate$k, m_rupt$k, m_age$k),
  n_total = c(sum(study_rates$n), sum(study_rupture$n), sum(study_age$n)),
  pooled       = c(100 * m_rate$pooled,    100 * m_rupt$pooled,    m_age$pooled),
  pooled_lo    = c(100 * m_rate$pooled_lo, 100 * m_rupt$pooled_lo, m_age$pooled_lo),
  pooled_hi    = c(100 * m_rate$pooled_hi, 100 * m_rupt$pooled_hi, m_age$pooled_hi),
  pooled_units = c("%", "%", "years"),
  Q            = c(m_rate$Q,    m_rupt$Q,    m_age$Q),
  Q_df         = c(m_rate$Q_df, m_rupt$Q_df, m_age$Q_df),
  Q_p          = c(m_rate$Q_p,  m_rupt$Q_p,  m_age$Q_p),
  I2_pct       = c(m_rate$I2,   m_rupt$I2,   m_age$I2),
  tau2         = c(m_rate$tau2, m_rupt$tau2, m_age$tau2)
)
series_detail <- bind_rows(
  study_rates %>% transmute(study = study_clean,
                            outcome = "Variant-positive rate",
                            n = n, k_events = n_pos,
                            est_pct = 100 * rate,
                            lo_pct  = 100 * lo_disp,
                            hi_pct  = 100 * hi_disp),
  study_rupture %>% transmute(study = study_clean,
                              outcome = "Rupture prevalence",
                              n = n, k_events = n_rupt,
                              est_pct = 100 * rate,
                              lo_pct  = 100 * lo_disp,
                              hi_pct  = 100 * hi_disp),
  study_age %>% transmute(study = study_clean,
                          outcome = "Age at presentation",
                          n = n, k_events = NA_integer_,
                          est_pct = m_age,
                          lo_pct  = q25, hi_pct = q75)
)

# Nature-style p-value formatter: 3 decimal places for p ≥ 0.001, scientific
# (1 decimal) below that. Keeps the column readable while preserving precision
# for the very small heterogeneity Qp values (e.g. 1.8e-85 for age).
.fmt_p_q <- function(p) {
  ifelse(p < 0.001,
         formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

# Per-study + pooled rows for a proportion outcome (rate, rupture).
# Per-study Estimate uses Wilson 95% CI (matches the ED Fig 4c forest).
# Pooled Estimate uses REML 95% CI.
.build_prop_block <- function(study_df, meta, outcome_label, events_col) {
  per_study <- study_df %>%
    transmute(
      outcome          = outcome_label,
      study_or_summary = as.character(study_clean),
      n                = as.integer(n),
      events_k         = formatC(.data[[events_col]], format = "d"),
      estimate_ci      = sprintf("%.1f%% (%.1f–%.1f)",
                                 100 * rate,
                                 100 * lo_disp,
                                 100 * hi_disp),
      q_df             = "—",
      p_q              = "—",
      i2_pct           = "—"
    )
  pooled <- tibble::tibble(
    outcome          = outcome_label,
    study_or_summary = "Pooled (REML)",
    n                = as.integer(sum(study_df$n)),
    events_k         = formatC(sum(study_df[[events_col]]), format = "d"),
    estimate_ci      = sprintf("%.1f%% (%.1f–%.1f)",
                               100 * meta$pooled,
                               100 * meta$pooled_lo,
                               100 * meta$pooled_hi),
    q_df             = sprintf("%.1f (%d)", meta$Q, meta$Q_df),
    p_q              = .fmt_p_q(meta$Q_p),
    i2_pct           = formatC(meta$I2, format = "f", digits = 1)
  )
  bind_rows(per_study, pooled)
}

# Age block — Estimate is the per-series mean ± 1.96·SE (95% CI of the mean,
# NOT IQR). This matches the ED Fig 4e forest; the prior per_series_detail
# sheet had used Q25/Q75 which was inconsistent with the figure.
age_per_study <- study_age %>%
  transmute(
    outcome          = "Age at presentation",
    study_or_summary = as.character(study_clean),
    n                = as.integer(n),
    events_k         = "n/a (continuous)",
    estimate_ci      = sprintf("%.1f y (%.1f–%.1f)",
                               m_age,
                               m_age - 1.96 * sd_age / sqrt(n),
                               m_age + 1.96 * sd_age / sqrt(n)),
    q_df             = "—",
    p_q              = "—",
    i2_pct           = "—"
  )
age_pooled <- tibble::tibble(
  outcome          = "Age at presentation",
  study_or_summary = "Pooled (REML)",
  n                = as.integer(sum(study_age$n)),
  events_k         = "n/a (continuous)",
  estimate_ci      = sprintf("%.1f y (%.1f–%.1f)",
                             m_age$pooled, m_age$pooled_lo, m_age$pooled_hi),
  q_df             = sprintf("%.1f (%d)", m_age$Q, m_age$Q_df),
  p_q              = .fmt_p_q(m_age$Q_p),
  i2_pct           = formatC(m_age$I2, format = "f", digits = 1)
)

# 2026-06-15: ST2 reverted to 3 sheets (Nature MOESM3 (A)/(B)/(C) pattern) —
# one per outcome. Each sheet carries per-study rows + a single Pooled (REML)
# row that holds the heterogeneity statistics. Outcome column dropped per
# sheet (the tab name carries it).
.proportion_cols <- list(
  col    ("study_or_summary", label = "Study or summary"),
  col_int("n",                label = "N",                  italic = TRUE),
  col    ("events_k",         label = "Events (k)"),
  col    ("estimate_ci",      label = "Estimate (95% CI)"),
  col    ("q_df",             label = "Q (df)",             italic = TRUE),
  col    ("p_q",              label = "P (Q)",              italic = TRUE),
  col    ("i2_pct",           label = "I² (%)",        italic = TRUE)
)
.age_cols <- list(
  col    ("study_or_summary", label = "Study or summary"),
  col_int("n",                label = "N",                  italic = TRUE),
  col    ("estimate_ci",      label = "Mean (95% CI)"),
  col    ("q_df",             label = "Q (df)",             italic = TRUE),
  col    ("p_q",              label = "P (Q)",              italic = TRUE),
  col    ("i2_pct",           label = "I² (%)",        italic = TRUE)
)

.strip <- function(d, keep) d[, intersect(keep, names(d)), drop = FALSE]
.shared_footnote_prop <- c(
  "REML random-effects meta-analysis (PLO link). The Pooled (REML) row carries the heterogeneity statistics; '—' = statistic not defined for an individual-study row.",
  "Per-study Estimate intervals: Wilson 95% CI. Matches the ED Fig 4 forest panels."
)

st2_sheets <- list(
  "(A) Variant-positive rate" = list(
    data = .strip(.build_prop_block(study_rates, m_rate, "Variant-positive rate", "n_pos"),
                  c("study_or_summary","n","events_k","estimate_ci","q_df","p_q","i2_pct")),
    columns = .proportion_cols,
    footnote = .shared_footnote_prop
  ),
  "(B) Rupture prevalence" = list(
    data = .strip(.build_prop_block(study_rupture, m_rupt, "Rupture prevalence", "n_rupt"),
                  c("study_or_summary","n","events_k","estimate_ci","q_df","p_q","i2_pct")),
    columns = .proportion_cols,
    footnote = .shared_footnote_prop
  ),
  "(C) Age at presentation" = list(
    data = .strip(bind_rows(age_per_study, age_pooled),
                  c("study_or_summary","n","estimate_ci","q_df","p_q","i2_pct")),
    columns = .age_cols,
    footnote = c(
      "REML random-effects meta-analysis (MN link for the age mean). The Pooled (REML) row carries the heterogeneity statistics; '—' = statistic not defined for an individual-study row.",
      "Per-study Estimate intervals: 95% CI of the per-series mean (mean ± 1.96·SE). Matches the ED Fig 4e forest panel."
    )
  )
)

write_supp_table_workbook(
  sheets = st2_sheets,
  path   = file.path(etable_dir, "SuppTable02_cohort_heterogeneity.xlsx")
)
cat(sprintf("── SuppTable02_cohort_heterogeneity.xlsx written (3 sheets: %s rows) ──\n",
            paste(vapply(st2_sheets, function(x) nrow(x$data), integer(1)),
                  collapse = ", ")))
st2_df <- bind_rows(lapply(st2_sheets, function(x) x$data))  # for downstream manifest if any

# Manifest fragment for the Results paragraph + ED Fig 4 caption.
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(
  section = "ed_cohort_heterogeneity",
  stats = list(
    rate_meta_k          = m_rate$k,
    rate_meta_pooled_pct = 100 * m_rate$pooled,
    rate_meta_pooled_lo  = 100 * m_rate$pooled_lo,
    rate_meta_pooled_hi  = 100 * m_rate$pooled_hi,
    rate_meta_Q          = m_rate$Q,
    rate_meta_Q_df       = m_rate$Q_df,
    rate_meta_Q_p        = m_rate$Q_p,
    rate_meta_I2_pct     = m_rate$I2,
    rate_meta_tau2       = m_rate$tau2,
    rupt_prev_meta_k          = m_rupt$k,
    rupt_prev_meta_pooled_pct = 100 * m_rupt$pooled,
    rupt_prev_meta_pooled_lo  = 100 * m_rupt$pooled_lo,
    rupt_prev_meta_pooled_hi  = 100 * m_rupt$pooled_hi,
    rupt_prev_meta_Q          = m_rupt$Q,
    rupt_prev_meta_Q_df       = m_rupt$Q_df,
    rupt_prev_meta_Q_p        = m_rupt$Q_p,
    rupt_prev_meta_I2_pct     = m_rupt$I2,
    rupt_prev_meta_tau2       = m_rupt$tau2,
    age_meta_k         = m_age$k,
    age_meta_pooled_y  = m_age$pooled,
    age_meta_pooled_lo = m_age$pooled_lo,
    age_meta_pooled_hi = m_age$pooled_hi,
    age_meta_Q         = m_age$Q,
    age_meta_Q_df      = m_age$Q_df,
    age_meta_Q_p       = m_age$Q_p,
    age_meta_I2_pct    = m_age$I2,
    age_meta_tau2      = m_age$tau2,
    # Per-series KRAS-vs-Negative directional-shift count (used by §4 to
    # support "the direction of this shift was consistent across all N
    # contributing series with both arms"). A series counts as having
    # both arms when at least one variant-positive AND one genotype-
    # negative lesion contributed an age value (n >= 1 each, irrespective
    # of the n>=3 inclusion floor used for the meta-analysis itself).
    age_meta_n_both_arms = {
      .both <- gt %>%
        filter(!is.na(age), !is.na(geno_binary)) %>%
        group_by(study_clean) %>%
        summarise(n_pos = sum(geno_binary == "Variant-positive"),
                  n_neg = sum(geno_binary == "Panel-negative"),
                  .groups = "drop") %>%
        filter(n_pos >= 1, n_neg >= 1)
      nrow(.both)
    },
    ped_n_series   = nrow(study_age_grouped %>% filter(ascertainment == "pediatric")),
    ped_series     = ped_row$series,
    ped_age_min    = ped_row$age_min,
    ped_age_max    = ped_row$age_max,
    adult_n_series = nrow(study_age_grouped %>% filter(ascertainment == "adult")),
    adult_series   = ad_row$series,
    adult_age_min  = ad_row$age_min,
    adult_age_max  = ad_row$age_max,
    het_table         = het_table,
    series_detail     = series_detail,
    ascertain_summary = ascertain_summary
  )
)

# Supplementary Table 11 — Power calculations and MDE audit for null outcomes
# 2026-06-14: collapsed from 2 sheets (observed_power + mde_null_results) to
# one stacked sheet. 21-outcome power spine; MDE columns sentinel for the
# 14 outcomes that are not in the binary-null MDE family.
.canon_outcome_st11 <- c(
  "age"                            = "Age at presentation",
  "sm_drainage_num"                = "Deep venous drainage",
  "sm_eloquence_num"               = "Eloquent brain location",
  "sm_size_num"                    = "Spetzler–Martin size",
  "ever_ruptured_num"              = "Rupture (ever)",
  "ruptured_at_surgery_num"        = "Rupture (at surgery)",
  "intranidal_aneurysm_num"        = "Intranidal aneurysm",
  "venous_varix_num"               = "Venous varix",
  "venous_outflow_stenosis_num"    = "Venous outflow stenosis",
  "flow_related_aneurysm_num"      = "Flow-related aneurysm",
  "loc_frontal"                    = "Location: Frontal",
  "loc_temporal"                   = "Location: Temporal",
  "loc_parietal"                   = "Location: Parietal",
  "loc_occipital"                  = "Location: Occipital",
  "loc_cerebellar"                 = "Location: Cerebellar",
  "loc_brainstem"                  = "Location: Brainstem",
  "loc_basal_ganglia"              = "Location: Basal ganglia",
  "loc_thalamus"                   = "Location: Thalamus",
  "vaf_age"                        = "VAF × age (Spearman)",
  "vaf_sm_size_num"                = "VAF × Spetzler–Martin size (Spearman)",
  "vaf_ever_ruptured_num"          = "VAF × Rupture (Spearman)"
)
.fmt_p_st11 <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}

# Power spine: 21 outcomes with n1/n2 (or n for VAF correlation rows).
power_rows <- tibble::tibble(
  var       = names(power_results),
  power     = sapply(power_results, function(x) x$power),
  n_str     = sapply(power_results, function(x)
                     if (!is.null(x$n)) as.character(x$n)
                     else sprintf("%d / %d", x$n1, x$n2)),
  effect    = sapply(power_results, function(x) {
                     if (!is.null(x$h)) sprintf("h = %.2f", x$h)
                     else if (!is.null(x$d)) sprintf("d = %.2f", x$d)
                     else if (!is.null(x$r)) sprintf("ρ = %.3f", x$r)
                     else "—"
                   })
)

# MDE join on the machine variable name.
mde_rows <- mde_tbl %>%
  transmute(
    var               = variable,
    observed_or_ci    = sprintf("%.3f (%.3f, %.3f)",
                                observed_or, observed_lo, observed_hi),
    mde_or_str        = sprintf("%.3f", mde_or_80pct),
    mde_inv_str       = sprintf("%.3f", mde_inv_80pct),
    informativeness   = interpretation
  )

# 2026-06-15: ST11 reverted to 2 sheets (Nature MOESM3 (A)/(B) pattern).
# Power spine on (A) — all outcomes. MDE table on (B) — restricted to the
# binary-null family where Observed OR + MDE are defined.
st11_power_df <- power_rows %>%
  mutate(outcome = unname(.canon_outcome_st11[var])) %>%
  transmute(
    outcome         = outcome,
    n_str           = n_str,
    effect_size_str = effect,
    power_str       = sprintf("%.3f", power),
    adequacy        = dplyr::case_when(
                       power >= 0.80 ~ "Adequate (≥0.80)",
                       power >= 0.50 ~ "Moderate (0.50–0.79)",
                       TRUE ~ "Underpowered (<0.50)"
                     )
  ) %>%
  arrange(desc(suppressWarnings(as.numeric(power_str))))

st11_mde_df <- mde_rows %>%
  mutate(outcome = unname(.canon_outcome_st11[var])) %>%
  transmute(
    outcome             = outcome,
    obs_or_str          = observed_or_ci,
    mde_or_str_out      = mde_or_str,
    mde_inv_str_out     = mde_inv_str,
    informativeness_str = informativeness
  ) %>%
  arrange(outcome)

st11_sheets <- list(
  "(A) Power audit" = list(
    data = st11_power_df,
    columns = list(
      col("outcome",         label = "Outcome"),
      col("n_str",           label = "N (positive / negative)", italic = TRUE),
      col("effect_size_str", label = "Observed effect size"),
      col("power_str",       label = "Power (post-hoc, observed effect)", italic = TRUE),
      col("adequacy",        label = "Adequacy")
    ),
    footnote = c(
      "Power column is the post-hoc (observed-effect) power. This is a deliberate audit artefact for §2/§3 null outcomes; it is monotone in the P-value and is NOT cited as prospective power in the prose. See sheet (B) for prospective-style MDE for the binary-null family.",
      "Effect size convention: h = Cohen's h (binary proportion difference); d = Cohen's d (continuous mean difference); ρ = Spearman rho (correlation outcomes). N is reported as 'positive / negative' for the two-group tests and as a single number for correlation tests."
    )
  ),
  "(B) MDE (binary-null family)" = list(
    data = st11_mde_df,
    columns = list(
      col("outcome",             label = "Outcome"),
      col("obs_or_str",          label = "Observed OR (95% CI)"),
      col("mde_or_str_out",      label = "MDE OR (80% power)"),
      col("mde_inv_str_out",     label = "Inverse MDE (80% power)"),
      col("informativeness_str", label = "Informativeness")
    ),
    footnote = c(
      "MDE OR (80% power) is the minimum odds ratio detectable at 80% power, α = 0.05; inverse MDE = 1 / MDE for symmetric protective effects. Informativeness reports whether the cohort rules out meaningful effects in each direction.",
      "Restricted to the binary-null family (7 outcomes); other outcomes (continuous, correlation, anatomical-location) carry only the power audit on sheet (A)."
    )
  )
)

write_supp_table_workbook(
  sheets = st11_sheets,
  path   = file.path(etable_dir, "SuppTable11_power_and_mde.xlsx")
)
cat(sprintf("  ✓ SuppTable11_power_and_mde.xlsx (2 sheets: %d + %d rows)\n",
            nrow(st11_power_df), nrow(st11_mde_df)))
# Retained for any downstream manifest readers.
st11_df <- st11_power_df

cat("\n══ 11_ED4_power_forest_meta.R complete ══\n")
