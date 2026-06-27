# 21_F3_rupture_score_panels.R — Fig 2 panels f + g (bedside score).
#
# Builds the two main-text panels that translate the Cox PH genotype +
# angioarchitecture model (anchor configuration:
# factor × 2) into a clinical risk stratification:
#
#   Panel f (km_by_score)            — KM rupture-free survival by integer
#                                       score (0/1/2), all eligible patients.
#   Panel g (rupture_lookup_heatmap) — heatmap of cumulative rupture
#                                       probability at reference ages, with
#                                       2.5/97.5 percentile bootstrap CIs
#                                       inside each tile.
#
# Producer also writes the fig3_score stats fragment via
# write_stats_section() into the stats manifest, so every cited number can be
# inline-r every cited number without re-running R.
#
# The equation card and the apparent/corrected/holdout C-index validation
# panel live in the supporting Extended Data figure; see
# 22_ST_rupture_score_card.R.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here)
  library(caret); library(survival); library(survminer); library(patchwork)
})

source(here("analysis", "helper_scripts", "utils.R"))  # provides MASTER_SEED
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

SCALE_FACTOR <- 2L
B_BOOT       <- 500L
REF_AGES     <- c(5, 10, 15, 20, 25, 30, 40)

set.seed(MASTER_SEED)

fig3_dir   <- here("results", "Figure2")
panel_f_dir <- panel_slot_dir("km_by_score")
panel_g_dir <- panel_slot_dir("rupture_lookup_heatmap")
for (d in c(panel_f_dir, panel_g_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
save_panel <- function(dir, name, plot, w, h)
  save_panel_impl(dir, name, plot, w, h, device = "cairo")

# ── Cohort + score derivation (fixed seed → identical to script 25 anchor) ──
df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(mutation_positive)) %>%
  mutate(
    geno_group = case_when(
      mutation == "KRAS G12D"                                    ~ "G12D",
      mutation == "KRAS G12V"                                    ~ "G12V",
      mutation_gene == "KRAS" & !(mutation %in% c("KRAS G12D", "KRAS G12V")) ~ "Other",
      mutation_gene == "BRAF"                                    ~ "Other",
      mutation_positive == FALSE                                 ~ "Negative",
      TRUE                                                        ~ NA_character_
    ),
    geno_G12D  = as.integer(geno_group == "G12D"),
    sex_female = as.integer(sex_f == "Female"),
    cox_time   = age,
    cox_event  = dplyr::case_when(
      rupture_category == "Ruptured at presentation" ~ 1L,
      rupture_category == "Never ruptured"           ~ 0L,
      TRUE                                           ~ NA_integer_
    )
  ) %>%
  filter(!is.na(cox_event), !is.na(cox_time))

ANCHOR_FEATURES <- c("geno_G12D", "sm_drainage_num", "sm_size_num")
coh_all <- df %>%
  select(all_of(c("cox_time", "cox_event", ANCHOR_FEATURES, "study_clean"))) %>%
  drop_na()
cat(sprintf("Eligible n = %d | events = %d\n",
            nrow(coh_all), sum(coh_all$cox_event)))

# Anchor Cox on all eligible data
anchor_form <- as.formula(
  paste("survival::Surv(cox_time, cox_event) ~",
        paste(sprintf("`%s`", ANCHOR_FEATURES), collapse = " + ")))
mod_all <- survival::coxph(anchor_form, data = coh_all, x = TRUE, y = TRUE)
betas   <- coef(mod_all)
ci_all  <- summary(mod_all)$conf.int
pvals   <- summary(mod_all)$coefficients[, "Pr(>|z|)"]

# Schoenfeld-residual proportional-hazards check for the anchor model.
# Per the 2026-04-29 audit, the rupture-score Cox shows a mild PH
# violation (typically driven by KRAS G12D and SM size, biologically
# expected because G12D's effect concentrates at younger ages). Expose
# the global + per-covariate P so the prose / Methods can cite them.
.zph_anchor          <- survival::cox.zph(mod_all)
score_zph_global_p   <- unname(.zph_anchor$table["GLOBAL", "p"])
.zph_anchor_per_cov  <- .zph_anchor$table[setdiff(rownames(.zph_anchor$table),
                                                  "GLOBAL"), "p", drop = TRUE]
score_zph_g12d_p     <- unname(.zph_anchor_per_cov["geno_G12D"])
score_zph_drainage_p <- unname(.zph_anchor_per_cov["sm_drainage_num"])
score_zph_size_p     <- unname(.zph_anchor_per_cov["sm_size_num"])
cat(sprintf("  Schoenfeld global p (anchor) = %.3g\n", score_zph_global_p))

# Integer card (scale × round)
integer_pts <- as.integer(round(betas[ANCHOR_FEATURES] * SCALE_FACTOR))
names(integer_pts) <- ANCHOR_FEATURES
nonzero <- integer_pts[integer_pts != 0]
score_fn <- function(data) {
  as.integer(rowSums(
    sapply(names(nonzero), function(v) data[[v]] * nonzero[v])
  ))
}
coh_all$int_score <- score_fn(coh_all)
score_levels <- sort(unique(coh_all$int_score))
cat("Integer score distribution:\n")
print(table(coh_all$int_score, coh_all$cox_event,
            dnn = c("score", "event")))

# Refit Cox on integer score (factor) for risk lookup
mod_card <- survival::coxph(
  survival::Surv(cox_time, cox_event) ~ factor(int_score),
  data = coh_all, x = TRUE, y = TRUE
)

# Stratified holdout (matches script 25 setup)
stratified_split_cox <- function(data, p = 0.8, strat = "study_clean",
                                 min_cell = 2L, seed = MASTER_SEED) {
  strat_vec <- paste(data[[strat]], data$cox_event, sep = "|")
  tab <- table(strat_vec)
  tiny <- names(tab)[tab < min_cell]
  if (length(tiny)) {
    strat_vec[strat_vec %in% tiny] <- paste0("small_studies|",
                                             data$cox_event[strat_vec %in% tiny])
  }
  set.seed(seed)
  idx <- caret::createDataPartition(factor(strat_vec), p = p,
                                     list = FALSE)[, 1]
  list(train = data[idx, ], test = data[-idx, ])
}
sp <- stratified_split_cox(coh_all)

c_index <- function(time, event, lp) {
  suppressWarnings(survival::concordance(
    survival::Surv(time, event) ~ lp, reverse = TRUE)$concordance)
}
sp$test$int_score <- score_fn(sp$test)
c_holdout_lp  <- c_index(sp$test$cox_time, sp$test$cox_event,
                          predict(mod_all, newdata = sp$test, type = "lp"))
c_holdout_int <- c_index(sp$test$cox_time, sp$test$cox_event,
                          as.numeric(sp$test$int_score))

# ── Bootstrap (Harrell) for apparent/corrected C + lookup-cell CIs ──────────
apparent_c <- c_index(coh_all$cox_time, coh_all$cox_event,
                      predict(mod_all, type = "lp"))
set.seed(MASTER_SEED)
n_all <- nrow(coh_all)

boot_run <- function() {
  idx <- sample.int(n_all, replace = TRUE)
  D_b <- coh_all[idx, ]
  mod_b <- tryCatch(survival::coxph(anchor_form, data = D_b),
                    error = function(e) NULL)
  if (is.null(mod_b)) return(NULL)
  c_train <- c_index(D_b$cox_time, D_b$cox_event,
                     predict(mod_b, newdata = D_b, type = "lp"))
  c_orig  <- c_index(coh_all$cox_time, coh_all$cox_event,
                     predict(mod_b, newdata = coh_all, type = "lp"))
  mod_card_b <- tryCatch(
    survival::coxph(survival::Surv(cox_time, cox_event) ~ factor(int_score),
                    data = D_b),
    error = function(e) NULL)
  risk_b <- matrix(NA_real_, nrow = length(score_levels),
                   ncol = length(REF_AGES))
  if (!is.null(mod_card_b)) {
    boot_levels <- intersect(score_levels, unique(D_b$int_score))
    sf_b <- tryCatch(
      survival::survfit(mod_card_b,
                        newdata = data.frame(int_score = boot_levels)),
      error = function(e) NULL)
    if (!is.null(sf_b)) {
      for (j in seq_along(REF_AGES)) {
        i <- suppressWarnings(max(which(sf_b$time <= REF_AGES[j]),
                                  na.rm = TRUE))
        if (is.finite(i) && i >= 1) {
          row_idx <- match(boot_levels, score_levels)
          risk_b[row_idx, j] <- 1 - sf_b$surv[i, ]
        }
      }
    }
  }
  list(c_train = c_train, c_orig = c_orig, risk = risk_b)
}

.boot_cache_path <- here("results", "stats", "_boot_cache_29_age_score.rds")
.boot_cache_key  <- paste(
  tools::md5sum(here("data", "processed", "bAVM_analysis_ready.rds")),
  tools::md5sum(here("analysis", "01_main_analysis", "21_F3_rupture_score_panels.R")),
  B_BOOT, MASTER_SEED, sep = "|")
if (file.exists(.boot_cache_path) &&
    identical(readRDS(.boot_cache_path)$key, .boot_cache_key)) {
  boot <- readRDS(.boot_cache_path)$boot
  cat(sprintf("Bootstrap cache hit — %d replicates loaded (skipping recompute)\n",
              length(boot)))
} else {
  cat(sprintf("Running B = %d bootstrap replicates...\n", B_BOOT))
  boot <- lapply(seq_len(B_BOOT), function(b) boot_run())
  boot <- boot[!vapply(boot, is.null, logical(1))]
  saveRDS(list(key = .boot_cache_key, boot = boot), .boot_cache_path)
}
boot <- boot[!vapply(boot, is.null, logical(1))]
c_train_b <- vapply(boot, `[[`, numeric(1), "c_train")
c_orig_b  <- vapply(boot, `[[`, numeric(1), "c_orig")
optimism  <- mean(c_train_b - c_orig_b)
opt_ci    <- quantile(c_train_b - c_orig_b, c(0.025, 0.975), na.rm = TRUE)
corrected_c <- apparent_c - optimism

risk_stack <- array(unlist(lapply(boot, `[[`, "risk")),
                    dim = c(length(score_levels), length(REF_AGES),
                            length(boot)))
risk_mean <- apply(risk_stack, c(1, 2), mean, na.rm = TRUE)
risk_lo   <- apply(risk_stack, c(1, 2), quantile, 0.025, na.rm = TRUE)
risk_hi   <- apply(risk_stack, c(1, 2), quantile, 0.975, na.rm = TRUE)

# AJK #C38: youngest REF_AGES point where the score 0 vs. score 2 bootstrap
# 95% CIs do not overlap. Drives the abstract "(age 15)" parenthetical and
# the §4 "from age 15 onward" claim. If the bootstrap separation shifts on
# regen, the prose updates automatically.
.score_separation_age_y <- {
  i0 <- match(0, score_levels); i2 <- match(2, score_levels)
  sep <- risk_lo[i2, ] > risk_hi[i0, ]
  if (any(sep, na.rm = TRUE)) REF_AGES[min(which(sep))] else NA_integer_
}

# Point-estimate lookup (from full-cohort survfit on integer-score Cox)
sf_pt <- survival::survfit(mod_card,
                           newdata = data.frame(int_score = score_levels))
risk_pt <- sapply(REF_AGES, function(a) {
  i <- max(which(sf_pt$time <= a), na.rm = TRUE)
  if (!is.finite(i) || i < 1) rep(0, ncol(sf_pt$surv)) else
    1 - sf_pt$surv[i, ]
})
if (is.null(dim(risk_pt))) risk_pt <- matrix(risk_pt, nrow = 1)
colnames(risk_pt) <- paste0("age_", REF_AGES)

# ═════════════════════════════════════════════════════════════════════════════
# PANEL F — KM rupture-free survival by integer score
# ═════════════════════════════════════════════════════════════════════════════
# Score-tier KM palette: sourced from PAL_SCORE in utils.R (Tier C — three
# distinct hues outside Tier A variant colours and Tier B Mut+ green; score
# 0 reuses grey as the null tier, score 1 = teal, score 2 = purple).
SCORE_PAL <- PAL_SCORE[as.character(score_levels)]

km_fit <- survival::survfit(
  survival::Surv(cox_time, cox_event) ~ factor(int_score), data = coh_all)
pF <- survminer::ggsurvplot(
  legend = "bottom", km_fit, data = coh_all,
  conf.int = TRUE, conf.int.alpha = 0.2,
  # v6.7 (2026-05-20): linewidth + ribbon alpha unified across Fig 1
  # KM/density/regression line geoms.
  size     = 0.5,
  palette = unname(SCORE_PAL),
  risk.table = TRUE, risk.table.height = 0.25,
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
  xlab = "Age (years)", ylab = "Proportion rupture-free",
  legend.title = "Score",
  legend.labs = paste0("Score ", score_levels,
                       "  (n = ", as.integer(table(coh_all$int_score)), ")"),
  break.time.by = 10, xlim = c(0, KM_AGE_XLIM_MAX),  # see KM_AGE_XLIM_MAX in utils.R
  censor.size = 1.5,  # smaller censor ticks; alpha applied below
  ggtheme = theme_avm(),
  tables.theme = theme_cleantable()
)
# Dim censor marks (ggsurvplot exposes no censor.alpha) — walk $plot$layers
# and reduce alpha on the GeomPoint layer. Pattern matches 09_F1_km_age.R.
pF$plot$layers <- lapply(pF$plot$layers, function(l) {
  if (inherits(l$geom, "GeomPoint")) l$aes_params$alpha <- 0.4
  l
})
# save_km_panel splits the ggsurvplot block into curve + table RDS so the
# native-size composer (compose_figure.R) can re-theme each at print
# footprint. Standalone PDF/PNG still saved at 6.5x6 for review.
save_km_panel(panel_f_dir, "km_by_score", pF, w = 10, h = 8)

# Reverse-Kaplan–Meier (score-stratified) retired 2026-05-17 — see
# `09_F1_km_age.R` and `results/stats/reverse_censoring_audit.md` for
# the audit. Reverse-KM omnibus log-rank across score tiers gave
# P = 0.28 (balanced follow-up), so the panel was uninformative.

# ═════════════════════════════════════════════════════════════════════════════
# PANEL G — Score × age cumulative rupture probability heatmap with CIs
# ═════════════════════════════════════════════════════════════════════════════
grid_long <- tibble::tibble(
  score = rep(score_levels, times = length(REF_AGES)),
  ref_age = rep(REF_AGES, each = length(score_levels)),
  pct = c(round(100 * risk_mean, 1)),
  lo  = c(round(100 * risk_lo, 1)),
  hi  = c(round(100 * risk_hi, 1))
)
# Contrast-aware text colour: tiles whose fill is darker than the
# midpoint of the white→purple ramp get white text, lighter tiles get
# black. Threshold tuned to the gradient's perceived midpoint (~55% of
# the max).
.fill_max  <- max(grid_long$hi, na.rm = TRUE)
.txt_thresh <- 0.55 * .fill_max
grid_long$text_col <- ifelse(grid_long$pct >= .txt_thresh, "white", "black")

pG <- ggplot(grid_long,
             aes(x = factor(ref_age), y = factor(score), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  # 2026-05-19: strip the (lo–hi) CI line per design — just keep the
  # point estimate. CIs live in Supplementary Table N (rupture_lookup
  # by score x age) so the figure doesn't carry them twice.
  geom_text(aes(label = sprintf("%.0f%%", pct),
                colour = text_col),
            # v6.2 (2026-05-21): size 4.0 mm (~11 pt) was way over Nature's
            # 5-7 pt ceiling at native cell footprint. 2.0 mm ~ 5.7 pt.
            size = NM$text$body_mm, lineheight = 0.95, fontface = "bold") +
  # 2026-05-19 (revert): keep the white -> score-2 purple ramp so the
  # heatmap stays visually anchored to panel H (score 2 = highest
  # risk = darkest purple cell). The colorbar is rendered horizontal
  # with only the min/max ticks so it consumes less horizontal real
  # estate at composite scale.
  scale_fill_gradient(
    low    = "white",
    high   = unname(PAL_SCORE[["2"]]),
    name   = "Rupture risk (%)",
    # v6.2 (2026-05-21): cap at 80 instead of .fill_max (which reaches
    # ~84). Data tops out around 74; capping at 80 trims the empty
    # tail from the colorbar and gives a cleaner round scale.
    limits = c(0, 80),
    breaks = c(0, 80),
    labels = c("0", "80"),
    guide  = guide_colorbar(
      direction        = "horizontal",
      barwidth         = unit(4.5, "cm"),
      barheight        = unit(0.35, "cm"),
      title.position   = "top",
      title.hjust      = 0.5,
      ticks            = FALSE,
      frame.colour     = "grey30",
      frame.linewidth  = 0.3
    )) +
  scale_colour_identity() +
  scale_y_discrete(limits = as.character(rev(score_levels))) +
  labs(x = "Age (years)", y = "Integer score") +
  theme_avm() +
  theme(panel.grid       = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_text(size = rel(0.9)),
        legend.text      = element_text(size = rel(0.85)))
  # Font sizes inherit from theme_avm() (utils.R typography contract).
save_panel(panel_g_dir, "rupture_lookup_heatmap", pG, w = 9, h = 4.5)

# ── Stats fragment ──────────────────────────────────────────────────────────
pull_age <- function(score, age, mat) round(100 * mat[
  match(score, score_levels),
  match(age, REF_AGES)], 1)

write_stats_section("fig3_score", list(
  n_eligible = nrow(coh_all),
  n_events   = sum(coh_all$cox_event),
  n_train    = nrow(sp$train),
  n_test     = nrow(sp$test),
  anchor_features = ANCHOR_FEATURES,
  hr_g12d        = unname(ci_all["geno_G12D",       "exp(coef)"]),
  hr_g12d_lo     = unname(ci_all["geno_G12D",       "lower .95"]),
  hr_g12d_hi     = unname(ci_all["geno_G12D",       "upper .95"]),
  hr_g12d_p      = unname(pvals["geno_G12D"]),
  hr_drainage    = unname(ci_all["sm_drainage_num", "exp(coef)"]),
  hr_drainage_lo = unname(ci_all["sm_drainage_num", "lower .95"]),
  hr_drainage_hi = unname(ci_all["sm_drainage_num", "upper .95"]),
  hr_drainage_p  = unname(pvals["sm_drainage_num"]),
  hr_size        = unname(ci_all["sm_size_num",     "exp(coef)"]),
  hr_size_lo     = unname(ci_all["sm_size_num",     "lower .95"]),
  hr_size_hi     = unname(ci_all["sm_size_num",     "upper .95"]),
  hr_size_p      = unname(pvals["sm_size_num"]),
  zph_global_p   = score_zph_global_p,
  zph_g12d_p     = score_zph_g12d_p,
  zph_drainage_p = score_zph_drainage_p,
  zph_size_p     = score_zph_size_p,
  scale_factor   = SCALE_FACTOR,
  pts_g12d       = unname(integer_pts["geno_G12D"]),
  pts_drainage   = unname(integer_pts["sm_drainage_num"]),
  pts_size       = unname(integer_pts["sm_size_num"]),
  score_min      = min(score_levels),
  score_max      = max(score_levels),
  score1_hr      = unname(summary(mod_card)$conf.int["factor(int_score)1", "exp(coef)"]),
  score1_hr_lo   = unname(summary(mod_card)$conf.int["factor(int_score)1", "lower .95"]),
  score1_hr_hi   = unname(summary(mod_card)$conf.int["factor(int_score)1", "upper .95"]),
  score1_p       = unname(summary(mod_card)$coefficients["factor(int_score)1", "Pr(>|z|)"]),
  score2_hr      = unname(summary(mod_card)$conf.int["factor(int_score)2", "exp(coef)"]),
  score2_hr_lo   = unname(summary(mod_card)$conf.int["factor(int_score)2", "lower .95"]),
  score2_hr_hi   = unname(summary(mod_card)$conf.int["factor(int_score)2", "upper .95"]),
  score2_p       = unname(summary(mod_card)$coefficients["factor(int_score)2", "Pr(>|z|)"]),
  c_holdout_lp     = c_holdout_lp,
  c_holdout_int    = c_holdout_int,
  c_holdout_delta  = c_holdout_lp - c_holdout_int,
  boot_B           = length(boot),
  c_apparent       = apparent_c,
  c_corrected      = corrected_c,
  optimism_mean    = optimism,
  optimism_ci_lo   = unname(opt_ci[1]),
  optimism_ci_hi   = unname(opt_ci[2]),
  lookup_grid_point = as.data.frame(round(100 * risk_pt, 1)) |>
                      tibble::rownames_to_column("score_idx") |>
                      dplyr::mutate(score = score_levels) |>
                      dplyr::select(-score_idx),
  lookup_grid_ci    = grid_long,
  score_separation_age_y = .score_separation_age_y,
  risk_age20_score0_pct = pull_age(0, 20, risk_mean),
  risk_age20_score0_lo  = pull_age(0, 20, risk_lo),
  risk_age20_score0_hi  = pull_age(0, 20, risk_hi),
  risk_age20_score2_pct = pull_age(2, 20, risk_mean),
  risk_age20_score2_lo  = pull_age(2, 20, risk_lo),
  risk_age20_score2_hi  = pull_age(2, 20, risk_hi),
  risk_age10_score0_pct = pull_age(0, 10, risk_mean),
  risk_age10_score2_pct = pull_age(2, 10, risk_mean),
  risk_age30_score0_pct = pull_age(0, 30, risk_mean),
  risk_age30_score2_pct = pull_age(2, 30, risk_mean),
  n_score0       = sum(coh_all$int_score == 0),
  n_score1       = sum(coh_all$int_score == 1),
  n_score2       = sum(coh_all$int_score == 2),
  events_score0  = sum(coh_all$cox_event[coh_all$int_score == 0]),
  events_score1  = sum(coh_all$cox_event[coh_all$int_score == 1]),
  events_score2  = sum(coh_all$cox_event[coh_all$int_score == 2])
))

cat(sprintf("\n✓ Wrote panels to %s\n", fig3_dir))
cat(sprintf("✓ Apparent C = %.3f | optimism = %.3f | corrected C = %.3f\n",
            apparent_c, optimism, corrected_c))
