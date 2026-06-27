# 22_ST_rupture_score_card.R — bedside rupture-age score Supplementary Table.
#
# 2026-04-26: the former ed_rupture_score_card ED figure (panels b/c =
# Cox HR forest + bootstrap validation lollipop) was retired. The
# numerical content is small (3 anchor-model HR rows + 6 validation
# scalars) and read more cleanly as a Supplementary Table than as a
# figure. This script now writes only that supplementary table:
#
#   results/SupplementaryTables/SuppTable14_score_card.xlsx
#     Sheet 1 "Anchor model"      — Cox HR for the three anchor features
#                                   (G12D, deep drainage, SM size) plus
#                                   the integer points each contributes
#                                   to the bedside score under the
#                                   scale-factor-2 rounding.
#     Sheet 2 "Bootstrap validation" — apparent / optimism-corrected /
#                                       holdout C-indices and the
#                                       bootstrap optimism distribution
#                                       summary.
#
# Producer 21_F3_rupture_score_panels.R remains the source of truth for the
# `fig3_score` stats manifest fragment that the §4 prose consumes; this
# script does NOT write a stats fragment because every value it computes
# is also computed there with the same seed.
#
# Score equation card (former panel A) was retired the same day; the
# equation is fully described in the §4 prose.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(here)
  library(caret); library(survival); library(tibble)
})

source(here("analysis", "helper_scripts", "utils.R"))  # provides MASTER_SEED

SCALE_FACTOR <- 2L
B_BOOT       <- 500L

set.seed(MASTER_SEED)

supp_dir <- here("results", "SupplementaryTables")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

# ── Cohort + anchor refit (same seed as 21_F3_rupture_score_panels.R) ──────────────
df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(mutation_positive)) %>%
  mutate(
    geno_G12D  = as.integer(mutation == "KRAS G12D"),
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

anchor_form <- as.formula(
  paste("survival::Surv(cox_time, cox_event) ~",
        paste(sprintf("`%s`", ANCHOR_FEATURES), collapse = " + ")))
mod_all <- survival::coxph(anchor_form, data = coh_all, x = TRUE, y = TRUE)
betas   <- coef(mod_all)
ci_all  <- summary(mod_all)$conf.int
pvals   <- summary(mod_all)$coefficients[, "Pr(>|z|)"]

integer_pts <- as.integer(round(betas * SCALE_FACTOR))
names(integer_pts) <- ANCHOR_FEATURES

# Stratified holdout split (matches script 25/26/29)
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
c_holdout_lp <- c_index(sp$test$cox_time, sp$test$cox_event,
                         predict(mod_all, newdata = sp$test, type = "lp"))

# ── Bootstrap (Harrell) — same logic as scripts 26 / 29 ────────────────────
apparent_c <- c_index(coh_all$cox_time, coh_all$cox_event,
                      predict(mod_all, type = "lp"))
set.seed(MASTER_SEED)
n_all <- nrow(coh_all)
boot_c <- function() {
  idx <- sample.int(n_all, replace = TRUE)
  D_b <- coh_all[idx, ]
  mod_b <- tryCatch(survival::coxph(anchor_form, data = D_b),
                    error = function(e) NULL)
  if (is.null(mod_b)) return(NULL)
  list(
    train = c_index(D_b$cox_time, D_b$cox_event,
                    predict(mod_b, newdata = D_b, type = "lp")),
    orig  = c_index(coh_all$cox_time, coh_all$cox_event,
                    predict(mod_b, newdata = coh_all, type = "lp"))
  )
}
.boot_cache_path <- here("results", "stats", "_boot_cache_30_score_card.rds")
.boot_cache_key  <- paste(
  tools::md5sum(here("data", "processed", "bAVM_analysis_ready.rds")),
  tools::md5sum(here("analysis", "01_main_analysis", "22_ST_rupture_score_card.R")),
  B_BOOT, MASTER_SEED, sep = "|")
if (file.exists(.boot_cache_path) &&
    identical(readRDS(.boot_cache_path)$key, .boot_cache_key)) {
  bres <- readRDS(.boot_cache_path)$bres
  cat(sprintf("Bootstrap cache hit — %d replicates loaded (skipping recompute)\n",
              length(bres)))
} else {
  cat(sprintf("Running B = %d bootstrap replicates...\n", B_BOOT))
  bres <- lapply(seq_len(B_BOOT), function(b) boot_c())
  bres <- bres[!vapply(bres, is.null, logical(1))]
  saveRDS(list(key = .boot_cache_key, bres = bres), .boot_cache_path)
}
bres <- bres[!vapply(bres, is.null, logical(1))]
c_train_b <- vapply(bres, `[[`, numeric(1), "train")
c_orig_b  <- vapply(bres, `[[`, numeric(1), "orig")
optimism  <- mean(c_train_b - c_orig_b)
opt_ci    <- quantile(c_train_b - c_orig_b, c(0.025, 0.975), na.rm = TRUE)
corrected_c <- apparent_c - optimism

# ── Sheet 1: Anchor model HRs ───────────────────────────────────────────────
anchor_sheet <- tibble(
  feature           = c("KRAS G12D",
                        "Deep venous drainage",
                        "Spetzler–Martin size (per grade)"),
  hazard_ratio      = unname(ci_all[ANCHOR_FEATURES, "exp(coef)"]),
  ci_low_95         = unname(ci_all[ANCHOR_FEATURES, "lower .95"]),
  ci_high_95        = unname(ci_all[ANCHOR_FEATURES, "upper .95"]),
  p_value           = unname(pvals[ANCHOR_FEATURES]),
  integer_points    = unname(integer_pts[ANCHOR_FEATURES]),
  scale_factor      = SCALE_FACTOR,
  n                 = nrow(coh_all),
  events            = sum(coh_all$cox_event)
)

# ── Sheet 2: Bootstrap validation ───────────────────────────────────────────
validation_sheet <- tibble(
  metric    = c("Apparent C-index (full cohort)",
                "Mean optimism (Harrell bootstrap)",
                "Optimism-corrected C-index",
                "Holdout C-index (stratified 20% split)"),
  value     = c(apparent_c, optimism, corrected_c, c_holdout_lp),
  ci_low_95 = c(NA_real_,
                unname(opt_ci[1]),
                apparent_c - unname(opt_ci[2]),
                NA_real_),
  ci_high_95 = c(NA_real_,
                 unname(opt_ci[2]),
                 apparent_c - unname(opt_ci[1]),
                 NA_real_),
  bootstrap_B = c(NA_integer_, length(bres), length(bres), NA_integer_),
  notes = c("Apparent C is uncorrected; serves only as the baseline for the optimism estimate.",
            "Mean of (train C - test C) across bootstrap replicates.",
            "Apparent C - mean optimism. CIs are translated from the bootstrap optimism quantiles.",
            "Single 80/20 stratified split (stratified by series + event).")
)

st_path <- file.path(supp_dir, "SuppTable14_score_card.xlsx")
# 2026-06-14: collapsed 2 sheets (Anchor model + Bootstrap validation) into
# one stacked sheet with a Block column. Sentinel '—' for fields not
# defined per row type (e.g. P-value on validation-metric rows; bootstrap
# CIs on the apparent/holdout rows that aren't bootstrap-derived).
source(here("analysis", "helper_scripts", "supp_table_writer.R"))

.fmt_p_st14 <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}

# 2026-06-15: ST14 back to 2 sheets (Nature MOESM3 (A)/(B) pattern).
# (A) Anchor Cox model — per-feature HR, integer points. (B) Bootstrap
# validation metrics with optimism correction.
anchor_sheet_df <- anchor_sheet %>%
  transmute(
    feature       = feature,
    n             = as.integer(n),
    events        = as.integer(events),
    estimate      = sprintf("%.3f", hazard_ratio),
    ci_str        = sprintf("(%.3f, %.3f)", ci_low_95, ci_high_95),
    p_str         = .fmt_p_st14(p_value),
    notes         = sprintf("Integer points = %d; scale factor = %d",
                            as.integer(integer_points), as.integer(scale_factor))
  )

validation_sheet_df <- validation_sheet %>%
  transmute(
    metric      = metric,
    estimate    = sprintf("%.3f", value),
    ci_str      = ifelse(is.na(ci_low_95), "—",
                          sprintf("(%.3f, %.3f)", ci_low_95, ci_high_95)),
    notes       = ifelse(is.na(bootstrap_B),
                          notes,
                          sprintf("%s (B = %d)", notes, as.integer(bootstrap_B)))
  )

st14_sheets <- list(
  "(A) Anchor model (Cox)" = list(
    data = anchor_sheet_df,
    columns = list(
      col    ("feature",  label = "Feature"),
      col_int("n",        label = "N",       italic = TRUE),
      col_int("events",   label = "Events"),
      col    ("estimate", label = "HR",      italic = TRUE),
      col    ("ci_str",   label = "95% CI"),
      col    ("p_str",    label = "P",       italic = TRUE),
      col    ("notes",    label = "Notes")
    ),
    footnote = "Three pre-specified Cox proportional-hazards covariates and their integer-point representations in the rupture score card. Estimate = hazard ratio. N = cohort size; Events = ruptures in that cohort."
  ),
  "(B) Bootstrap validation" = list(
    data = validation_sheet_df,
    columns = list(
      col("metric",   label = "Metric"),
      col("estimate", label = "C-index", italic = TRUE),
      col("ci_str",   label = "95% CI"),
      col("notes",    label = "Notes")
    ),
    footnote = c(
      "Apparent C is the in-sample C-index (no overfitting correction); Mean optimism is the average (train C − test C) across B Harrell bootstrap replicates; Optimism-corrected C = Apparent C − Mean optimism; Holdout C is from a stratified 80/20 split (stratified by series + event).",
      "'—' = the Apparent and Holdout C-index rows are point estimates without bootstrap CIs."
    )
  )
)

write_supp_table_workbook(
  sheets = st14_sheets,
  path   = st_path
)
cat(sprintf("  ✓ %s (2 sheets: %d + %d rows)\n",
            basename(st_path),
            nrow(anchor_sheet_df), nrow(validation_sheet_df)))
st14_df <- bind_rows(
  anchor_sheet_df %>% transmute(item = feature, n_str = as.character(n),
                                events_str = as.character(events),
                                estimate, ci_str, p_str, notes_str = notes),
  validation_sheet_df %>% transmute(item = metric, n_str = "—", events_str = "—",
                                    estimate, ci_str, p_str = "—",
                                    notes_str = notes)
)
cat(sprintf("\n✓ Wrote %s\n", st_path))
cat(sprintf("✓ Apparent C = %.3f | corrected C = %.3f | holdout C = %.3f\n",
            apparent_c, corrected_c, c_holdout_lp))
