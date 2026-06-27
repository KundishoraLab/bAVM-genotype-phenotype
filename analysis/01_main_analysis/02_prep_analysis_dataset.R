# 02_prep_analysis_dataset.R — Prepare analysis-ready dataset from cleaned master
#
# Input:  data/processed/bAVM_genopheno_clean.rds
# Output: data/processed/bAVM_analysis_ready.rds
#
# Steps:
#   1. Load cleaned data
#   2. Create genotype grouping variables (Tier 1 binary + Tier 2 per-variant)
#   3. Extract recurrence from CHOP notes → recurrence_num
#   4. Verify/clean growing_num
#   5. Define analysable subsets (genotyped-only filter)
#   6. Generate data completeness summary + Extended Data Fig. 1 (missingness heatmap)
#   7. Export analysis-ready dataset
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)

source(here("analysis", "helper_scripts", "utils.R"))

# ── 1. Load cleaned data ────────────────────────────────────────────────────

master <- readRDS(here("data", "processed", "bAVM_genopheno_clean.rds"))
cat(sprintf("Loaded master dataset: %d rows × %d columns\n", nrow(master), ncol(master)))

# ── Subset to AVM-event level (one row per lesion event) ───────────────────
# Per DATA_DECISIONS.md §41, the analysis-ready dataset operates at the
# AVM-event level (n=473): each lesion event (primary or recurrence) gets
# one row. The cleaner preserves all sample-level rows in
# bAVM_genopheno_clean.rds (so the duplicate_block secondary block and the
# mosaic S2 samples remain accessible for sample-level VAF analyses), but
# the analysis-ready dataset filters to `is_event_primary_row == TRUE` so
# every event is counted exactly once.
#
# Effect on the 18 May 26 cohort: 476 → 473 rows.
#  - 1 duplicate_block S2 row dropped (CHOP62's second pathology block;
#    discarded VAF preserved in pub_data_note of the surviving S1)
#  - 2 intra_avm_mosaic S2 rows dropped (AVMUAB019, CHOP44) — both rows
#    share a lesion_id, so the AVM event is counted once at the cohort
#    level; the non-primary sample-level row remains in the cleaned rds
#    if a downstream analysis needs to look at both samples.
#  - Recurrent AVM rows (16) all kept — each is its own AVM event.
if ("is_event_primary_row" %in% names(master)) {
  n_before <- nrow(master)
  master <- master %>% filter(is_event_primary_row)
  cat(sprintf("  → AVM-event subset (is_event_primary_row==TRUE): %d → %d rows (-%d)\n",
              n_before, nrow(master), n_before - nrow(master)))
} else {
  warning("master is missing `is_event_primary_row` — running at sample level. ",
          "Regenerate via 01_clean_master.R.", call. = FALSE)
}

# ── 2. Create genotype grouping variables ────────────────────────────────────

# Tier 1: Binary (variant-positive vs "Panel-negative")
# Tier 2: Per-variant (G12D, G12V, other_KRAS, BRAF, "Negative")
# Pending patients get NA for both tiers (excluded from genotype analyses)
#
# Nomenclature note (2026-05-16): geno_binary now uses "Panel-negative" as
# the factor level directly — global rename from legacy "Genotype-negative"
# applied across all producers. The Tier-2 variant_group factor keeps the
# terser "Negative" level internally; render-layer helpers (see
# 09_F1_km_age.R::display_levels) relabel it as "Panel-negative" in legends
# / risk-table strata so the user-visible token is consistent everywhere.
# "Panel-negative" emphasises that institutional samples were screened
# only on the 5-allele multiplex dPCR (KRAS G12D/G12V/G12C/G12A + BRAF
# V600E), so a negative call does not exclude variants outside the
# assayed loci.

master <- master %>%
  mutate(
    # Tier 1: Binary genotype
    geno_binary = case_when(
      is.na(mutation_positive) ~ NA_character_,
      mutation_positive == TRUE ~ "Variant-positive",
      mutation_positive == FALSE ~ "Panel-negative"
    ),
    geno_binary = factor(geno_binary, levels = c("Variant-positive", "Panel-negative")),

    # Tier 2: Per-variant genotype
    geno_variant = case_when(
      is.na(mutation_positive) ~ NA_character_,
      mut_KRAS_G12D == 1 ~ "KRAS G12D",
      mut_KRAS_G12V == 1 ~ "KRAS G12V",
      mutation_gene == "KRAS" ~ "Other KRAS",
      mutation_gene == "BRAF" ~ "BRAF",
      mutation_positive == FALSE ~ "Negative"
    ),
    geno_variant = factor(geno_variant,
      levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF", "Negative")
    ),

    # Genotype status (3-level: includes pending)
    geno_status = case_when(
      is.na(mutation_positive) ~ "Pending",
      mutation_positive == TRUE ~ "Variant-positive",
      mutation_positive == FALSE ~ "Panel-negative"
    ),
    geno_status = factor(geno_status,
      levels = c("Variant-positive", "Panel-negative", "Pending")
    )
  )

# Validate grouping
cat("\n── Tier 1 (binary) ──\n")
print(table(master$geno_binary, useNA = "ifany"))

cat("\n── Tier 2 (per-variant) ──\n")
print(table(master$geno_variant, useNA = "ifany"))

cat("\n── Genotype status (3-level) ──\n")
print(table(master$geno_status, useNA = "ifany"))

# Cross-check: Other KRAS breakdown
cat("\n── Other KRAS detail ──\n")
master %>%
  filter(geno_variant == "Other KRAS") %>%
  select(patient_id, mutation, mutation_gene) %>%
  print()

# ── 2b. Rupture variable validation and previous_rupture construction ───────
# Per Hale:
#   ruptured_at_surgery = 1 → ruptured at presentation (active hemorrhage)
#   ever_ruptured = 1 → ruptured at some point prior to or at intervention
#   If ruptured_at_surgery=1 then ever_ruptured SHOULD be 1
#   previous_rupture = ever_ruptured=1 AND ruptured_at_surgery=0

cat("\n── Rupture variable validation ──\n")

# Post-condition check: 01_clean_master.R [17b]/§42 reconciles
# ruptured_at_surgery=1 ⇒ ever_ruptured=1 (lifetime coding, Hale v0 L25), so no
# rupt_surg=1 & ever=0 row should survive to this stage.
rupture_incon <- master %>%
  filter(ruptured_at_surgery_num == 1 & ever_ruptured_num == 0)
if (nrow(rupture_incon) > 0) {
  warning(sprintf("Rupture invariant violated: %d row(s) with rupt_surg=1 & ever=0 reached 00a (expected 0 after §42). IDs: %s",
                  nrow(rupture_incon), paste(rupture_incon$patient_id, collapse = ", ")))
} else {
  cat("✓ Invariant holds: ruptured_at_surgery=1 ⇒ ever_ruptured=1 (0 residual contradictions; §42)\n")
}

# Construct previous_rupture: a documented rupture BEFORE the surgical event.
# Because ever_ruptured is lifetime/inclusive (§42), a presentation-rupture
# patient (rupt_surg=1) carries ever=1 for that *presenting* bleed — which does
# NOT establish a *prior* rupture. The data cannot distinguish "presentation
# only" from "presentation + earlier bleed" within the rupt_surg=1 group, so we
# do not assert prior rupture for them. previous_rupture=1 ONLY when the rupture
# is unambiguously before surgery: ever=1 AND rupt_surg=0 (elective resection of
# a previously-ruptured AVM). All rupt_surg=1 rows → 0 here; the presentation
# event is captured by rupture_category == "Ruptured at presentation".
master <- master %>%
  mutate(
    previous_rupture = case_when(
      is.na(ever_ruptured_num) | is.na(ruptured_at_surgery_num) ~ NA_integer_,
      ever_ruptured_num == 1 & ruptured_at_surgery_num == 0 ~ 1L,  # prior rupture, elective surgery
      TRUE ~ 0L  # never ruptured, or presentation rupture (no documented prior)
    ),
    # 3-level rupture category
    rupture_category = case_when(
      is.na(ever_ruptured_num) & is.na(ruptured_at_surgery_num) ~ NA_character_,
      ruptured_at_surgery_num == 1 ~ "Ruptured at presentation",
      previous_rupture == 1 ~ "Prior rupture (not at surgery)",
      ever_ruptured_num == 0 & ruptured_at_surgery_num == 0 ~ "Never ruptured",
      TRUE ~ NA_character_
    ),
    rupture_category = factor(rupture_category,
      levels = c("Never ruptured", "Prior rupture (not at surgery)", "Ruptured at presentation"))
  )

cat("\n── Rupture category (3-level) ──\n")
print(table(master$rupture_category, useNA = "ifany"))

cat("\n── Rupture × genotype ──\n")
print(table(master$rupture_category, master$geno_binary, useNA = "ifany"))

# ── 3. Extract recurrence from CHOP notes ───────────────────────────────────

# 13 CHOP patients mention recurrence/residual in notes (all are pending genotype)
# Create a formal binary variable for future use

master <- master %>%
  mutate(
    recurrence_num = case_when(
      grepl("recur|regrow", notes, ignore.case = TRUE) ~ 1L,
      # "residual" only when paired with AVM context, not surgical residual
      grepl("residual AVM|residual.*resect", notes, ignore.case = TRUE) ~ 1L,
      TRUE ~ NA_integer_
    )
  )

# Review what we extracted
cat("\n── Recurrence extraction ──\n")
recur_patients <- master %>% filter(recurrence_num == 1)
cat(sprintf("Patients flagged as recurrence: %d\n", nrow(recur_patients)))
cat("Study breakdown:\n")
print(table(recur_patients$study_clean))
cat("Genotype status:\n")
print(table(recur_patients$geno_status))

# ── 4. Verify growing_num ───────────────────────────────────────────────────

cat("\n── Growth variable check ──\n")
cat("growing_num distribution:\n")
print(table(master$growing_num, master$geno_status, useNA = "ifany"))

cat(sprintf("\nGenotyped patients with growth data: %d / %d\n",
  sum(!is.na(master$growing_num) & !is.na(master$mutation_positive)),
  sum(!is.na(master$mutation_positive))
))

# ── 5. Pediatric/adult age grouping ─────────────────────────────────────────

master <- master %>%
  mutate(
    age_group = case_when(
      is.na(age) ~ NA_character_,
      age < 18 ~ "Pediatric (<18)",
      TRUE ~ "Adult (>=18)"
    ),
    age_group = factor(age_group, levels = c("Pediatric (<18)", "Adult (>=18)")),

    # Relabel sample type for clarity
    sample_type_clean = case_when(
      sample_type_f == "FFPE" ~ "FFPE",
      sample_type_f == "Tissue" ~ "Fresh/Frozen",
      sample_type_f == "Literature" ~ "Literature",
      TRUE ~ NA_character_
    ),
    sample_type_clean = factor(sample_type_clean, levels = c("FFPE", "Fresh/Frozen", "Literature"))
  )

cat("\n── Age group ──\n")
print(table(master$age_group, master$geno_status, useNA = "ifany"))

# ── 6. Define genotyped subset ──────────────────────────────────────────────

genotyped <- master %>% filter(!is.na(mutation_positive))
cat(sprintf("\n── Genotyped subset: %d patients ──\n", nrow(genotyped)))

# Data completeness summary for genotyped patients
analysis_vars <- c(
  "sm_grade", "sm_size_num", "sm_eloquence_num", "sm_drainage_num",
  "ever_ruptured_num", "location_codes", "n_high_risk_num",
  "intranidal_aneurysm_num", "venous_varix_num",
  "venous_outflow_stenosis_num", "flow_related_aneurysm_num",
  "vaf_prop", "age", "sex_f", "prior_seizure_num",
  "prior_radiation_num", "prior_embolization_num",
  "growing_num", "laterality_f", "compact_nidus_num"
)

completeness <- genotyped %>%
  summarise(across(all_of(analysis_vars), ~ sum(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_available") %>%
  mutate(
    n_total = nrow(genotyped),
    pct = round(100 * n_available / n_total, 1)
  ) %>%
  arrange(desc(pct))

cat("\n── Variable completeness (genotyped, n=", nrow(genotyped), ") ──\n")
print(as.data.frame(completeness), row.names = FALSE)

# ── 7. Extended Data Fig. 1: Missingness heatmap ───────────────────────────────────────

# Build study × variable missingness matrix for genotyped patients.
#
# Two variable-level notes:
#  - `vaf_prop` is NA by construction for mutation-negative cases (no variant
#    to quantify). Reporting VAF completeness against the full genotyped
#    cohort therefore conflates biology with missingness. We restrict the VAF
#    row to the variant-positive stratum and label it accordingly.
#  - `growing_num` has systematically missing data across nearly all studies
#    (longitudinal imaging unavailable for most contributing cohorts) and is
#    excluded from this figure.
key_vars <- c(
  "sm_size_num", "sm_drainage_num", "sm_eloquence_num", "sm_grade",
  "ever_ruptured_num", "location_codes", "n_high_risk_num",
  "vaf_prop", "age", "sex_f",
  "prior_seizure_num", "prior_radiation_num", "prior_embolization_num",
  "laterality_f"
)

# Pretty labels for the heatmap
var_labels <- c(
  sm_size_num = "SM Size", sm_drainage_num = "SM Drainage",
  sm_eloquence_num = "SM Eloquence", sm_grade = "SM Grade (composite)",
  ever_ruptured_num = "Rupture", location_codes = "Location",
  n_high_risk_num = "High-Risk Features",
  vaf_prop = "VAF (variant-positive only)",
  age = "Age at Presentation", sex_f = "Sex",
  prior_seizure_num = "Seizure Hx", prior_radiation_num = "Radiation Hx",
  prior_embolization_num = "Embolization Hx",
  laterality_f = "Laterality"
)

# Per-study completeness. For each variable we use all genotyped patients
# as the denominator EXCEPT `vaf_prop`, where the denominator is restricted
# to variant-positive cases within that study.
pct_available <- function(x, restrict = NULL) {
  if (!is.null(restrict)) {
    keep <- !is.na(restrict) & restrict == 1
    x <- x[keep]
  }
  if (length(x) == 0) return(NA_real_)
  round(100 * sum(!is.na(x)) / length(x), 1)
}

missingness_matrix <- genotyped %>%
  group_by(study_clean) %>%
  summarise(
    n = n(),
    across(all_of(setdiff(key_vars, "vaf_prop")),
           ~ round(100 * sum(!is.na(.)) / n(), 1)),
    vaf_prop = pct_available(vaf_prop, restrict = mutation_positive),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(key_vars),
    names_to = "variable",
    values_to = "pct_available"
  ) %>%
  mutate(
    var_label = var_labels[variable],
    var_label = factor(var_label, levels = rev(var_labels))
  )

# Order studies by descending mean completeness across the assessed
# variables, so the heatmap reads left-to-right from the most-complete
# series to the least-complete (per AJK v0 #C97).
missingness_matrix <- missingness_matrix %>%
  mutate(study_clean = forcats::fct_reorder(
    study_clean, pct_available, .fun = mean, na.rm = TRUE, .desc = TRUE
  ))

# Create output directory
output_dir <- here("results")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Extended Data Fig. 1: Missingness heatmap
# v6.19 (2026-05-20): adopt Fig 1 polish standards — theme_nature_panel
# (7 pt body / 6 pt ticks Helvetica), blanket no-grid, geom_text dropped
# 3.0 -> 2.4 to stay inside Nature's 5-7 pt body band at the 7.20 in
# canvas. PAL_HEAT_HIGH = W_BLUE (Wong palette).
efig1 <- ggplot(missingness_matrix,
    aes(x = study_clean, y = var_label, fill = pct_available)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct_available)),
    size = NM$text$body_mm, color = "black") +
  scale_fill_gradient(
    low = PAL_HEAT_LOW, high = PAL_HEAT_HIGH, limits = c(0, 100),
    name = "% Available"
  ) +
  labs(x = "Study", y = NULL) +
  theme_nature_panel() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

efig01_dir <- file.path(output_dir, "ExtendedData", "ed_completeness")
dir.create(efig01_dir, recursive = TRUE, showWarnings = FALSE)
# v6.19: canvas 10 x 8 -> 7.20 x 5.50 in (within Nature double-col spec).
# Direct ggsave (not save_panel) so we can also emit an SVG copy
# alongside the PDF/PNG for Illustrator round-tripping.
ggsave(file.path(efig01_dir, "missingness_heatmap.pdf"), efig1,
       width = 7.20, height = 5.50, device = cairo_pdf, family = NM$font_family)
ggsave(file.path(efig01_dir, "missingness_heatmap.png"), efig1,
       width = 7.20, height = 5.50, dpi = 300, type = "cairo")
ggsave(file.path(efig01_dir, "missingness_heatmap.svg"), efig1,
       width = 7.20, height = 5.50, device = svglite::svglite,
       fix_text_size = FALSE)
saveRDS(efig1, file.path(efig01_dir, "missingness_heatmap.rds"))

cat("\n── Extended Data Fig. 1 saved to ExtendedData/ed_completeness/ ──\n")

# ── 7b. efig01 manifest fragment (missingness stats) ───────────────────────
# Keys are scoped so prose in Section 1 (cohort assembly) can cite missingness
# with inline  `r stats$edfig01$<key>`  expressions.

# Missingness per variable (averaged across studies, weighted by study n)
miss_by_var <- missingness_matrix %>%
  dplyr::group_by(variable, var_label) %>%
  dplyr::summarise(
    pct_available_mean  = mean(pct_available, na.rm = TRUE),
    pct_missing_mean    = 100 - mean(pct_available, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(pct_missing_mean)

# Missingness per study (averaged across variables)
miss_by_study <- missingness_matrix %>%
  dplyr::group_by(study_clean, n) %>%
  dplyr::summarise(
    pct_available_mean  = mean(pct_available, na.rm = TRUE),
    pct_missing_mean    = 100 - mean(pct_available, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(pct_missing_mean)

# Complete-case counts across the key_vars set (on genotyped subset)
complete_rows <- genotyped %>%
  dplyr::select(dplyr::all_of(setdiff(key_vars, "vaf_prop"))) %>%
  stats::complete.cases() %>%
  sum()

pct_missing_overall <- mean(
  100 - missingness_matrix$pct_available,
  na.rm = TRUE
)

efig01_fragment <- list(
  # scope
  n_vars_harmonised             = length(key_vars),
  n_genotyped_in_heatmap        = nrow(genotyped),
  n_studies_in_heatmap          = dplyr::n_distinct(missingness_matrix$study_clean),
  # summary
  pct_missing_overall           = pct_missing_overall,
  pct_available_overall         = 100 - pct_missing_overall,
  n_complete_case_key_vars      = complete_rows,
  pct_complete_case_key_vars    = 100 * complete_rows / nrow(genotyped),
  # full tables for Supplementary Table export + prose lookup
  missingness_by_variable       = miss_by_var,
  missingness_by_study          = miss_by_study,
  missingness_matrix_long       = missingness_matrix,
  # quick lookups for prose (best/worst cells)
  best_variable                 = as.character(miss_by_var$var_label[1]),
  best_variable_pct_available   = miss_by_var$pct_available_mean[1],
  worst_variable                = as.character(miss_by_var$var_label[nrow(miss_by_var)]),
  worst_variable_pct_available  = miss_by_var$pct_available_mean[nrow(miss_by_var)],
  best_study                    = as.character(miss_by_study$study_clean[1]),
  best_study_pct_available      = miss_by_study$pct_available_mean[1],
  worst_study                   = as.character(miss_by_study$study_clean[nrow(miss_by_study)]),
  worst_study_pct_available     = miss_by_study$pct_available_mean[nrow(miss_by_study)]
)

# ── 8. Export analysis-ready dataset ─────────────────────────────────────────

saveRDS(master, here("data", "processed", "bAVM_analysis_ready.rds"))
cat(sprintf("\n── Saved analysis-ready dataset: %d rows × %d columns ──\n",
  nrow(master), ncol(master)))
cat("New columns added: geno_binary, geno_variant, geno_status, recurrence_num, age_group\n")

# ── 9. Export genotyped-only subset for convenience ──────────────────────────

saveRDS(genotyped, here("data", "processed", "bAVM_genotyped_only.rds"))
cat(sprintf("Saved genotyped-only subset: %d rows\n", nrow(genotyped)))

# ── 10. Emit manifest fragment AFTER the upstream .rds is written.  ─────────
# check_stats_manifest.R compares the fragment mtime against
# bAVM_analysis_ready.rds; writing the fragment before the .rds means the
# fragment is always 0–10ms stale under strict mode. Keep this block LAST so
# the fragment mtime is guaranteed to be strictly later than the upstream
# data.
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(section = "edfig01", stats = efig01_fragment)

cat("\n══ 02_prep_analysis_dataset.R complete ══\n")
