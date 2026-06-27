#!/usr/bin/env Rscript
# 01_clean_master.R — Column-by-column cleaning of bAVM genotype-phenotype master dataset
#
# Input:  data/raw/phi_safe/bAVMgeno-pheno_master18May26.xlsx (sheet: "Data Collection")
#         (2026-05-18 bump: replaced bAVMgeno-pheno_master_04_28_master.xlsx;
#          older masters left in phi_unsafe/+phi_safe/ for provenance.)
#         Masked by analysis/00_data_prep/00_mask_phi.R from data/raw/phi_unsafe/
#         (the latter is gitignored; the masked copy is committed).
# Output: data/processed/bAVM_genopheno_clean.rds
#         data/processed/bAVM_genopheno_clean.csv
#         data/processed/cleaning_report.html (via column-level summaries)
#
# Usage:
#   Rscript analysis/00_data_prep/01_clean_master.R
#   # or source() interactively

# ── Load helpers ─────────────────────────────────────────────────────────────
source(here::here("analysis", "helper_scripts", "utils.R"))
load_common_packages()
library(readxl)

cat("═══ Column-by-Column Data Cleaning ═══\n\n")

# ── Load masked review IDs (audit flags) ─────────────────────────────────────
# MRNs in this file are PHI-masked (see analysis/00_data_prep/00_mask_phi.R)
# so they match the masked MRN column in the phi_safe spreadsheet.
.review_ids_path <- here::here("data", "flags", "review_ids.csv")
if (file.exists(.review_ids_path)) {
  .review_ids <- readr::read_csv(.review_ids_path, show_col_types = FALSE)
  .conflict_mrns       <- .review_ids$value[.review_ids$flag_type == "conflict_mrn"]
  .neg_vaf_patient_ids <- .review_ids$value[.review_ids$flag_type == "neg_vaf_patient_id"]
} else {
  warning(
    "data/flags/review_ids.csv not found; MRN-/patient_id-specific audit flags ",
    "will be empty. Create this file locally (gitignored) to enable those flags."
  )
  .conflict_mrns       <- character(0)
  .neg_vaf_patient_ids <- character(0)
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 0: Load raw data
# ══════════════════════════════════════════════════════════════════════════════
raw <- read_excel(
  here::here("data", "raw", "phi_safe",
             "bAVMgeno-pheno_master18May26.xlsx"),
  sheet = "Data Collection"
)
cat(sprintf("Loaded: %d rows × %d columns\n\n", nrow(raw), ncol(raw)))

# ── Assign clean short column names ──────────────────────────────────────────
# Store original names for reference
original_names <- names(raw)

names(raw) <- c(
  "mrn",                  # [1]  MRN
  "study_code",           # [2]  Study (coded)
  "patient_id",           # [3]  Patient ID (Disease)
  "sample_type",          # [4]  Sample type
  "mutation_code",        # [5]  Mutation
  "vaf",                  # [6]  VAF
  "sex",                  # [7]  Sex
  "age_dx",               # [8]  Age at diagnosis
  "age_surgery",          # [9]  Age at surgery
  "race",                 # [10] Race
  "sm_size",              # [11] Spetzler Martin size
  "sm_eloquence",         # [12] Spetzler Martin eloquence
  "sm_drainage",          # [13] Spetzler Martin deep venous drainage
  "sm_total",             # [14] Spetzler Martin total score
  "compact_nidus",        # [15] Compact nidus
  "ruptured_at_surgery",  # [16] Ruptured at time of surgery
  "ever_ruptured",        # [17] Has the bAVM ever ruptured
  "growing",              # [18] Was the bAVM growing over time
  "location_codes",       # [19] Location (comma-separated codes)
  "laterality",           # [20] Laterality
  "prior_seizure",        # [21] History of prior seizure
  "prior_radiation",      # [22] Prior radiation to bAVM
  "radiation_shrink",     # [23] Did radiation shrink by > 50%
  "prior_embolization",   # [24] Embolization prior to surgery
  "intranidal_aneurysm",  # [25] Intranidal aneurysm
  "venous_varix",         # [26] Venous varix
  "venous_outflow_stenosis", # [27] Venous outflow stenosis
  "flow_related_aneurysm",   # [28] Flow-related aneurysm
  "n_high_risk_features",    # [29] # of high risk features
  "notes"                    # [30] Notes
)

# ── Re-inclusion: dPCR-confirmed CHOP7 / CHOP8 (KRAS G12D) ────────────────────
# CHOP7 ([CHOP7 specimen accession]) and CHOP8 ([CHOP8 specimen accession]) were removed by the CHOP team while
# genotype was pending (04_17 → 04_21) and are absent from the 18 May master.
# Run 20 dPCR (2026-05-18) confirmed KRAS G12D in both, and the new UAB vascular
# IRB (A. Hale email, 2026-05-27) covers them regardless of surgical-consent
# status, so they are re-included here. Curated values come from the flag CSV;
# its "code (label)" cells are reduced to the bare numeric codes used in the raw
# master. This reverses the earlier "do-not-re-add" note (DATA_DECISIONS.md §0).
.chop78_path <- here::here("data", "flags", "dropped_chop7_chop8_dpcr_confirmed.csv")
if (file.exists(.chop78_path)) {
  .c78  <- readr::read_csv(.chop78_path, show_col_types = FALSE)
  .code <- function(x) sub("\\s.*$", "", as.character(x))   # "2 (Female)" -> "2"

  # VAF must be computed the SAME way as the rest of the cohort and the Methods:
  # TERT-referenced, VAF = mutant / (mutant + TERT) copies. The flag CSV's
  # `dpcr_run20_vaf_pct` is mutant/valid-partitions (a different denominator), so
  # we recompute from the raw Run 20 QIAcuity export (G + TERT channels) for each
  # case's well at the instrument-applied thresholds.
  .run20_dir <- here::here("data", "raw",
    "bAVM FFPE Samples Run 20 051826_RFU_img1_G_Y_O_R_C_Fr_21_05_2026_09_26_16_UTC-04_00")
  .read_run20 <- function(ch) {
    f <- list.files(.run20_dir, pattern = sprintf("img1_%s_[0-9].*\\.csv$", ch),
                    full.names = TRUE)[1]
    d <- data.table::fread(f, skip = 1, sep = ",", header = TRUE, showProgress = FALSE)
    data.table::setnames(d, make.names(names(d)))
    d
  }
  .g20 <- .read_run20("G"); .fr20 <- .read_run20("Fr")
  .aligned_vaf <- vapply(as.character(.c78$dpcr_run20_well), function(w) {
    gp <- .g20[Well == w, sum(Is.positive == 1, na.rm = TRUE)]
    tp <- .fr20[Well == w, sum(Is.positive == 1, na.rm = TRUE)]
    round(100 * gp / (gp + tp), 4)
  }, numeric(1))
  cat(sprintf("  [re-inclusion] TERT-referenced VAF (mutant/(mut+TERT)): %s\n",
              paste(sprintf("%s=%.2f%%", .c78$patient_id, .aligned_vaf), collapse = "; ")))
  .reinc <- tibble::tibble(
    mrn                 = as.character(.c78$mrn),
    study_code          = "3",                                # CHOP
    patient_id          = as.character(.c78$patient_id),
    sample_type         = .code(.c78$sample_type),            # 1 = FFPE
    mutation_code       = "1",                                # KRAS G12D (dPCR Run 20)
    vaf                 = as.character(.aligned_vaf),         # TERT-referenced, cohort-aligned
    sex                 = .code(.c78$sex),
    age_dx              = as.character(.c78$age_at_diagnosis),
    age_surgery         = as.character(.c78$age_at_surgery),
    sm_size             = as.character(.c78$sm_size),
    sm_eloquence        = .code(.c78$sm_eloquence),
    sm_drainage         = .code(.c78$sm_deep_venous),
    sm_total            = as.character(.c78$sm_total),
    compact_nidus       = .code(.c78$compact_nidus),
    ruptured_at_surgery = .code(.c78$ruptured_at_surgery),
    ever_ruptured       = .code(.c78$ever_ruptured),
    location_codes      = .code(.c78$location),
    laterality          = .code(.c78$laterality),
    venous_varix        = .code(.c78$venous_varix),
    notes               = "Re-included 2026-05-27: dPCR-confirmed KRAS G12D (Run 20); new vascular IRB (Hale email)."
  )
  # Match each column's type to the raw frame (read_excel infers numeric for
  # some coded columns) so bind_rows does not error on type clashes.
  for (.col in intersect(names(.reinc), names(raw))) {
    .reinc[[.col]] <- if (is.numeric(raw[[.col]]))
        suppressWarnings(as.numeric(.reinc[[.col]])) else as.character(.reinc[[.col]])
  }
  raw <- dplyr::bind_rows(raw, .reinc)
  cat(sprintf("  [re-inclusion] +%d dPCR-confirmed CHOP cases (CHOP7/CHOP8); now %d rows\n",
              nrow(.reinc), nrow(raw)))
}

# ── Coerce mutation_code to character ─────────────────────────────────────────
# In the 04_17 master this column contained mixed strings ("Tissue available")
# and numeric codes, so readxl returned character. In the 04_21 master the
# free-text entries were removed, so readxl now infers numeric (double).
# Downstream logic compares against string codes ("1", "2", ...), so coerce
# here to keep a single code path.
if (!is.character(raw[["mutation_code"]])) {
  raw[["mutation_code"]] <- as.character(raw[["mutation_code"]])
}

# ══════════════════════════════════════════════════════════════════════════════
# HELPER: Standardize missing values
# Many columns use "-", "N/A", or empty string alongside true NA.
# After this step, all missing values are represented as R's native NA.
# ══════════════════════════════════════════════════════════════════════════════
na_strings <- c("-", "N/A", "n/a", "")

to_na <- function(x) {
  x[x %in% na_strings] <- NA
  x
}

# Report what will be standardized
cat("Standardizing missing values → NA\n")
total_converted <- 0L
for (col in names(raw)) {
  if (is.character(raw[[col]])) {
    n_dash  <- sum(raw[[col]] == "-", na.rm = TRUE)
    n_na    <- sum(raw[[col]] == "N/A", na.rm = TRUE)
    n_other <- sum(raw[[col]] %in% c("n/a", ""), na.rm = TRUE)
    n_total <- n_dash + n_na + n_other
    if (n_total > 0) {
      detail <- paste(c(
        if (n_dash > 0) sprintf('"-"=%d', n_dash),
        if (n_na > 0)   sprintf('"N/A"=%d', n_na),
        if (n_other > 0) sprintf('other=%d', n_other)
      ), collapse = ", ")
      cat(sprintf("  %-30s %d values → NA (%s)\n", col, n_total, detail))
      total_converted <- total_converted + n_total
    }
  }
}
cat(sprintf("  Total: %d non-standard missing values → NA across %d columns\n\n",
            total_converted,
            sum(sapply(names(raw), function(col)
              is.character(raw[[col]]) && any(raw[[col]] %in% na_strings, na.rm = TRUE)))))

# Apply to all character columns
raw <- raw %>%
  mutate(across(where(is.character), to_na))

# ── Pre-exclusion per-study counts (for CONSORT Tier 1) ──────────────────────
# Fig 1A's source-cohort boxes need pre-exclusion counts (i.e., before the
# spinal-AVM rows below are dropped) so that Tier 1 totals add up to the Tier 2
# pre-exclusion total. Spine cases all carry an NP- prefix and would be imputed
# to UAB by the §[2] block below; we apply the same prefix-imputation logic
# here so spine cases are attributed to UAB pre-exclusion. Rows whose
# Study (coded) is NA AND whose patient_id has no recognized prefix
# (currently the single "or 1923741/1923740" orphan) are bucketed as
# "Unassigned" so they remain visible in the figure as a pending-resolution
# case rather than silently dropping out of Tier 1.
.pre_excl_prefix_to_study <- c(
  "BCH" = 1, "AVMUAB" = 2, "NP" = 2, "B" = 2,
  "CHOP" = 3, "Nikolaev" = 4, "Priemer" = 5,
  "Hong" = 6, "Goss" = 7, "Gao" = 8
)
.pre_excl_study_labels <- c(
  "1" = "BCH", "2" = "UAB", "3" = "CHOP", "4" = "Nikolaev",
  "5" = "Priemer", "6" = "Hong", "7" = "Goss", "8" = "Gao"
)
.pre_excl_pid_prefix <- stringr::str_extract(raw[["patient_id"]], "^[A-Za-z]+")
.pre_excl_study_code <- as.character(raw[["study_code"]])
.pre_excl_study_code <- ifelse(
  is.na(.pre_excl_study_code) & .pre_excl_pid_prefix %in% names(.pre_excl_prefix_to_study),
  as.character(.pre_excl_prefix_to_study[.pre_excl_pid_prefix]),
  .pre_excl_study_code
)
.pre_excl_study <- .pre_excl_study_labels[.pre_excl_study_code]
.pre_excl_study[is.na(.pre_excl_study)] <- "Unassigned"
n_per_study_pre_exclusion <- as.list(table(.pre_excl_study))
rm(.pre_excl_prefix_to_study, .pre_excl_study_labels,
   .pre_excl_pid_prefix, .pre_excl_study_code, .pre_excl_study)

# ── Exclude spinal AVMs (brain AVM study only) ────────────────────────────────
# The 04_17 master relabeled 4 UAB pathology cases as spinal AVMs with "(spine)"
# suffix in patient_id. These are excluded from the brain AVM cohort.
n_pre_exclusion <- nrow(raw)
.exclusion_log <- list()  # accumulator for cohort_exclusions.rds

spine_idx <- grepl("(spine)", raw[["patient_id"]], fixed = TRUE)
n_spine <- sum(spine_idx)
if (n_spine > 0) {
  cat(sprintf("Excluding %d spinal AVM(s) from brain AVM cohort:\n", n_spine))
  for (pid in raw[["patient_id"]][spine_idx]) cat(sprintf("  %s\n", pid))
  raw <- raw[!spine_idx, ]
  cat(sprintf("Remaining: %d rows\n\n", nrow(raw)))
}
.exclusion_log[["Spinal AVM"]] <- as.integer(n_spine)

# Persist the exclusion bookkeeping so 03_prep_cohort_counts.R (and the Fig 1A
# CONSORT producer) can render pre/post-exclusion totals from live data
# instead of hardcoded numbers.
exclusion_log <- list(
  n_pre_exclusion           = as.integer(n_pre_exclusion),
  n_post_exclusion          = as.integer(nrow(raw)),
  n_excluded                = as.integer(n_pre_exclusion - nrow(raw)),
  reasons                   = .exclusion_log,
  n_per_study_pre_exclusion = n_per_study_pre_exclusion
)
.exclusion_out_dir <- here::here("data", "processed")
dir.create(.exclusion_out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(exclusion_log,
        file.path(.exclusion_out_dir, "cohort_exclusions.rds"))
rm(n_pre_exclusion, .exclusion_log, exclusion_log, .exclusion_out_dir,
   n_per_study_pre_exclusion)

# ══════════════════════════════════════════════════════════════════════════════
# COLUMN-BY-COLUMN CLEANING
# ══════════════════════════════════════════════════════════════════════════════

df <- raw  # working copy
df$row_num <- seq_len(nrow(df))

# ══════════════════════════════════════════════════════════════════════════════
# PUBLISHED DATA INTEGRATION — fill gaps from manuscript tables
# Sources: Gao Suppl Table 1, Hong Table 1, Nikolaev Tables S6/S8
# Track all fills & inconsistencies in .pub_note (appended to review_note later)
# ══════════════════════════════════════════════════════════════════════════════
cat("═══ Published Data Integration ═══\n")
.pub_note <- rep(NA_character_, nrow(df))

# ── Fix Excel date-serial corruption in location_codes ──────────────────────
# Comma-separated codes "7,8" etc. were auto-formatted as dates by Excel
excel_date_fixes <- c(
  "45846" = "7,8",    # 7/8  → parietal + occipital
  "45664" = "1,7",    # 1/7  → temporal + parietal
  "45670" = "1,13",   # 1/13 → temporal + sylvian_fissure
  "45851" = "7,13"    # 7/13 → parietal + sylvian_fissure
)
n_date_fix <- 0
for (bad_val in names(excel_date_fixes)) {
  hits <- which(!is.na(df$location_codes) & df$location_codes == bad_val)
  if (length(hits) > 0) {
    cat(sprintf("  location_codes: %s → '%s' (%d rows)\n",
                bad_val, excel_date_fixes[bad_val], length(hits)))
    note_txt <- sprintf("location_codes: Excel date '%s' corrected to '%s'",
                        bad_val, excel_date_fixes[bad_val])
    for (h in hits) {
      .pub_note[h] <- if (is.na(.pub_note[h])) note_txt else paste(.pub_note[h], note_txt, sep = "; ")
    }
    df$location_codes[hits] <- excel_date_fixes[bad_val]
    n_date_fix <- n_date_fix + length(hits)
  }
}
# Fix literal "zero" → "0"
zero_hits <- which(!is.na(df$location_codes) & grepl("zero", df$location_codes))
if (length(zero_hits) > 0) {
  df$location_codes[zero_hits] <- gsub("zero", "0", df$location_codes[zero_hits])
  cat(sprintf("  location_codes: 'zero' → '0' (%d rows)\n", length(zero_hits)))
  note_txt <- "location_codes: 'zero' corrected to '0'"
  for (h in zero_hits) {
    .pub_note[h] <- if (is.na(.pub_note[h])) note_txt else paste(.pub_note[h], note_txt, sep = "; ")
  }
  n_date_fix <- n_date_fix + length(zero_hits)
}
cat(sprintf("  Total location_codes fixes: %d\n", n_date_fix))

# ── Fill Gao intranidal_aneurysm from Suppl Table 1 (Associated aneurysm) ──
gao_pub_path <- here::here("data", "processed", "published", "gao_published.csv")
if (file.exists(gao_pub_path)) {
  gao_pub <- read_csv(gao_pub_path, show_col_types = FALSE)
  gao_fill <- gao_pub %>%
    filter(!is.na(intranidal_aneurysm_pub)) %>%
    select(patient_id, intranidal_aneurysm_pub, aneurysm_count_pub)

  n_before <- sum(is.na(df$intranidal_aneurysm[df$study_code == 8 & !is.na(df$study_code)]))
  for (i in seq_len(nrow(gao_fill))) {
    idx <- which(df$patient_id == gao_fill$patient_id[i] & is.na(df$intranidal_aneurysm))
    if (length(idx) > 0) {
      val <- gao_fill$intranidal_aneurysm_pub[i]
      df$intranidal_aneurysm[idx] <- as.character(val)
      note_txt <- sprintf("intranidal_aneurysm filled=%d (Gao Suppl Table 1)", val)
      .pub_note[idx] <- ifelse(is.na(.pub_note[idx]), note_txt,
                               paste(.pub_note[idx], note_txt, sep = "; "))
    }
  }
  n_after <- sum(is.na(df$intranidal_aneurysm[df$study_code == 8 & !is.na(df$study_code)]))
  cat(sprintf("  Gao intranidal_aneurysm: filled %d/%d (from Suppl Table 1 'Associated aneurysm')\n",
              n_before - n_after, n_before))
} else {
  cat("  ⚠ gao_published.csv not found — skipping Gao fills\n")
}

# ── Per-row manual corrections (single registry) ─────────────────────────────
# Every hand-curated per-patient value override lives in ONE place:
#   data/flags/manual_corrections.csv  (patient_id, field, from, to,
#                                        pub_data_note, decision_ref)
# applied here in one loop instead of scattered one-off blocks. Each row sets
# `field` = `to` for the named `patient_id`, but only where the current value
# matches `from` (use the literal NA for "fill only when missing"), and records
# `pub_data_note` for provenance. To add a future per-row fix, append a CSV row
# — do NOT add a new code block. (Rule-based corrections such as the SM-total
# component recompute (§40) and the rupture reconciliation [17b] stay as code;
# published-source fills come from data/raw/<study>/. See data/flags/README.md.)
.mc_path <- here::here("data", "flags", "manual_corrections.csv")
if (file.exists(.mc_path)) {
  .mc <- readr::read_csv(.mc_path, show_col_types = FALSE,
                         col_types = readr::cols(.default = "c"))
  for (.i in seq_len(nrow(.mc))) {
    .pid <- .mc$patient_id[.i]; .fld <- .mc$field[.i]
    .frm <- .mc$from[.i];       .to  <- .mc$to[.i]
    if (!.fld %in% names(df)) {
      warning("manual_corrections.csv: unknown field '", .fld, "' for ", .pid, " — skipped")
      next
    }
    .match_from <- if (is.na(.frm) || .frm == "NA") is.na(df[[.fld]]) else df[[.fld]] == .frm
    .idx <- which(df$patient_id == .pid & .match_from)
    if (length(.idx) > 0) {
      df[[.fld]][.idx] <- .to
      .note_txt <- .mc$pub_data_note[.i]
      .pub_note[.idx] <- ifelse(is.na(.pub_note[.idx]), .note_txt,
                                paste(.pub_note[.idx], .note_txt, sep = "; "))
      cat(sprintf("  [manual_corrections] %s %s: %s → %s\n",
                  .pid, .fld, ifelse(is.na(.frm) || .frm == "NA", "NA", .frm), .to))
    }
  }
} else {
  warning("data/flags/manual_corrections.csv not found; per-row manual corrections skipped")
}

# ── Cohort SM-total recomputation ────────────────────────────────────────────
# Policy (DATA_DECISIONS.md §40, 2026-05-19): when sm_total ≠ sm_size +
# sm_eloquence + sm_drainage in a cohort sample (sample_type != 3 /
# Literature), overwrite the reported total with the component sum.
# Components are the radiologist/surgeon-graded inputs; the total is
# derivative arithmetic, so the component sum is the higher-evidence value.
# Distinct from §17/§21 where we trust the *published* Hong total over its
# components — that policy is preserved by the sample_type filter below.
# Provenance is recorded per row in `pub_data_note` ("sm_total corrected
# X→Y (components: a+b+c = Y; trust components per DATA_DECISIONS §40)").
#
# NB: this block runs BEFORE the [11]-[14] column-rename pass below, so the
# columns are still `sm_eloquence` and `sm_drainage` (not the later
# `sm_eloq_num` / `sm_drain_num`). If this block is ever moved past the
# rename pass, the column references will need to be updated accordingly.
.sm_size_int  <- suppressWarnings(as.integer(df$sm_size))
.sm_eloq_int  <- suppressWarnings(as.integer(df$sm_eloquence))
.sm_drain_int <- suppressWarnings(as.integer(df$sm_drainage))
.sm_total_int <- suppressWarnings(as.integer(df$sm_total))
.sm_computed  <- .sm_size_int + .sm_eloq_int + .sm_drain_int
.sample_type_int <- suppressWarnings(as.integer(df$sample_type))
.is_cohort <- !is.na(.sample_type_int) & .sample_type_int != 3   # 3 = Literature
.sm_mismatch_idx <- which(.is_cohort &
                          !is.na(.sm_computed) & !is.na(.sm_total_int) &
                          .sm_computed != .sm_total_int)
for (i in .sm_mismatch_idx) {
  note_txt <- sprintf(
    "sm_total corrected %s→%d (components: %s+%s+%s = %d; trust components per DATA_DECISIONS §40)",
    df$sm_total[i], .sm_computed[i],
    df$sm_size[i], df$sm_eloquence[i], df$sm_drainage[i], .sm_computed[i])
  df$sm_total[i] <- as.character(.sm_computed[i])
  .pub_note[i] <- ifelse(is.na(.pub_note[i]), note_txt,
                          paste(.pub_note[i], note_txt, sep = "; "))
  cat(sprintf("  Cohort SM total corrected: %s — total %d→%d (components sum)\n",
              df$patient_id[i], .sm_total_int[i], .sm_computed[i]))
}
if (length(.sm_mismatch_idx) > 0) {
  cat(sprintf("  → %d cohort row(s) had sm_total recomputed from components\n",
              length(.sm_mismatch_idx)))
}

# ── Flag Nikolaev_AVM_R6: ever_ruptured inconsistency ───────────────────────
# Master=0 but Nikolaev Table S8 says rupture=yes
r6_idx <- which(df$patient_id == "Nikolaev_AVM_R6")
if (length(r6_idx) > 0) {
  note_txt <- "ever_ruptured inconsistency: master=0 but Nikolaev Table S8 says rupture=yes"
  .pub_note[r6_idx] <- ifelse(is.na(.pub_note[r6_idx]), note_txt,
                              paste(.pub_note[r6_idx], note_txt, sep = "; "))
  cat("  ⚠ Nikolaev_AVM_R6: ever_ruptured master=0 vs Table S8=yes (flagged)\n")
}

# ── Fill Nikolaev prior_seizure from Table S8 (Finnish, Epilepsy at presentation)
nik_pub_path <- here::here("data", "processed", "published", "nikolaev_published.csv")
if (file.exists(nik_pub_path)) {
  nik_pub <- read_csv(nik_pub_path, show_col_types = FALSE)
  nik_fill <- nik_pub %>%
    filter(!is.na(prior_seizure_pub)) %>%
    select(patient_id, prior_seizure_pub)

  n_before <- sum(is.na(df$prior_seizure[grepl("Nikolaev_AVM_R", df$patient_id)]))
  for (i in seq_len(nrow(nik_fill))) {
    idx <- which(df$patient_id == nik_fill$patient_id[i] & is.na(df$prior_seizure))
    if (length(idx) > 0) {
      val <- nik_fill$prior_seizure_pub[i]
      df$prior_seizure[idx] <- as.character(val)
      note_txt <- sprintf("prior_seizure filled=%d (Nikolaev Table S8)", val)
      .pub_note[idx] <- ifelse(is.na(.pub_note[idx]), note_txt,
                               paste(.pub_note[idx], note_txt, sep = "; "))
    }
  }
  n_after <- sum(is.na(df$prior_seizure[grepl("Nikolaev_AVM_R", df$patient_id)]))
  cat(sprintf("  Nikolaev Finnish prior_seizure: filled %d/%d (from Table S8 'Epilepsy at presentation')\n",
              n_before - n_after, n_before))
} else {
  cat("  ⚠ nikolaev_published.csv not found — skipping Nikolaev fills\n")
}

# ── Fill Priemer prior_seizure from published 'Presentation' column ────────
# Priemer et al. reports presenting symptom; "Seizure" → prior_seizure=1
# Cannot infer prior_seizure=0 from other presentations (hemorrhage, headache
# as presenting symptom does not rule out seizure history)
pri_pub_path <- here::here("data", "processed", "published", "priemer_published.csv")
if (file.exists(pri_pub_path)) {
  pri_pub <- read_csv(pri_pub_path, show_col_types = FALSE)
  pri_seizure <- pri_pub %>%
    filter(presentation == "Seizure") %>%
    pull(patient_id)

  n_before <- sum(is.na(df$prior_seizure[df$patient_id %in% pri_pub$patient_id]))
  n_filled <- 0
  for (pid in pri_seizure) {
    idx <- which(df$patient_id == pid & is.na(df$prior_seizure))
    if (length(idx) > 0) {
      df$prior_seizure[idx] <- "1"
      note_txt <- "prior_seizure filled=1 (Priemer: presenting symptom was seizure)"
      .pub_note[idx] <- ifelse(is.na(.pub_note[idx]), note_txt,
                               paste(.pub_note[idx], note_txt, sep = "; "))
      n_filled <- n_filled + length(idx)
    }
  }
  cat(sprintf("  Priemer prior_seizure: filled %d (seizure presentation → 1; %d remain NA — cannot infer 0)\n",
              n_filled, n_before - n_filled))
} else {
  cat("  ⚠ priemer_published.csv not found — skipping Priemer fills\n")
}

n_pub_flagged <- sum(!is.na(.pub_note))
cat(sprintf("  Published data notes: %d rows annotated\n", n_pub_flagged))

cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# IDENTIFIER CLEANING (columns 1–3): mrn, study_code, patient_id
# Original columns are NEVER modified. All fixes go into *_clean columns.
# ══════════════════════════════════════════════════════════════════════════════

# ── [1] mrn — kept as-is (character identifier, PHI) ─────────────────────────
cat("[1/30] mrn: character identifier — no transformation\n")
cat(sprintf("  %d non-NA, %d NA\n", sum(!is.na(df$mrn)), sum(is.na(df$mrn))))

# ── [2] study_code → study_code_clean (impute 6 NA rows from patient_id prefix)
cat("[2/30] study_code → study_code_clean, study, study_clean\n")
study_labels <- c(
  "1" = "BCH", "2" = "UAB", "3" = "CHOP", "4" = "Nikolaev",
  "5" = "Priemer", "6" = "Hong", "7" = "Goss", "8" = "Gao"
)

# Prefix → study_code mapping (derived from complete rows)
# VCHOPUC added 2026-05-19 for the 3 CHOP orphan rows in the 18 May 26
# master (VCHOPUC066 + 2 numeric pathology-id rows). See DATA_DECISIONS §43.
prefix_to_study <- c(
  "BCH" = 1, "AVMUAB" = 2, "NP" = 2, "B" = 2,
  "CHOP" = 3, "VCHOPUC" = 3,
  "Nikolaev" = 4, "Priemer" = 5,
  "Hong" = 6, "Goss" = 7, "Gao" = 8
)

df <- df %>%
  mutate(
    pid_prefix = str_extract(patient_id, "^[A-Za-z]+"),
    # CHOP's institutional pathology numbering uses a 7-digit serial in the
    # 192xxxx range (visible as paired accessions in [specimen accession]-... rows: e.g.
    # `[specimen accession]-93 ([paired specimen accessions])`). The 18 May 26 master added two orphan
    # rows with only the bare 192xxxx number as patient_id (1927848 and
    # 1923749). They imputeto CHOP (study 3) by this pattern. Same rule
    # for `7316-` parenthetical (CHOP-CAG identifier prefix).
    study_code_clean = case_when(
      !is.na(study_code) ~ study_code,
      pid_prefix %in% names(prefix_to_study) ~ prefix_to_study[pid_prefix],
      str_detect(patient_id, "^192[0-9]{4}") ~ 3,           # CHOP pathology serial
      str_detect(patient_id, "7316-[0-9]+") ~ 3,            # CHOP-CAG identifier
      # CHOP46 phenotype-ambiguous orphan: patient_id "or 1923741/1923740"
      # is the secondary specimen of the CHOP46-1 AVM ambiguity documented
      # in DATA_DECISIONS §8. The phenotype data could not be confidently
      # linked to this specimen so demographics/location were cleared to
      # NA, but the specimen itself is a CHOP-genotyped FFPE sample and
      # belongs to the CHOP institutional cohort. The phenotype-NA status
      # remains encoded by phenotype_pending = TRUE elsewhere; no analysis
      # depends on the missing demographics.
      str_detect(patient_id, "1923741|1923740") ~ 3,
      TRUE ~ NA_real_
    ),
    study_code_imputed = is.na(study_code) & !is.na(study_code_clean),
    study = factor(study_labels[as.character(study_code)],
                   levels = unname(study_labels)),
    study_clean = factor(study_labels[as.character(study_code_clean)],
                         levels = unname(study_labels))
  )

n_imputed <- sum(df$study_code_imputed)
cat(sprintf("  Imputed %d NA study_codes from patient_id prefix:\n", n_imputed))
imputed_rows <- df %>% filter(study_code_imputed)
for (i in seq_len(nrow(imputed_rows))) {
  cat(sprintf("    Row %d: patient_id=%s → study=%s\n",
              imputed_rows$row_num[i], imputed_rows$patient_id[i],
              as.character(imputed_rows$study_clean[i])))
}
cat(sprintf("  %s\n", paste(capture.output(table(df$study_clean, useNA = "ifany")), collapse = "\n  ")))

# ── [3] patient_id → patient_id_clean (standardize formatting) ───────────────
cat("[3/30] patient_id → patient_id_clean\n")

df <- df %>%
  mutate(
    # Detect " (spatial)" suffix (case-insensitive) BEFORE strip — 7 BCH
    # rows in the 18 May 26 master gained this annotation to flag
    # spatial-transcriptomics samples. Boolean surfaced as
    # `spatial_sample` so downstream producers can subset on it. See
    # DATA_DECISIONS §44.
    spatial_sample = str_detect(patient_id,
                                 regex("\\s*\\(spatial\\)\\s*$", ignore_case = TRUE)),
    patient_id_clean = str_replace(patient_id,
                                    regex("\\s*\\(spatial\\)\\s*$", ignore_case = TRUE), ""),
    # Replace en-dash (–), em-dash (—) with underscore
    patient_id_clean = str_replace_all(patient_id_clean, "[\u2013\u2014]", "_"),
    # Replace semicolons with underscore
    patient_id_clean = str_replace_all(patient_id_clean, ";", "_"),
    # Replace commas with underscore
    patient_id_clean = str_replace_all(patient_id_clean, ",", "_"),
    # Collapse multiple spaces to single space, then trim
    patient_id_clean = str_squish(patient_id_clean),
    # Replace remaining spaces with underscore
    patient_id_clean = str_replace_all(patient_id_clean, " ", "_")
  )

.n_spatial <- sum(df$spatial_sample, na.rm = TRUE)
if (.n_spatial > 0) {
  cat(sprintf("  spatial_sample = TRUE: %d row(s) (suffix stripped from patient_id_clean)\n",
              .n_spatial))
}

# Report formatting changes
changed <- df %>% filter(patient_id != patient_id_clean)
cat(sprintf("  Formatting standardized in %d patient_ids:\n", nrow(changed)))
for (i in seq_len(nrow(changed))) {
  cat(sprintf("    Row %3d: \"%s\" → \"%s\"\n",
              changed$row_num[i], changed$patient_id[i], changed$patient_id_clean[i]))
}

n_dup <- sum(duplicated(df$patient_id))
cat(sprintf("  %d unique IDs, %d duplicates\n", n_distinct(df$patient_id), n_dup))

# ── Derive UDN (Unified Designation Number) and sample_uid ────────────────────
# UDN = patient-level ID with embedded verification info:
#   Format: {STUDY}{NNN}_{SEX}{AGE}
#     STUDY = abbreviated study name (BCH, UAB, CHOP, NIK, PRI, HON, GOS, GAO)
#     NNN   = zero-padded sequential number within study
#     SEX   = F (female), M (male), U (unknown/NA)
#     AGE   = integer age at diagnosis, zero-padded 2 digits; XX if NA
#
# sample_uid = row-level ID: {UDN}_S{N} where N = sample number within patient
#
# Patient grouping logic:
#   - Rows sharing the same MRN → same patient → same UDN
#   - Rows without MRN → each gets a unique UDN
#   - For multi-sample patients, sex/age_dx are taken from the first row
#     (confirmed identical across all duplicate-MRN pairs in the data)
cat("\n  Deriving UDN (unified designation number) and sample_uid...\n")

# Study abbreviation map
study_abbrev <- c(
  "1" = "BCH", "2" = "UAB", "3" = "CHOP", "4" = "NIK",
  "5" = "PRI", "6" = "HON", "7" = "GOS", "8" = "GAO"
)

# Pre-extract sex and age_dx as numeric for UDN encoding
# (these columns are still raw character/numeric at this point)
df <- df %>%
  mutate(
    .sex_raw = as.integer(sex),
    .age_dx_raw = round(as.numeric(age_dx), 0)
  )

# Step 1: Assign patient group IDs based on MRN linkage
df <- df %>%
  mutate(.patient_group = NA_integer_)

pg_counter <- 0L

# First pass: rows with MRN — group by MRN
mrn_groups <- df %>%
  filter(!is.na(mrn)) %>%
  distinct(mrn) %>%
  pull(mrn)

for (m in mrn_groups) {
  pg_counter <- pg_counter + 1L
  df$.patient_group[df$mrn == m & !is.na(df$mrn)] <- pg_counter
}

# Second pass: rows without MRN — each gets unique group
for (i in which(is.na(df$.patient_group))) {
  pg_counter <- pg_counter + 1L
  df$.patient_group[i] <- pg_counter
}

# Step 2: Build UDN string for each patient group
# Use first row's values for study, sex, age_dx within each group
patient_info <- df %>%
  group_by(.patient_group) %>%
  summarise(
    .study_code = first(na.omit(study_code_clean)),
    .sex_val    = first(na.omit(.sex_raw)),
    .age_val    = first(na.omit(.age_dx_raw)),
    .keep = "none"
  ) %>%
  ungroup()

# Assign within-study sequential numbers
patient_info <- patient_info %>%
  mutate(.study_abbrev = study_abbrev[as.character(.study_code)]) %>%
  group_by(.study_abbrev) %>%
  mutate(.seq = row_number()) %>%
  ungroup() %>%
  mutate(
    .sex_code = case_when(
      .sex_val == 1 ~ "F",
      .sex_val == 2 ~ "M",
      TRUE ~ "U"
    ),
    .age_code = ifelse(is.na(.age_val), "XX", sprintf("%02d", as.integer(.age_val))),
    udn = sprintf("%s%03d_%s%s", .study_abbrev, .seq, .sex_code, .age_code)
  )

# Join UDN back to main dataframe
df <- df %>%
  left_join(patient_info %>% select(.patient_group, udn), by = ".patient_group")

# Step 3: Assign sample number within each UDN and build sample_uid
df <- df %>%
  group_by(udn) %>%
  mutate(
    n_samples = n(),
    sample_num = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    sample_uid = sprintf("%s_S%d", udn, sample_num),
    multi_sample_flag = (n_samples > 1)
  )

# Step 3b: Multi-sample MRN classification + lesion_id (2026-05-19)
# ── Policy (DATA_DECISIONS.md §41) ─────────────────────────────────────────
# Every multi-sample MRN group is classified into one of five classes,
# loaded from a curator-maintained CSV at data/flags/multi_sample_classification.csv:
#
#   single_sample           : default (one row per MRN)
#   intra_avm_mosaic        : same surgery, multiple samples, mutation discord
#                             within one AVM (e.g. AVMUAB019 FFPE+G12D vs Tissue+neg).
#                             Keep all rows; rows share a single lesion_id.
#   duplicate_block         : same surgery, multiple FFPE blocks, identical
#                             clinical fields (e.g. CHOP62 VCHOPUC005_1 / _2).
#                             COLLAPSE to the first sample row; discarded sample's
#                             VAF is appended to pub_data_note for provenance.
#   recurrent_same_mut      : two surgeries at different ages with the same
#                             mutation; each AVM event is its own lesion.
#   recurrent_mut_discord   : two surgeries at different ages with DIFFERENT
#                             mutations (independent clonal events); each AVM
#                             event is its own lesion. CHOP91, CHOP23 in 18 May 26.
#   excluded_spinal         : both samples are spinal AVMs; rows already
#                             filtered at top of cleaner so the class is
#                             documentation-only.
#
# lesion_id encoding:
#   recurrent_*    → lesion_id = sample_uid   (each row = its own AVM event)
#   everything else → lesion_id = udn         (one lesion per patient)
#
# is_event_primary_row: TRUE for the canonical row per lesion. Used by
# 02_prep_analysis_dataset.R to subset to AVM-event-level (n=473).
.msc_path <- here::here("data", "flags", "multi_sample_classification.csv")
if (file.exists(.msc_path)) {
  .msc <- readr::read_csv(.msc_path, show_col_types = FALSE)
  # Join by masked MRN (cleaner runs against phi_safe/, so df$mrn already masked).
  df <- df %>%
    left_join(.msc %>% select(mrn_masked, class, collapse) %>%
                rename(mrn = mrn_masked,
                       multi_sample_class = class,
                       multi_sample_collapse = collapse),
              by = "mrn") %>%
    mutate(
      multi_sample_class = ifelse(is.na(multi_sample_class),
                                   "single_sample", multi_sample_class),
      multi_sample_collapse = ifelse(is.na(multi_sample_collapse),
                                      FALSE, multi_sample_collapse),
      lesion_id = ifelse(multi_sample_class %in%
                            c("recurrent_same_mut", "recurrent_mut_discord"),
                          sample_uid, udn)
    )
  cat("\n  Multi-sample classification (from data/flags/multi_sample_classification.csv):\n")
  cat(sprintf("    %s\n", paste(capture.output(table(df$multi_sample_class, useNA = "ifany")), collapse = "\n    ")))

  # Duplicate_block annotation: record the secondary block's VAF in the
  # surviving S1 row's pub_data_note (provenance for the collapse that
  # happens downstream in 02_prep_analysis_dataset.R). NO row drop here — the
  # cleaner keeps every sample-level row so .pub_note's length stays in
  # sync with df. Subsetting to one row per lesion event happens via
  # `is_event_primary_row` in the prepare-analysis step.
  .collapse_pairs <- df %>%
    filter(multi_sample_collapse, n_samples > 1) %>%
    select(udn, sample_num, vaf, sample_uid) %>%
    group_by(udn) %>%
    summarise(
      vaf_s1 = vaf[sample_num == 1][1],
      vaf_other = paste(vaf[sample_num > 1], collapse = ","),
      sample_uid_other = paste(sample_uid[sample_num > 1], collapse = ","),
      .groups = "drop"
    )
  if (nrow(.collapse_pairs) > 0) {
    for (i in seq_len(nrow(.collapse_pairs))) {
      udn_i <- .collapse_pairs$udn[i]
      note_txt <- sprintf(
        "duplicate_block collapse: S1 kept (vaf=%s); secondary %s (vaf=%s) suppressed from event-level analyses per DATA_DECISIONS §41",
        as.character(.collapse_pairs$vaf_s1[i]),
        .collapse_pairs$sample_uid_other[i],
        .collapse_pairs$vaf_other[i])
      idx <- which(df$udn == udn_i & df$sample_num == 1)
      .pub_note[idx] <- ifelse(is.na(.pub_note[idx]), note_txt,
                                paste(.pub_note[idx], note_txt, sep = "; "))
      cat(sprintf("    duplicate_block annotation: %s — S1 kept as primary event row, %s suppressed downstream\n",
                  udn_i, .collapse_pairs$sample_uid_other[i]))
    }
  }
  # is_event_primary_row: for duplicate_block, only S1 is primary (S2+
  # secondary blocks are non-primary and will be filtered out at the
  # event level). For every other class (single, mosaic, recurrent), the
  # first sample row per lesion_id is the canonical row. For recurrent
  # classes lesion_id = sample_uid so every row is its own primary.
  df <- df %>%
    group_by(lesion_id) %>%
    mutate(is_event_primary_row = (sample_num == min(sample_num))) %>%
    ungroup() %>%
    mutate(
      is_event_primary_row = ifelse(multi_sample_collapse & sample_num > 1,
                                     FALSE, is_event_primary_row)
    )
  cat(sprintf("    lesion_id: %d distinct AVM events across %d rows (%d primary, %d non-primary)\n",
              n_distinct(df$lesion_id), nrow(df),
              sum(df$is_event_primary_row), sum(!df$is_event_primary_row)))
} else {
  warning("data/flags/multi_sample_classification.csv not found — ",
          "multi_sample_class + lesion_id columns will not be created. ",
          "Create the file locally to enable the AVM-event-level subsetting.",
          call. = FALSE)
  df <- df %>% mutate(
    multi_sample_class = "single_sample",
    multi_sample_collapse = FALSE,
    lesion_id = udn,
    is_event_primary_row = (sample_num == min(sample_num))
  )
}

# Step 4: Flag cases requiring investigator review with explanatory notes
# ── Data-quality flags ──
# Case 1: specific MRN(s) in .conflict_mrns — conflicting mutation results
#         (FFPE positive vs Tissue negative). Actual IDs loaded from
#         data/flags/review_ids.csv (gitignored).
# Case 2: specific patient_id(s) in .neg_vaf_patient_ids — negative dPCR call
#         (code 12) but non-zero VAF. Actual IDs loaded from sidecar.
# Case 3: 6 rows with NA study_code — missing all demographics
# Case 4: 4 rows with negative age gap (age_surgery < age_dx)
# ── Clinical flags from notes (CHOP chart review) ──
# Case 5: SM grade UNGRADEABLE
# Case 6: Possible non-bAVM (cervical cord, cavernoma)
# Case 7: Syndromic / hereditary (HHT, SDHB)
# Case 8: Other notable clinical context from notes
.age_dx_tmp  <- as.numeric(df$age_dx)
.age_surg_tmp <- as.numeric(df$age_surgery)
.age_gap <- .age_surg_tmp - .age_dx_tmp

# Build note-based flags by scanning the notes column (case-insensitive)
.notes_lower <- tolower(replace_na(df$notes, ""))

df <- df %>%
  mutate(
    .flag_conflict  = (mrn %in% .conflict_mrns & !is.na(mrn)),
    # neg_vaf case: was dPCR neg (code 12) in v0 but corrected to G12D (code 1) in 04_17 master.
    # Flag retained for audit trail; the correction is now in the source data.
    # Specific patient_id loaded from data/flags/review_ids.csv (gitignored).
    .flag_neg_vaf   = (patient_id %in% .neg_vaf_patient_ids & mutation_code == "12"),
    # Fire only when study_code_clean (post-imputation) is still unresolved;
    # rows whose raw study_code is NA but whose patient_id prefix mapped to a
    # known study are correctly attributed and should not be flagged.
    .flag_no_study  = is.na(study_code_clean),
    .flag_neg_gap   = (!is.na(.age_gap) & .age_gap < -0.01),
    .flag_ungradeable = str_detect(.notes_lower, "ungradeable"),
    .flag_non_bavm    = str_detect(.notes_lower, "cervical cord|cavernoma|not sure if (true )?brain avm|not sure if avm"),
    .flag_syndromic   = str_detect(.notes_lower, "\\bhht\\b|\\bsdhb\\b|hereditary"),
    .flag_recurrence  = str_detect(.notes_lower, "recur"),
    .flag_prior_resect = str_detect(.notes_lower, "prior resection|resection .* (osh|china|elsewhere)"),
    .flag_no_resection = str_detect(.notes_lower, "not sure (if|which).*resect|no op note|never.*resect|don.t see a second operation"),
    .flag_hong25 = (patient_id == "Hong_AVM_25"),
    review_note = case_when(
      # Data-quality flags (highest priority)
      .flag_conflict & .flag_neg_gap ~
        "Conflicting mutation: FFPE KRAS+ vs Tissue dPCR-; Negative age gap (surgery < diagnosis)",
      .flag_conflict ~
        "Conflicting mutation: FFPE KRAS G12D positive vs Tissue dPCR negative (same patient)",
      .flag_neg_vaf ~
        "Negative dPCR (code 12) but VAF=6.83%; likely data entry error",
      .flag_no_study ~
        "Missing study_code; missing all demographics",
      .flag_neg_gap ~
        sprintf("Negative age gap: dx=%.2f surg=%.2f (gap=%.2f yr)",
                .age_dx_tmp[row_num], .age_surg_tmp[row_num], .age_gap[row_num]),
      # Manuscript-cited flags
      .flag_hong25 ~
        "Hong Table 1: sm_size/sm_eloquence NA — nidus not measured due to emergency surgery (Patient 25)",
      # Clinical flags from notes
      .flag_non_bavm & .flag_ungradeable ~
        sprintf("Possible non-bAVM + UNGRADEABLE SM: %s", notes),
      .flag_non_bavm ~
        sprintf("Possible non-bAVM: %s", notes),
      .flag_ungradeable ~
        sprintf("UNGRADEABLE SM grade: %s", notes),
      .flag_syndromic ~
        sprintf("Syndromic/hereditary: %s", notes),
      .flag_no_resection ~
        sprintf("Uncertain resection history: %s", notes),
      .flag_recurrence & .flag_prior_resect ~
        sprintf("Recurrence + prior outside resection: %s", notes),
      .flag_recurrence ~
        sprintf("Recurrence case: %s", notes),
      .flag_prior_resect ~
        sprintf("Prior outside resection: %s", notes),
      TRUE ~ NA_character_
    ),
    needs_review = !is.na(review_note)
  ) %>%
  select(-starts_with(".flag_"))

rm(.age_dx_tmp, .age_surg_tmp, .age_gap, .notes_lower)

# ── Append published data notes to review_note ──────────────────────────────
# .pub_note was built during Published Data Integration (fills, corrections,
# inconsistencies). Store as pub_data_note column and append to review_note.
df$pub_data_note <- .pub_note
has_pub <- !is.na(.pub_note)
has_rev <- !is.na(df$review_note)
# Both: append pub note after existing review note
df$review_note[has_pub & has_rev] <- paste(df$review_note[has_pub & has_rev],
                                           .pub_note[has_pub & has_rev], sep = " | PUB: ")
# Pub only: create new review note with PUB prefix
df$review_note[has_pub & !has_rev] <- paste0("PUB: ", .pub_note[has_pub & !has_rev])
# Update needs_review
df$needs_review <- !is.na(df$review_note)
rm(.pub_note, has_pub, has_rev)

# Report flagged rows by category
n_flagged <- sum(df$needs_review, na.rm = TRUE)
cat(sprintf("  Flagged needs_review: %d rows total\n", n_flagged))
flagged <- df %>% filter(needs_review)
# Categorize for summary
cat("\n  ── Flag summary by category ──\n")
cats <- flagged %>%
  mutate(category = case_when(
    str_detect(review_note, "Conflicting mutation") ~ "Conflicting mutation",
    str_detect(review_note, "Negative dPCR")        ~ "Negative mutation + VAF",
    str_detect(review_note, "Missing study_code")    ~ "Missing study_code",
    str_detect(review_note, "Negative age gap")      ~ "Negative age gap",
    str_detect(review_note, "non-bAVM")              ~ "Possible non-bAVM",
    str_detect(review_note, "UNGRADEABLE")           ~ "Ungradeable SM",
    str_detect(review_note, "Syndromic")             ~ "Syndromic/hereditary",
    str_detect(review_note, "Uncertain resection")   ~ "Uncertain resection",
    str_detect(review_note, "Recurrence")            ~ "Recurrence case",
    str_detect(review_note, "Prior outside")         ~ "Prior outside resection",
    str_detect(review_note, "^PUB:")                 ~ "Published data fill/flag only",
    TRUE ~ "Other"
  ))
for (cat_name in unique(cats$category)) {
  n <- sum(cats$category == cat_name)
  rows <- cats %>% filter(category == cat_name)
  cat(sprintf("    %s (%d):\n", cat_name, n))
  for (i in seq_len(nrow(rows))) {
    cat(sprintf("      Row %3d: %s\n", rows$row_num[i], rows$patient_id[i]))
  }
}

# Report
cat(sprintf("  Total patients (unique UDN): %d\n", n_distinct(df$udn)))
cat(sprintf("  Total rows: %d\n", nrow(df)))
multi_sample <- df %>% filter(n_samples > 1)
cat(sprintf("  Multi-sample patients: %d patients (%d rows)\n",
            n_distinct(multi_sample$udn), nrow(multi_sample)))

if (nrow(multi_sample) > 0) {
  cat("  Multi-sample detail:\n")
  for (u in unique(multi_sample$udn)) {
    rows <- multi_sample %>% filter(udn == u)
    review_tag <- ifelse(any(rows$needs_review), " ⚠ NEEDS REVIEW", "")
    cat(sprintf("    %s (MRN %s): samples %s%s\n",
                u, rows$mrn[1],
                paste(rows$sample_uid, collapse = ", "),
                review_tag))
    for (j in seq_len(nrow(rows))) {
      cat(sprintf("      S%d: patient_id=%s | sample_type=%s | mutation=%s\n",
                  rows$sample_num[j], rows$patient_id[j],
                  rows$sample_type[j],
                  ifelse(is.na(rows$mutation_code[j]), "NA", rows$mutation_code[j])))
    }
  }
}

# Clean up temporary columns
df <- df %>%
  select(-starts_with("."))

# ── [4] sample_type → sample_type_f (labeled factor) ─────────────────────────
cat("[4/30] sample_type → sample_type_f\n")
sample_labels <- c("1" = "FFPE", "2" = "Tissue", "3" = "Literature",
                   "4" = "Coil biopsy", "5" = "Liquid biopsy")
df <- df %>%
  mutate(
    sample_type_f = factor(sample_labels[as.character(sample_type)],
                           levels = unname(sample_labels))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$sample_type_f, useNA = "ifany")), collapse = "\n  ")))

# ── [5] mutation_code → mutation (labeled factor) + mutation_gene + mutation_variant
cat("[5/30] mutation_code → mutation, mutation_gene, mutation_variant\n")
mutation_labels <- c(
  "1"  = "KRAS G12D",
  "2"  = "KRAS G12V",
  "3"  = "KRAS G12C",
  "4"  = "KRAS Q61H",
  "5"  = "KRAS c.191_196dup",
  "6"  = "BRAF V600E",
  "7"  = "BRAF Q636X",
  "8"  = "KRAS G12D/V negative",
  "9"  = "No pathogenic mutation",
  "10" = "KRAS negative",
  "11" = "Multi-gene panel negative",
  "12" = "Multiplex dPCR negative",
  "13" = "KRAS G12A"
)

# Handle "Tissue available" as a special case (tissue exists but not yet tested)
df <- df %>%
  mutate(
    mutation_code_clean = case_when(
      mutation_code == "Tissue available" ~ NA_character_,
      TRUE ~ mutation_code
    ),
    tissue_available_not_tested = (mutation_code == "Tissue available"),
    mutation = factor(mutation_labels[mutation_code_clean],
                      levels = unname(mutation_labels)),
    # Derive gene and positive/negative status
    mutation_positive = case_when(
      mutation_code_clean %in% c("1","2","3","4","5","6","7","13") ~ TRUE,
      mutation_code_clean %in% c("8","9","10","11","12")      ~ FALSE,
      TRUE ~ NA
    ),
    mutation_gene = case_when(
      mutation_code_clean %in% c("1","2","3","4","5","13") ~ "KRAS",
      mutation_code_clean %in% c("6","7")             ~ "BRAF",
      mutation_positive == FALSE                      ~ "Negative",
      TRUE ~ NA_character_
    ),
    mutation_gene = factor(mutation_gene, levels = c("KRAS", "BRAF", "Negative"))
  )

cat(sprintf("  %s\n", paste(capture.output(table(df$mutation, useNA = "ifany")), collapse = "\n  ")))
cat(sprintf("  Positive: %d, Negative: %d, NA: %d\n",
            sum(df$mutation_positive == TRUE, na.rm = TRUE),
            sum(df$mutation_positive == FALSE, na.rm = TRUE),
            sum(is.na(df$mutation_positive))))
cat(sprintf("  Tissue available (not yet tested): %d\n",
            sum(df$tissue_available_not_tested, na.rm = TRUE)))

# Binary indicator columns for each specific variant
df <- df %>%
  mutate(
    mut_KRAS_G12D          = as.integer(mutation_code_clean == "1"),
    mut_KRAS_G12V          = as.integer(mutation_code_clean == "2"),
    mut_KRAS_G12C          = as.integer(mutation_code_clean == "3"),
    mut_KRAS_Q61H          = as.integer(mutation_code_clean == "4"),
    mut_KRAS_c191_196dup   = as.integer(mutation_code_clean == "5"),
    mut_BRAF_V600E         = as.integer(mutation_code_clean == "6"),
    mut_BRAF_Q636X         = as.integer(mutation_code_clean == "7"),
    mut_KRAS_G12A          = as.integer(mutation_code_clean == "13")
  )
cat("  Variant binary indicators: mut_KRAS_G12D, mut_KRAS_G12V, mut_KRAS_G12C,\n")
cat("    mut_KRAS_Q61H, mut_KRAS_c191_196dup, mut_KRAS_G12A, mut_BRAF_V600E, mut_BRAF_Q636X\n")

# ── [6] vaf → vaf_prop (normalized to proportion 0–1) ────────────────────────
# All centers (incl. CHOP from the 04_28 master onward) report VAF as PERCENTAGE
# (0–100); divide by 100 to land on the genomics-standard proportion (0–1).
#
# Until the 04_21 master, CHOP rows were entered as proportion (e.g. 0.01 = 1%)
# and the cleaning script branched on study_code_clean == 3 to skip the divide.
# In the 04_28 master those rows were re-entered in percent form, so the
# branch was dropped. The verification block below hard-fails the build if a
# CHOP row's raw VAF looks like a proportion (< 0.05) — which would indicate
# a regression to the old convention and silently inflate downstream numbers
# by 100×.
cat("[6/30] vaf → vaf_prop (normalized proportion)\n")

df <- df %>%
  mutate(
    vaf_raw  = round(as.numeric(vaf), 4),
    vaf_prop = round(if_else(is.na(vaf_raw), NA_real_, vaf_raw / 100), 6),
    vaf_pct  = round(vaf_prop * 100, 4)  # convenience: percentage form
  )

# Sanity guard: scoped to CHOP because that center is the only one whose
# entry convention has historically flipped (proportion in 04_21 → percent
# in 04_28). Any non-NA CHOP raw VAF below 0.1 indicates a regression to
# proportion-form entry and would silently inflate downstream stats by 100×.
# Other centers' legitimate low-VAF measurements (e.g. Gao 0.04%) are not
# scoped here; their entry conventions have been stable across releases.
.chop_vaf_regression <- df %>%
  filter(!is.na(vaf_raw), study_code_clean == 3, vaf_raw < 0.1) %>%
  select(row_num, study_clean, vaf_raw)
if (nrow(.chop_vaf_regression) > 0) {
  msg <- paste0(
    "Possible CHOP VAF scale-convention regression: ",
    nrow(.chop_vaf_regression),
    " CHOP row(s) have raw VAF < 0.1, which suggests proportion-form entry ",
    "(reverting to the pre-04_28 convention). Confirm the source ",
    "spreadsheet enters CHOP VAF as percent (0–100), or update ",
    "analysis/00_data_prep/01_clean_master.R to handle a mixed convention ",
    "explicitly. Offending rows:\n",
    paste(utils::capture.output(print(.chop_vaf_regression, n = Inf)),
          collapse = "\n")
  )
  stop(msg, call. = FALSE)
}

cat(sprintf("  Raw range: %.4f – %.4f (mixed scales)\n",
            min(df$vaf_raw, na.rm = TRUE), max(df$vaf_raw, na.rm = TRUE)))
cat(sprintf("  Normalized (proportion): %.6f – %.6f | Median: %.6f\n",
            min(df$vaf_prop, na.rm = TRUE), max(df$vaf_prop, na.rm = TRUE),
            median(df$vaf_prop, na.rm = TRUE)))
cat(sprintf("  Normalized (percentage): %.4f%% – %.4f%% | Median: %.4f%%\n",
            min(df$vaf_pct, na.rm = TRUE), max(df$vaf_pct, na.rm = TRUE),
            median(df$vaf_pct, na.rm = TRUE)))
cat(sprintf("  NA: %d\n", sum(is.na(df$vaf_prop))))

# Report by center after normalization
cat("  Per-center normalized (proportion):\n")
for (sc in sort(unique(df$study_code_clean[!is.na(df$vaf_prop)]))) {
  sub <- df %>% filter(study_code_clean == sc & !is.na(vaf_prop))
  sname <- study_labels[as.character(sc)]
  cat(sprintf("    %s (n=%d): %.6f – %.6f, median %.6f\n",
              sname, nrow(sub), min(sub$vaf_prop), max(sub$vaf_prop), median(sub$vaf_prop)))
}

# Flag: negative mutation with VAF (see data/flags/review_ids.csv, gitignored)
neg_with_vaf <- df %>%
  filter(mutation_code_clean %in% c("8","9","10","11","12") & !is.na(vaf_prop))
if (nrow(neg_with_vaf) > 0) {
  cat(sprintf("  ⚠ %d negative-mutation row(s) with VAF (flagged needs_review):\n", nrow(neg_with_vaf)))
  for (i in seq_len(nrow(neg_with_vaf))) {
    cat(sprintf("    Row %d: %s | mut=%s | VAF=%.4f%% | needs_review=TRUE\n",
                neg_with_vaf$row_num[i], neg_with_vaf$patient_id[i],
                neg_with_vaf$mutation_code[i], neg_with_vaf$vaf_pct[i]))
  }
}

# ── [7] sex → sex_f (labeled factor) ─────────────────────────────────────────
cat("[7/30] sex → sex_f\n")
df <- df %>%
  mutate(
    sex_f = factor(sex, levels = c(1, 2), labels = c("Female", "Male"))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$sex_f, useNA = "ifany")), collapse = "\n  ")))

# ── [8] age_dx → age_dx_numeric ──────────────────────────────────────────────
cat("[8/30] age_dx → age_dx_numeric\n")
df <- df %>%
  mutate(
    age_dx_numeric = round(as.numeric(age_dx), 2)
  )
cat(sprintf("  Range: %.2f – %.2f yrs | Median: %.2f | NA: %d\n",
            min(df$age_dx_numeric, na.rm = TRUE),
            max(df$age_dx_numeric, na.rm = TRUE),
            median(df$age_dx_numeric, na.rm = TRUE),
            sum(is.na(df$age_dx_numeric))))

# ── [9] age_surgery → age_surgery_numeric ─────────────────────────────────────
cat("[9/30] age_surgery → age_surgery_numeric\n")
df <- df %>%
  mutate(
    age_surgery_numeric = round(as.numeric(age_surgery), 2)
  )
cat(sprintf("  Range: %.2f – %.2f yrs | Median: %.2f | NA: %d\n",
            min(df$age_surgery_numeric, na.rm = TRUE),
            max(df$age_surgery_numeric, na.rm = TRUE),
            median(df$age_surgery_numeric, na.rm = TRUE),
            sum(is.na(df$age_surgery_numeric))))

# ── [9b] Consensus age: best available age per patient ──────────────────────
# Priority: age_surgery > age_dx (per Hale: "take age of surgery where present")
# Nikolaev has dx only; Priemer/Hong/Goss have surgery only; BCH/CHOP/Gao/UAB have both
cat("[9b] age (consensus), age_source\n")
df <- df %>%
  mutate(
    age = case_when(
      !is.na(age_surgery_numeric) ~ age_surgery_numeric,
      !is.na(age_dx_numeric) ~ age_dx_numeric,
      TRUE ~ NA_real_
    ),
    age_source = case_when(
      !is.na(age_surgery_numeric) ~ "surgery",
      !is.na(age_dx_numeric) ~ "dx",
      TRUE ~ NA_character_
    )
  )

n_dx   <- sum(df$age_source == "dx", na.rm = TRUE)
n_surg <- sum(df$age_source == "surgery", na.rm = TRUE)
n_na   <- sum(is.na(df$age))
cat(sprintf("  Source: dx=%d, surgery=%d, NA=%d\n", n_dx, n_surg, n_na))
cat(sprintf("  Range: %.2f – %.2f yrs | Median: %.2f\n",
            min(df$age, na.rm = TRUE), max(df$age, na.rm = TRUE),
            median(df$age, na.rm = TRUE)))

# Report by study
cat("  Per-study age source:\n")
age_by_study <- df %>%
  group_by(study_clean) %>%
  summarise(n = n(),
            from_dx = sum(age_source == "dx", na.rm = TRUE),
            from_surg = sum(age_source == "surgery", na.rm = TRUE),
            missing = sum(is.na(age)), .groups = "drop")
for (i in seq_len(nrow(age_by_study))) {
  cat(sprintf("    %-10s n=%3d | dx=%3d surg=%3d NA=%d\n",
              age_by_study$study_clean[i], age_by_study$n[i],
              age_by_study$from_dx[i], age_by_study$from_surg[i], age_by_study$missing[i]))
}

# ── [9c] dx_to_surgery_yrs: time between diagnosis and surgery ──────────────
cat("[9c] dx_to_surgery_yrs\n")
df <- df %>%
  mutate(
    dx_to_surgery_yrs = round(age_surgery_numeric - age_dx_numeric, 2)
  )
n_valid <- sum(!is.na(df$dx_to_surgery_yrs))
cat(sprintf("  n=%d | Range: %.2f – %.2f | Median: %.2f | NA: %d\n",
            n_valid,
            min(df$dx_to_surgery_yrs, na.rm = TRUE),
            max(df$dx_to_surgery_yrs, na.rm = TRUE),
            median(df$dx_to_surgery_yrs, na.rm = TRUE),
            sum(is.na(df$dx_to_surgery_yrs))))

# ── [10] race → race_f (labeled factor) ──────────────────────────────────────
cat("[10/30] race → race_f\n")
race_labels <- c(
  "1" = "White", "2" = "Black/African American", "3" = "Hispanic/Latino",
  "4" = "Asian", "5" = "Other", "6" = "Unable to collect"
)
df <- df %>%
  mutate(
    race_f = factor(race_labels[as.character(race)],
                    levels = unname(race_labels))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$race_f, useNA = "ifany")), collapse = "\n  ")))

# Binary indicator columns for each race category
df <- df %>%
  mutate(
    race_white             = as.integer(race == 1),
    race_black             = as.integer(race == 2),
    race_hispanic          = as.integer(race == 3),
    race_asian             = as.integer(race == 4),
    race_other             = as.integer(race == 5),
    race_unable_to_collect = as.integer(race == 6)
  )
cat("  Race binary indicators: race_white, race_black, race_hispanic, race_asian, race_other, race_unable_to_collect\n")
cat(sprintf("    White=%d, Black=%d, Hispanic=%d, Asian=%d, Other=%d, Unable=%d, NA=%d\n",
    sum(df$race_white == 1, na.rm=TRUE), sum(df$race_black == 1, na.rm=TRUE),
    sum(df$race_hispanic == 1, na.rm=TRUE), sum(df$race_asian == 1, na.rm=TRUE),
    sum(df$race_other == 1, na.rm=TRUE), sum(df$race_unable_to_collect == 1, na.rm=TRUE),
    sum(is.na(df$race_white))))

# ── [11] sm_size → sm_size_f (labeled factor + numeric) ──────────────────────
cat("[11/30] sm_size → sm_size_num, sm_size_f\n")
df <- df %>%
  mutate(
    sm_size_num = as.integer(sm_size),
    sm_size_f = factor(sm_size_num, levels = 1:3,
                       labels = c("<3 cm", "3-6 cm", ">6 cm"))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$sm_size_f, useNA = "ifany")), collapse = "\n  ")))

# ── [12] sm_eloquence → sm_eloquence_num (binary) ────────────────────────────
cat("[12/30] sm_eloquence → sm_eloquence_num\n")
df <- df %>%
  mutate(sm_eloquence_num = as.integer(sm_eloquence))
cat(sprintf("  0 (non-eloquent): %d | 1 (eloquent): %d | NA: %d\n",
            sum(df$sm_eloquence_num == 0, na.rm = TRUE),
            sum(df$sm_eloquence_num == 1, na.rm = TRUE),
            sum(is.na(df$sm_eloquence_num))))

# ── [13] sm_drainage → sm_drainage_num (binary) ──────────────────────────────
cat("[13/30] sm_drainage → sm_drainage_num\n")
df <- df %>%
  mutate(sm_drainage_num = as.integer(sm_drainage))
cat(sprintf("  0 (superficial): %d | 1 (deep): %d | NA: %d\n",
            sum(df$sm_drainage_num == 0, na.rm = TRUE),
            sum(df$sm_drainage_num == 1, na.rm = TRUE),
            sum(is.na(df$sm_drainage_num))))

# ── [14] sm_total → sm_total_num + sm_grade (Roman numeral factor) ───────────
cat("[14/30] sm_total → sm_total_num, sm_grade\n")
df <- df %>%
  mutate(
    sm_total_num = as.integer(sm_total),
    sm_grade = factor(sm_total_num, levels = 1:5,
                      labels = c("I", "II", "III", "IV", "V"))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$sm_grade, useNA = "ifany")), collapse = "\n  ")))

# Validate: sm_total should equal sm_size + sm_eloquence + sm_drainage
df <- df %>%
  mutate(
    sm_computed = sm_size_num + sm_eloquence_num + sm_drainage_num,
    sm_mismatch = !is.na(sm_total_num) & !is.na(sm_computed) & (sm_total_num != sm_computed)
  )
n_mismatch <- sum(df$sm_mismatch, na.rm = TRUE)
if (n_mismatch > 0) {
  cat(sprintf("  ⚠ WARNING: %d rows where SM total ≠ size + eloquence + drainage\n", n_mismatch))
} else {
  cat("  ✓ SM total validated: all match size + eloquence + drainage\n")
}

# Flag SM mismatches in needs_review with source-specific notes
# Literature rows: trust published total; our cohort (FFPE/Tissue): needs chart re-review
if (n_mismatch > 0) {
  df <- df %>%
    mutate(
      .sm_note = case_when(
        sm_mismatch & sample_type_f == "Literature" ~
          sprintf("SM mismatch (published): reported=%d computed=%d (sz=%d elo=%d dr=%d); trust published total",
                  sm_total_num, sm_computed, sm_size_num, sm_eloquence_num, sm_drainage_num),
        sm_mismatch ~
          sprintf("SM mismatch (our cohort): reported=%d computed=%d (sz=%d elo=%d dr=%d); needs chart re-review",
                  sm_total_num, sm_computed, sm_size_num, sm_eloquence_num, sm_drainage_num),
        TRUE ~ NA_character_
      ),
      review_note = case_when(
        !is.na(review_note) & !is.na(.sm_note) ~
          paste0(review_note, "; ", .sm_note),
        is.na(review_note) & !is.na(.sm_note) ~
          .sm_note,
        TRUE ~ review_note
      ),
      needs_review = !is.na(review_note)
    ) %>%
    select(-.sm_note)

  mm_rows <- df %>% filter(sm_mismatch)
  for (i in seq_len(nrow(mm_rows))) {
    cat(sprintf("    Row %3d: %s — %s\n",
                mm_rows$row_num[i], mm_rows$patient_id[i],
                ifelse(mm_rows$sample_type_f[i] == "Literature", "Literature (trust published)", "Our cohort (needs re-review)")))
  }
}

# ── [15] compact_nidus → compact_nidus_num (binary 0/1) ───────────────────────
# Code 2 = "Not available" → recode to NA to keep column strictly binary
cat("[15/30] compact_nidus → compact_nidus_num\n")
df <- df %>%
  mutate(
    compact_nidus_num = case_when(
      compact_nidus %in% c("0", "1") ~ as.integer(compact_nidus),
      TRUE ~ NA_integer_
    )
  )
n_2_cn <- sum(df$compact_nidus == "2", na.rm = TRUE)
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d (includes %d recoded from code 2 'not available')\n",
            sum(df$compact_nidus_num == 0, na.rm = TRUE),
            sum(df$compact_nidus_num == 1, na.rm = TRUE),
            sum(is.na(df$compact_nidus_num)), n_2_cn))

# ── [16] ruptured_at_surgery → ruptured_at_surgery_num (binary) ──────────────
cat("[16/30] ruptured_at_surgery → ruptured_at_surgery_num\n")
df <- df %>%
  mutate(ruptured_at_surgery_num = as.integer(ruptured_at_surgery))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$ruptured_at_surgery_num == 0, na.rm = TRUE),
            sum(df$ruptured_at_surgery_num == 1, na.rm = TRUE),
            sum(is.na(df$ruptured_at_surgery_num))))

# ── [17] ever_ruptured → ever_ruptured_num (binary 0/1) ────────────────────────
# Code 2 (1 row, UAB) is undefined for this binary column → recode to NA
cat("[17/30] ever_ruptured → ever_ruptured_num\n")
df <- df %>%
  mutate(
    ever_ruptured_num = case_when(
      ever_ruptured %in% c("0", "1") ~ as.integer(ever_ruptured),
      TRUE ~ NA_integer_
    )
  )
n_2_er <- sum(df$ever_ruptured == "2", na.rm = TRUE)
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d (includes %d recoded from code 2)\n",
            sum(df$ever_ruptured_num == 0, na.rm = TRUE),
            sum(df$ever_ruptured_num == 1, na.rm = TRUE),
            sum(is.na(df$ever_ruptured_num)), n_2_er))

# ── [17b] Rupture-flag reconciliation (DATA_DECISIONS §42; Hale v0 notes L25) ──
# Hale's coding rule:
#   "If ruptured_at_surgery is 1, then ever_ruptured should be 1."
# `ever_ruptured` is collected as LIFETIME rupture (inclusive of the acute
# presenting bleed), NOT the "prior to this surgery only" reading once floated
# in §23 — that strictly-prior concept is the DERIVED `previous_rupture`
# variable built in 02_prep_analysis_dataset.R, not `ever_ruptured`. Cohort-wide
# evidence agrees: 211/216 (98%) of rupt_surg=1 patients already carry ever=1.
#
# We therefore apply the implication as a RULE rather than a row-specific patch:
# any row with ruptured_at_surgery_num == 1 whose ever_ruptured_num is not
# already 1 — i.e. it is 0, or NA after the code-2 recode in [17] — is
# reconciled to 1. This supersedes §23 ("flag only; do not override"). Each
# change is recorded in pub_data_note. On the 18 May 26 master this reconciles
# 5 BCH rows (ever=0) + 1 UAB row (ever was code 2 → NA): BCH100/153/190/197/214
# and NP-24-37XXX.
.rupt_needs_fix <- which(df$ruptured_at_surgery_num == 1 &
                          !is.na(df$ruptured_at_surgery_num) &
                          (df$ever_ruptured_num == 0 | is.na(df$ever_ruptured_num)))
if (length(.rupt_needs_fix) > 0) {
  .n_from_0  <- sum(df$ever_ruptured_num[.rupt_needs_fix] == 0, na.rm = TRUE)
  .n_from_na <- sum(is.na(df$ever_ruptured_num[.rupt_needs_fix]))
  # NB: .pub_note vector was rm()'d after being copied into df$pub_data_note
  # earlier in the script. Update df$pub_data_note directly here.
  note_txt <- "ever_ruptured set to 1 (ruptured_at_surgery=1 ⇒ ever_ruptured=1; lifetime coding per Hale v0 L25 / DATA_DECISIONS §42)"
  for (i in .rupt_needs_fix) {
    df$pub_data_note[i] <- ifelse(is.na(df$pub_data_note[i]), note_txt,
                                   paste(df$pub_data_note[i], note_txt, sep = "; "))
    cat(sprintf("  [17b] Rupture reconcile: %s — ever_ruptured %s→1 (rupt_surg=1)\n",
                df$patient_id[i],
                ifelse(is.na(df$ever_ruptured_num[i]), "NA", "0")))
  }
  df$ever_ruptured_num[.rupt_needs_fix] <- 1L
  cat(sprintf("  → %d row(s) reconciled to ever_ruptured=1 (%d from 0, %d from NA/code-2) per §42\n",
              length(.rupt_needs_fix), .n_from_0, .n_from_na))
}

# ── [18] growing → growing_num (binary 0/1) ────────────────────────────────────
# Code 2 (145 rows) = "Not assessed" → recode to NA to keep column strictly binary
cat("[18/30] growing → growing_num\n")
df <- df %>%
  mutate(
    growing_num = case_when(
      growing %in% c("0", "1") ~ as.integer(growing),
      TRUE ~ NA_integer_
    )
  )
n_2_gr <- sum(df$growing == "2", na.rm = TRUE)
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d (includes %d recoded from code 2 'not assessed')\n",
            sum(df$growing_num == 0, na.rm = TRUE),
            sum(df$growing_num == 1, na.rm = TRUE),
            sum(is.na(df$growing_num)), n_2_gr))

# Flag residual clinical-column issues in needs_review.
# The former "rupt_surg=1 but ever=0" contradiction flag (the §23 "prior rupture
# only" reading) is RETIRED: those rows are now reconciled to ever=1 in [17b]
# per §42, so the contradiction cannot reach downstream. We keep a transparency
# note for any row whose ever_ruptured was the undefined code 2 and was therefore
# imputed (rather than read directly) when reconciled in [17b].
df <- df %>%
  mutate(
    .clin_note = case_when(
      ever_ruptured == "2" ~
        "ever_ruptured original code 2 (undefined); reconciled to 1 because ruptured_at_surgery=1 (§42)",
      TRUE ~ NA_character_
    ),
    review_note = case_when(
      !is.na(review_note) & !is.na(.clin_note) ~
        paste0(review_note, "; ", .clin_note),
      is.na(review_note) & !is.na(.clin_note) ~
        .clin_note,
      TRUE ~ review_note
    ),
    needs_review = !is.na(review_note)
  ) %>%
  select(-.clin_note)

# Post-condition: the §42 rule must leave zero rupt_surg=1 & ever=0 rows.
n_er2 <- sum(df$ever_ruptured == "2", na.rm = TRUE)
n_rupt_contra <- sum(df$ruptured_at_surgery_num == 1 & df$ever_ruptured_num == 0, na.rm = TRUE)
cat(sprintf("  ever_ruptured code-2 rows reconciled & flagged: %d | residual rupt contradictions (must be 0): %d\n",
            n_er2, n_rupt_contra))
if (n_rupt_contra > 0) warning("Rupture invariant violated: ", n_rupt_contra,
                               " rows still have rupt_surg=1 & ever=0 after §42 reconciliation")

# ── [19] location_codes → one-hot location columns ──────────────────────────
cat("[19/30] location_codes → individual location flags\n")
location_map <- c(
  "0"  = "frontal",
  "1"  = "temporal",
  "3"  = "insular",
  "4"  = "basal_ganglia",
  "5"  = "thalamus",
  "6"  = "periventricular",
  "7"  = "parietal",
  "8"  = "occipital",
  "9"  = "cerebellar",
  "10" = "brainstem",
  "11" = "corpus_callosum",
  "12" = "cingulate",
  "13" = "sylvian_fissure",
  "14" = "dura"
)

# Initialize location columns
for (loc_name in location_map) {
  df[[paste0("loc_", loc_name)]] <- 0L
}

# Parse comma-separated codes
for (i in seq_len(nrow(df))) {
  codes_str <- df$location_codes[i]
  if (is.na(codes_str)) next
  codes <- trimws(unlist(strsplit(codes_str, ",")))
  for (code in codes) {
    loc_name <- location_map[code]
    if (!is.na(loc_name)) {
      df[[paste0("loc_", loc_name)]][i] <- 1L
    }
  }
}

# Count locations per patient
df <- df %>%
  mutate(
    n_locations = rowSums(across(starts_with("loc_")), na.rm = TRUE),
    n_locations = ifelse(is.na(location_codes), NA_integer_, n_locations)
  )

loc_totals <- df %>%
  summarise(across(starts_with("loc_"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "location", values_to = "n") %>%
  arrange(desc(n))
cat("  Location frequencies:\n")
for (j in seq_len(nrow(loc_totals))) {
  cat(sprintf("    %s: %d\n", gsub("loc_", "", loc_totals$location[j]), loc_totals$n[j]))
}
cat(sprintf("  Multi-location AVMs: %d\n", sum(df$n_locations > 1, na.rm = TRUE)))

# ── [20] laterality → laterality_f ──────────────────────────────────────────
cat("[20/30] laterality → laterality_num, laterality_f\n")
# Data dictionary: R=1, L=2, bihemispheric/midline=3; data also has "0"
df <- df %>%
  mutate(
    laterality_num = as.integer(laterality),
    laterality_f = factor(laterality_num, levels = c(0, 1, 2, 3),
                          labels = c("Unknown/Not specified", "Right", "Left",
                                     "Bihemispheric/Midline"))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$laterality_f, useNA = "ifany")), collapse = "\n  ")))

# ── [21] prior_seizure → prior_seizure_num (binary) ─────────────────────────
cat("[21/30] prior_seizure → prior_seizure_num\n")
df <- df %>%
  mutate(prior_seizure_num = as.integer(prior_seizure))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$prior_seizure_num == 0, na.rm = TRUE),
            sum(df$prior_seizure_num == 1, na.rm = TRUE),
            sum(is.na(df$prior_seizure_num))))

# ── [22] prior_radiation → prior_radiation_num (binary) ─────────────────────
cat("[22/30] prior_radiation → prior_radiation_num\n")
df <- df %>%
  mutate(prior_radiation_num = as.integer(prior_radiation))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$prior_radiation_num == 0, na.rm = TRUE),
            sum(df$prior_radiation_num == 1, na.rm = TRUE),
            sum(is.na(df$prior_radiation_num))))

# ── [23] radiation_shrink → radiation_shrink_f ───────────────────────────────
cat("[23/30] radiation_shrink → radiation_shrink_num, radiation_shrink_f\n")
df <- df %>%
  mutate(
    radiation_shrink_num = as.integer(radiation_shrink),
    radiation_shrink_f = factor(radiation_shrink_num, levels = c(0, 1, 3),
                                labels = c("No", "Yes", "N/A (no radiation)"))
  )
cat(sprintf("  %s\n", paste(capture.output(table(df$radiation_shrink_f, useNA = "ifany")), collapse = "\n  ")))

# ── [24] prior_embolization → prior_embolization_num (binary) ────────────────
cat("[24/30] prior_embolization → prior_embolization_num\n")
df <- df %>%
  mutate(prior_embolization_num = as.integer(prior_embolization))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$prior_embolization_num == 0, na.rm = TRUE),
            sum(df$prior_embolization_num == 1, na.rm = TRUE),
            sum(is.na(df$prior_embolization_num))))

# ── [25] intranidal_aneurysm → intranidal_aneurysm_num (binary) ─────────────
cat("[25/30] intranidal_aneurysm → intranidal_aneurysm_num\n")
df <- df %>%
  mutate(intranidal_aneurysm_num = as.integer(intranidal_aneurysm))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$intranidal_aneurysm_num == 0, na.rm = TRUE),
            sum(df$intranidal_aneurysm_num == 1, na.rm = TRUE),
            sum(is.na(df$intranidal_aneurysm_num))))

# ── [26] venous_varix → venous_varix_num (binary) ───────────────────────────
cat("[26/30] venous_varix → venous_varix_num\n")
df <- df %>%
  mutate(venous_varix_num = as.integer(venous_varix))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$venous_varix_num == 0, na.rm = TRUE),
            sum(df$venous_varix_num == 1, na.rm = TRUE),
            sum(is.na(df$venous_varix_num))))

# ── [27] venous_outflow_stenosis → venous_outflow_stenosis_num (binary) ─────
cat("[27/30] venous_outflow_stenosis → venous_outflow_stenosis_num\n")
df <- df %>%
  mutate(venous_outflow_stenosis_num = as.integer(venous_outflow_stenosis))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$venous_outflow_stenosis_num == 0, na.rm = TRUE),
            sum(df$venous_outflow_stenosis_num == 1, na.rm = TRUE),
            sum(is.na(df$venous_outflow_stenosis_num))))

# ── [28] flow_related_aneurysm → flow_related_aneurysm_num (binary) ────────
cat("[28/30] flow_related_aneurysm → flow_related_aneurysm_num\n")
df <- df %>%
  mutate(flow_related_aneurysm_num = as.integer(flow_related_aneurysm))
cat(sprintf("  0 (no): %d | 1 (yes): %d | NA: %d\n",
            sum(df$flow_related_aneurysm_num == 0, na.rm = TRUE),
            sum(df$flow_related_aneurysm_num == 1, na.rm = TRUE),
            sum(is.na(df$flow_related_aneurysm_num))))

# ── [29] n_high_risk_features → n_high_risk_num ─────────────────────────────
cat("[29/30] n_high_risk_features → n_high_risk_num\n")
df <- df %>%
  mutate(n_high_risk_num = as.integer(n_high_risk_features))

# Fix 1: Gao uses binary 1=yes/2=no for "has associated aneurysm" — recode 2→0
n_gao_recode <- sum(df$n_high_risk_num == 2 & df$study_clean == "Gao", na.rm = TRUE)
df <- df %>%
  mutate(n_high_risk_num = ifelse(study_clean == "Gao" & n_high_risk_num == 2,
                                  0L, n_high_risk_num))
cat(sprintf("  Gao binary recode: %d rows 2→0 (1=yes/2=no → 1/0 count)\n", n_gao_recode))

# Validate: should equal sum of cols 25-28 (intranidal + varix + stenosis + flow)
df <- df %>%
  mutate(
    high_risk_computed = intranidal_aneurysm_num + venous_varix_num +
                         venous_outflow_stenosis_num + flow_related_aneurysm_num,
    high_risk_mismatch = !is.na(n_high_risk_num) & !is.na(high_risk_computed) &
                         (n_high_risk_num != high_risk_computed)
  )
n_hr_mismatch <- sum(df$high_risk_mismatch, na.rm = TRUE)
cat(sprintf("  Range: %d – %d | NA: %d\n",
            min(df$n_high_risk_num, na.rm = TRUE),
            max(df$n_high_risk_num, na.rm = TRUE),
            sum(is.na(df$n_high_risk_num))))
if (n_hr_mismatch > 0) {
  cat(sprintf("  ⚠ WARNING: %d rows where # high risk ≠ sum of cols 25-28\n", n_hr_mismatch))
} else {
  cat("  ✓ High risk count validated against component columns\n")
}

# Fix 2: Fill n_high_risk from component sum where all 4 are non-NA but count is NA
fill_idx <- which(is.na(df$n_high_risk_num) & !is.na(df$high_risk_computed))
n_hr_filled <- length(fill_idx)
if (n_hr_filled > 0) {
  df$n_high_risk_num[fill_idx] <- df$high_risk_computed[fill_idx]
  cat(sprintf("  Computed n_high_risk from component sum: filled %d rows\n", n_hr_filled))
} else {
  cat("  No rows needed n_high_risk computed from components\n")
}

# ── Genotype-only "phenotype_pending" flag (DATA_DECISIONS §43) ─────────────
# Rows that have mutation+VAF data but lack demographics, SM grading, and
# clinical phenotype (the 3 CHOP orphan rows in the 18 May 26 master:
# 1927848, 1923749 ([specimen accession]), VCHOPUC066). The flag lets downstream
# producers include these rows in raw genotype counts but auto-exclude
# them from phenotype-correlation analyses (KM curves, SM-by-genotype
# forests, etc.) without per-script special casing.
df <- df %>%
  mutate(
    phenotype_pending = !is.na(mutation_code_clean) &
                         is.na(age_dx_numeric) &
                         is.na(age_surgery_numeric) &
                         is.na(sm_total_num) &
                         is.na(ever_ruptured_num)
  )
.n_pending <- sum(df$phenotype_pending, na.rm = TRUE)
if (.n_pending > 0) {
  cat(sprintf("  phenotype_pending = TRUE: %d row(s) (genotype-only; see DATA_DECISIONS §43)\n",
              .n_pending))
  .pending_rows <- df %>% filter(phenotype_pending)
  for (i in seq_len(nrow(.pending_rows))) {
    cat(sprintf("    %s — mutation=%s, vaf=%s\n",
                .pending_rows$patient_id[i],
                as.character(.pending_rows$mutation_code_clean[i]),
                as.character(.pending_rows$vaf_raw[i])))
  }
}

# ── [30] notes — kept as-is ──────────────────────────────────────────────────
cat("[30/30] notes: free text — no transformation\n")
cat(sprintf("  %d rows with notes\n", sum(!is.na(df$notes))))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Select and order final clean columns
# ══════════════════════════════════════════════════════════════════════════════
cat("\n═══ Assembling clean dataset ═══\n")

df_clean <- df %>%
  select(
    # Identifiers (original + clean)
    row_num, udn, sample_uid, sample_num,
    multi_sample_flag, multi_sample_class, lesion_id, is_event_primary_row,
    spatial_sample, phenotype_pending,
    needs_review, review_note, pub_data_note,
    mrn, patient_id, patient_id_clean,
    study_code, study_code_clean, study, study_clean,
    sample_type, sample_type_f,

    # Genotype
    mutation_code, mutation_code_clean, mutation, mutation_positive,
    mutation_gene, tissue_available_not_tested,
    starts_with("mut_"),
    vaf_raw, vaf_prop, vaf_pct,

    # Demographics
    sex, sex_f, age_dx_numeric, age_surgery_numeric, age, age_source, dx_to_surgery_yrs,
    race, race_f, starts_with("race_"),

    # Spetzler-Martin grading
    sm_size_num, sm_size_f, sm_eloquence_num, sm_drainage_num,
    sm_total_num, sm_grade, sm_computed, sm_mismatch,

    # Nidus characteristics
    compact_nidus_num,

    # Clinical history
    ruptured_at_surgery_num, ever_ruptured_num,
    growing_num,

    # Location
    location_codes, starts_with("loc_"), n_locations,
    laterality_num, laterality_f,

    # Prior treatments
    prior_seizure_num, prior_radiation_num,
    radiation_shrink_num, radiation_shrink_f,
    prior_embolization_num,

    # High-risk features
    intranidal_aneurysm_num, venous_varix_num,
    venous_outflow_stenosis_num, flow_related_aneurysm_num,
    n_high_risk_num, high_risk_computed, high_risk_mismatch,

    # Notes
    notes
  )

cat(sprintf("Clean dataset: %d rows × %d columns\n", nrow(df_clean), ncol(df_clean)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Summary report
# ══════════════════════════════════════════════════════════════════════════════
cat("\n═══ Missing Data Summary ═══\n")
missing_summary <- tibble(
  column = names(df_clean),
  n_missing = sapply(df_clean, function(x) sum(is.na(x))),
  pct_missing = round(100 * n_missing / nrow(df_clean), 1)
) %>%
  filter(n_missing > 0) %>%
  arrange(desc(pct_missing))

for (i in seq_len(nrow(missing_summary))) {
  cat(sprintf("  %-35s %3d (%4.1f%%)\n",
              missing_summary$column[i],
              missing_summary$n_missing[i],
              missing_summary$pct_missing[i]))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Save
# ══════════════════════════════════════════════════════════════════════════════
out_dir <- here::here("data", "processed")
saveRDS(df_clean, file.path(out_dir, "bAVM_genopheno_clean.rds"))
writexl::write_xlsx(df_clean, file.path(out_dir, "bAVM_genopheno_clean.xlsx"))

cat(sprintf("\n═══ Done ═══\n"))
cat(sprintf("Saved: %s\n", file.path(out_dir, "bAVM_genopheno_clean.rds")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "bAVM_genopheno_clean.xlsx")))
