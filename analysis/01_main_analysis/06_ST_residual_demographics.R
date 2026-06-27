# 06_ST_residual_demographics.R
# ─────────────────────────────────────────────────────────────────────────────
# Computes the handful of prose-level stats that aren't produced by any
# other figure script: anatomical coarse-grouping tests (left-vs-right
# hemisphere, 3-level laterality, supratentorial-vs-infratentorial),
# sex and race Fisher tests, the age-by-phenotype interaction min-P
# read off Supp Table 9, and the VAF-nested-model F-test reported in
# Results §4. Dumps everything as one fragment ("demographics") so
# the manifest builder folds it into the manuscript stats.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(here); library(readxl); library(writexl)
})
# 2026-06-14 fix: utils.R was never sourced, so the MASTER_SEED reference
# inside the race Fisher tryCatch threw "object not found" and the
# tryCatch silently returned NA_real_ for race_p (which then surfaced as
# a blank cell in SuppTable10).
source(here("analysis", "helper_scripts", "utils.R"))

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(geno_binary))

# ── 1. Hemisphere (Left vs Right, drops Unknown + Bihemispheric) ────────────
hem_df <- df %>%
  filter(laterality_f %in% c("Left", "Right"))
hem_p <- fisher.test(table(hem_df$laterality_f, hem_df$geno_binary))$p.value
hem_n <- nrow(hem_df)

# ── 2. Three-level laterality (L / R / Bihemispheric) ───────────────────────
lat3_df <- df %>%
  filter(laterality_f %in% c("Left", "Right", "Bihemispheric/Midline"))
lat3_p <- fisher.test(table(lat3_df$laterality_f, lat3_df$geno_binary))$p.value
lat3_n <- nrow(lat3_df)

# ── 3. Supratentorial vs infratentorial ─────────────────────────────────────
# Any supratentorial lobe flagged -> supra; otherwise any infra flag -> infra.
# Lesions with both are supra (dominant coverage); lesions with neither drop.
supra_cols <- c("loc_frontal", "loc_temporal", "loc_parietal", "loc_occipital",
                "loc_insular", "loc_basal_ganglia", "loc_thalamus",
                "loc_periventricular", "loc_corpus_callosum", "loc_cingulate",
                "loc_sylvian_fissure")
infra_cols <- c("loc_cerebellar", "loc_brainstem")

tent_df <- df %>%
  mutate(
    any_supra = rowSums(across(all_of(supra_cols)), na.rm = TRUE) > 0,
    any_infra = rowSums(across(all_of(infra_cols)), na.rm = TRUE) > 0,
    region = case_when(
      any_supra             ~ "Supratentorial",
      !any_supra & any_infra ~ "Infratentorial",
      TRUE                  ~ NA_character_
    )
  ) %>%
  filter(!is.na(region))
tent_p <- fisher.test(table(tent_df$region, tent_df$geno_binary))$p.value
tent_n <- nrow(tent_df)

# ── 4. Sex (male fraction by genotype) ──────────────────────────────────────
sex_df <- df %>% filter(!is.na(sex_f), sex_f %in% c("Male", "Female"))
sex_tab <- table(sex_df$sex_f, sex_df$geno_binary)
sex_p <- fisher.test(sex_tab)$p.value
sex_n <- sum(sex_tab)
male_pct_mut <- 100 * sex_tab["Male", "Variant-positive"] /
  sum(sex_tab[, "Variant-positive"])
male_pct_neg <- 100 * sex_tab["Male", "Panel-negative"] /
  sum(sex_tab[, "Panel-negative"])

# ── 5. Race (Fisher among lesions with recorded race) ───────────────────────
race_df <- df %>% filter(!is.na(race_f), race_f != "Unknown/Not specified")
race_p <- tryCatch({
  set.seed(MASTER_SEED)   # audit F13: canonical seed from utils.R
  fisher.test(table(race_df$race_f, race_df$geno_binary),
              simulate.p.value = TRUE, B = 20000)$p.value
}, error = function(e) NA_real_)
race_n <- nrow(race_df)
# Per-arm race ascertainment denominators (audit C8): asymmetric
# missingness can bias the Fisher test, so the prose / Methods should
# disclose how many of each genotype arm contributed.
race_n_mut    <- sum(race_df$geno_binary == "Variant-positive", na.rm = TRUE)
race_n_neg    <- sum(race_df$geno_binary == "Panel-negative", na.rm = TRUE)
race_n_total  <- sum(!is.na(df$geno_binary))
race_pct_mut  <- 100 * race_n_mut /
  sum(df$geno_binary == "Variant-positive", na.rm = TRUE)
race_pct_neg  <- 100 * race_n_neg /
  sum(df$geno_binary == "Panel-negative", na.rm = TRUE)

# ── 6. Age × phenotype interaction min P (from Supp Table 9 producer) ───────
# The interaction tests are computed in 20_ST_supp_tables_8_9_10.R and written
# to results/SupplementaryTables/SuppTable04_age_phenotype_interactions.xlsx.
#
# 2026-06-14: ST4 is now a single mega-sheet with an Analysis column
# distinguishing the Cox model (HR-coefficient rows) from the two
# interaction framings (β-coefficient rows). Filter by Analysis prefix
# "Phenotype ←" / "Age ←" to pick out interaction rows. Tolerates the
# legacy two-sheet schema so a rerun in the wrong order doesn't silently
# NA the prose stat.
age_int_path <- here("results", "SupplementaryTables",
                     "SuppTable04_age_phenotype_interactions.xlsx")
.age_int <- tryCatch({
  shts <- readxl::excel_sheets(age_int_path)
  to_num <- function(x) suppressWarnings(as.numeric(x))
  if (length(shts) == 1L) {
    # New single-sheet schema: one sheet with an Analysis column.
    rows <- read_excel(age_int_path, sheet = 1, col_types = "text")
    rows <- rows[grepl("^Phenotype ←|^Age ←", rows$Analysis), , drop = FALSE]
  } else {
    # Legacy multi-sheet schema (pre-2026-06-14): "... Geno x ..." per sheet.
    shts <- shts[grepl("Geno", shts)]
    rows <- do.call(rbind, lapply(shts, function(s)
      read_excel(age_int_path, sheet = s, col_types = "text")))
  }
  p_col   <- intersect(c("P",     "interaction_p"),     names(rows))[1]
  fdr_col <- intersect(c("FDR P", "interaction_p_fdr"), names(rows))[1]
  list(p = min(to_num(rows[[p_col]]),   na.rm = TRUE),
       q = min(to_num(rows[[fdr_col]]), na.rm = TRUE))
}, error = function(e) list(p = NA_real_, q = NA_real_))
age_pheno_interaction_min_p <- .age_int$p
age_pheno_interaction_min_q <- .age_int$q

# ── 7. VAF nested-model F-test (Results §4 VAF dose-response) ───────────────
# Two specifications, both surfaced (audit 12 found a prose↔producer
# mismatch where §4 said "binary genotype" but the producer fit was on
# variant_group = G12D vs G12V):
#
#   (a) "Binary genotype": age ~ geno_binary  vs  age ~ geno_binary + vaf_pct
#       Negatives carry vaf_pct = 0 by construction. This is the test the
#       Results prose actually claims, and gives the most generalisable
#       answer ("does VAF add information beyond presence/absence of any
#       mutation?").
#
#   (b) Within-KRAS variant-resolved: age ~ variant_group  vs  + vaf_pct
#       Restricted to G12D + G12V Mut+ lesions (other KRAS / BRAF strata
#       have VAF n < 10). Tests whether VAF adds information beyond
#       knowing G12D-vs-G12V.
#
# We keep both because they answer different questions; the prose at
# §4 ¶6 should cite the binary version (vaf_nested_p) to match its
# "binary genotype" wording.

# (a) Binary genotype LHS
vaf_age_binary_df <- df %>%
  filter(!is.na(age), !is.na(geno_binary)) %>%
  mutate(vaf_pct = ifelse(is.na(vaf_prop), 0, vaf_prop * 100))
m1_bin <- lm(age ~ geno_binary,              data = vaf_age_binary_df)
m2_bin <- lm(age ~ geno_binary + vaf_pct,    data = vaf_age_binary_df)
aov_bin <- anova(m1_bin, m2_bin)
vaf_nested_f <- aov_bin$F[2]
vaf_nested_p <- aov_bin$`Pr(>F)`[2]

# (b) Within-KRAS variant-resolved LHS (audit-secondary, surfaced for
#     reproducibility of the previously-cited number)
vaf_age_df <- df %>%
  filter(!is.na(age), !is.na(vaf_prop),
         geno_binary == "Variant-positive",
         mutation %in% c("KRAS G12D", "KRAS G12V")) %>%
  mutate(vaf_pct = vaf_prop * 100,
         variant_group = factor(mutation))
m1 <- lm(age ~ variant_group,              data = vaf_age_df)
m2 <- lm(age ~ variant_group + vaf_pct,    data = vaf_age_df)
aov_comp <- anova(m1, m2)
vaf_nested_within_kras_f <- aov_comp$F[2]
vaf_nested_within_kras_p <- aov_comp$`Pr(>F)`[2]

# ── emit supplementary table: coarse anatomical + demographic null tests ───
fmt_p_disp <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}
source(here("analysis", "helper_scripts", "supp_table_writer.R"))
.fmt_p_st10 <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001,
                formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3)))
}
.race_p_vec <- c(hem_p, lat3_p, tent_p, sex_p, race_p)
st10_df <- tibble::tibble(
  variable   = c("Hemisphere (Left vs Right)",
                 "Laterality (Left / Right / Bihemispheric–Midline)",
                 "Supratentorial vs Infratentorial",
                 "Sex (Male vs Female)",
                 "Self-reported race"),
  comparison = rep("Variant-positive vs Panel-negative", 5),
  n          = as.integer(c(hem_n, lat3_n, tent_n, sex_n, race_n)),
  test       = c("Fisher exact (2 × 2)",
                 "Fisher exact (3 × 2)",
                 "Fisher exact (2 × 2)",
                 "Fisher exact (2 × 2)",
                 "Fisher exact (simulated, B = 20,000)"),
  p_str      = .fmt_p_st10(.race_p_vec)
)
coarse_path <- here("results", "SupplementaryTables",
                    "SuppTable10_coarse_groupings.xlsx")
write_supp_table(
  data    = st10_df,
  path    = coarse_path,
  sheet   = "Coarse anatomical & demographic",
  columns = list(
    col    ("variable",   label = "Variable"),
    col    ("comparison", label = "Comparison"),
    col_int("n",          label = "N",    italic = TRUE),
    col    ("test",       label = "Test"),
    col    ("p_str",      label = "P",    italic = TRUE)
  ),
  footnote = c(
    "Five residual coarse-grouping null tests not covered by the per-lobe interaction table (ST09). Each compares variant-positive vs panel-negative bAVMs.",
    "P-values from Fisher's exact test; the self-reported race test uses a simulated null (B = 20,000) because the 6 × 2 contingency table has cells too sparse for the deterministic Fisher computation."
  )
)
cat(sprintf("→ %s\n", coarse_path))

# ── BH-FDR across the 5 demographic tests (audit B4) ───────────────────────
# Audit found §4 ¶9 cited raw P, with laterality3 raw P = 0.031 sitting
# inside a "no association" frame. Surface BH-q so prose can cite the
# corrected value alongside (or in place of) raw P.
.demo_p   <- c(hemisphere = hem_p, laterality3 = lat3_p,
               supratentorial = tent_p, sex = sex_p, race = race_p)
.demo_q   <- p.adjust(.demo_p, method = "BH")

# ── emit fragment ───────────────────────────────────────────────────────────
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))
write_stats_section(
  section = "demographics",
  stats = list(
    hemisphere_fisher_p         = hem_p,
    hemisphere_fisher_q         = unname(.demo_q["hemisphere"]),
    hemisphere_fisher_n         = hem_n,
    laterality3_fisher_p        = lat3_p,
    laterality3_fisher_q        = unname(.demo_q["laterality3"]),
    laterality3_fisher_n        = lat3_n,
    supratentorial_fisher_p     = tent_p,
    supratentorial_fisher_q     = unname(.demo_q["supratentorial"]),
    supratentorial_fisher_n     = tent_n,
    sex_fisher_p                = sex_p,
    sex_fisher_q                = unname(.demo_q["sex"]),
    sex_fisher_n                = sex_n,
    sex_male_pct_mut            = male_pct_mut,
    sex_male_pct_neg            = male_pct_neg,
    race_fisher_p               = race_p,
    race_fisher_q               = unname(.demo_q["race"]),
    race_fisher_n               = race_n,
    race_n_mut                  = race_n_mut,
    race_n_neg                  = race_n_neg,
    race_pct_ascertained_mut    = race_pct_mut,
    race_pct_ascertained_neg    = race_pct_neg,
    demographics_q_min          = min(.demo_q, na.rm = TRUE),
    demographics_q_max          = max(.demo_q, na.rm = TRUE),
    age_pheno_interaction_min_p = age_pheno_interaction_min_p,
    age_pheno_interaction_min_q = age_pheno_interaction_min_q,
    vaf_nested_f                = vaf_nested_f,
    vaf_nested_p                = vaf_nested_p,
    vaf_nested_within_kras_f    = vaf_nested_within_kras_f,
    vaf_nested_within_kras_p    = vaf_nested_within_kras_p
  )
)

cat(sprintf("✓ residual demographics stats written (n=%d, geno-paired)\n",
            nrow(df)))
