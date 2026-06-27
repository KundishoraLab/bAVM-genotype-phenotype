# =============================================================================
# stats_schema.R — declarative required-keys for the manuscript stats manifest
# -----------------------------------------------------------------------------
# Single source of truth for what keys every section of the manuscript stats
# manifest MUST contain. check_stats_manifest.R reads this and errors if any
# key is missing. Grows section-by-section as Phase 2 wires each analysis
# script's write_stats_section() calls.
#
# Phase 1 seed: only `cohort` and `meta` are required. Everything else is
# listed commented-out, to be uncommented as Phase 2 migrates each script.
# =============================================================================

stats_schema <- list(

  # ---- cohort-level counts (always required) ------------------------------
  cohort = c(
    "n_total_harmonised",
    "n_tissue_tested",
    "n_bch", "n_chop", "n_uab",
    "n_nikolaev", "n_priemer", "n_hong", "n_goss", "n_gao"
  ),

  # ---- meta (always required; written by build_stats_manifest.R) ---------
  meta = c("generated_at", "git_sha", "n_sections"),

  # ---- VAF × age-at-rupture (mirror of Fig 2e) + prevalence-vs-timing ----
  # Producer: analysis/01_main_analysis/15_ED7_vaf_age_at_rupture.R
  # Added 2026-04-25 per Hale comment: Andy asked whether the inverse VAF
  # × age slope mirrors at the rupture endpoint. Pooled and per-variant
  # slopes are provided plus the side-by-side prevalence (Fig 3d) versus
  # timing (Fig 2d) summary that supports the §2→§3 transition sentence.
  vaf_rupture = c(
    "n_rupt_pooled",
    "slope_pooled", "slope_ci_lo_pooled", "slope_ci_hi_pooled",
    "slope_p_pooled", "rho_pooled", "rho_p_pooled",
    "n_g12d", "slope_g12d", "slope_g12d_p",
    "n_g12v", "slope_g12v", "slope_g12v_p",
    "prev_rupt_mut_pct", "prev_rupt_neg_pct",
    "median_age_rupt_mut", "median_age_rupt_neg",
    "median_age_pres_mut", "median_age_pres_neg"
  ),

  # ---- residual demographics / coarse-anatomy tests ----------------------
  # Producer: analysis/01_main_analysis/06_ST_residual_demographics.R
  # Covers the prose-level stats that don't belong to any figure producer:
  # hemisphere / 3-level laterality / supratentorial-vs-infratentorial
  # Fisher Ps, sex and race Fisher Ps (§3), the age-by-phenotype
  # interaction min P (§4), and the VAF nested-model F-test (§4 dose-
  # response claim).
  demographics = c(
    "hemisphere_fisher_p", "hemisphere_fisher_q", "hemisphere_fisher_n",
    "laterality3_fisher_p", "laterality3_fisher_q", "laterality3_fisher_n",
    "supratentorial_fisher_p", "supratentorial_fisher_q",
    "supratentorial_fisher_n",
    "sex_fisher_p", "sex_fisher_q", "sex_fisher_n",
    "sex_male_pct_mut", "sex_male_pct_neg",
    "race_fisher_p", "race_fisher_q", "race_fisher_n",
    "race_n_mut", "race_n_neg",
    "race_pct_ascertained_mut", "race_pct_ascertained_neg",
    "demographics_q_max",
    "age_pheno_interaction_min_p",
    "age_pheno_interaction_min_q",
    "vaf_nested_f", "vaf_nested_p",
    "vaf_nested_within_kras_f", "vaf_nested_within_kras_p"
  ),

  # ---- Figure 1 (Section 1 — cohort + variant landscape + VAF) -----------
  # Producer: analysis/01_main_analysis/08_F1_cohort_variants.R
  # Scope: every number Results Section 1 could cite. Keys intentionally
  # expansive; add new ones here when prose calls for them.
  fig1 = c(
    # cohort denominators
    "cohort_n_total", "cohort_n_tissue_tested", "cohort_n_not_tested",
    "cohort_n_mut_pos", "cohort_n_mut_neg",
    "cohort_pct_tissue_tested", "cohort_mut_rate_pct",
    # variant counts (fine-grained, variant-positive subset)
    "n_kras_g12d", "n_kras_g12v", "n_kras_g12c", "n_kras_g12a",
    "n_kras_q61h", "n_kras_dup", "n_braf_v600e", "n_braf_q636x",
    "n_kras_any", "n_braf_any",
    "n_kras_codon12", "pct_kras_codon12_of_mutpos",
    "n_institutional_cohort", "n_published_cohort",
    # variant counts (geno_variant 5-level factor)
    "n_g12d_grp", "n_g12v_grp", "n_other_kras_grp",
    "n_braf_grp", "n_neg_grp",
    # variant proportions among mut+
    "pct_g12d_of_mutpos", "pct_g12v_of_mutpos", "pct_otherkras_of_mutpos",
    "pct_kras_of_mutpos", "pct_braf_of_mutpos",
    # sample-type breakdown
    "n_ffpe", "n_fresh", "n_literature", "n_institutional_tissue",
    # VAF overall
    "vaf_n", "vaf_median", "vaf_q25", "vaf_q75",
    "vaf_min", "vaf_max", "vaf_mean", "vaf_sd",
    # VAF per-variant
    "vaf_median_g12d", "vaf_q25_g12d", "vaf_q75_g12d",
    "vaf_median_g12v", "vaf_q25_g12v", "vaf_q75_g12v",
    "vaf_median_otherkras", "vaf_median_braf",
    "vaf_n_g12d", "vaf_n_g12v", "vaf_n_otherkras", "vaf_n_braf",
    # VAF tests
    "vaf_kw_stat", "vaf_kw_df", "vaf_kw_p",
    "vaf_pw_g12d_vs_g12v",
    "vaf_pw_g12d_vs_otherkras", "vaf_pw_g12v_vs_otherkras",
    # per-study statistics
    "chisq_study_genotype_stat", "chisq_study_genotype_p",
    "vaf_by_study_kw_stat", "vaf_by_study_kw_p",
    # full tables (data.frame objects)
    "mut_rate_by_study", "vaf_per_variant",
    "sample_type_breakdown", "variant_counts_full"
  ),

  # ---- Extended Data Fig. 01 (data completeness heatmap) ---------------------------
  # Producer: analysis/01_main_analysis/02_prep_analysis_dataset.R
  edfig01 = c(
    "n_vars_harmonised", "n_genotyped_in_heatmap", "n_studies_in_heatmap",
    "pct_missing_overall", "pct_available_overall",
    "n_complete_case_key_vars", "pct_complete_case_key_vars",
    "missingness_by_variable", "missingness_by_study",
    "missingness_matrix_long",
    "best_variable", "best_variable_pct_available",
    "worst_variable", "worst_variable_pct_available",
    "best_study", "best_study_pct_available",
    "worst_study", "worst_study_pct_available"
  ),

  # ---- Extended Data Fig. 02 (study-level heterogeneity: mut+ rate + VAF) ----------
  # Producer: analysis/01_main_analysis/08_F1_cohort_variants.R
  edfig02 = c(
    "n_studies", "n_tissue_tested", "n_mut_pos", "pooled_mut_rate_pct",
    # panel A
    "highest_rate_study", "highest_rate_pct",
    "highest_rate_n", "highest_rate_n_mut",
    "lowest_rate_study", "lowest_rate_pct",
    "lowest_rate_n", "lowest_rate_n_mut",
    "chisq_study_genotype_stat", "chisq_study_genotype_p",
    "mut_rate_by_study",
    # panel B
    "vaf_by_study_kw_stat", "vaf_by_study_kw_df", "vaf_by_study_kw_p",
    "vaf_by_study_summary"
  ),

  # ---- ED Fig 4 panels c–e (per-series REML meta-analyses) ---------------
  # Producer: analysis/01_main_analysis/11_ED4_power_forest_meta.R
  # Replaces the prior inverse-variance pooling: real metafor::rma REML
  # fits with Cochran Q, I^2, and tau^2 for each of the three outcomes.
  ed_cohort_heterogeneity = c(
    # Variant-positive rate meta (logit-link proportion)
    "rate_meta_k", "rate_meta_pooled_pct",
    "rate_meta_pooled_lo", "rate_meta_pooled_hi",
    "rate_meta_Q", "rate_meta_Q_df", "rate_meta_Q_p",
    "rate_meta_I2_pct", "rate_meta_tau2",
    # Rupture-prevalence meta (logit-link proportion)
    "rupt_prev_meta_k", "rupt_prev_meta_pooled_pct",
    "rupt_prev_meta_pooled_lo", "rupt_prev_meta_pooled_hi",
    "rupt_prev_meta_Q", "rupt_prev_meta_Q_df", "rupt_prev_meta_Q_p",
    "rupt_prev_meta_I2_pct", "rupt_prev_meta_tau2",
    # Age-at-presentation meta (sample-mean)
    "age_meta_k", "age_meta_pooled_y",
    "age_meta_pooled_lo", "age_meta_pooled_hi",
    "age_meta_Q", "age_meta_Q_df", "age_meta_Q_p",
    "age_meta_I2_pct", "age_meta_tau2",
    # Count of contributing series with both KRAS-mut and Negative arms
    # (used by §4 prose for "all N series with both arms")
    "age_meta_n_both_arms",
    # Pediatric-vs-adult ascertainment groups (data-driven, cut: per-series
    # mean age < 25 y) — supports the Results-paragraph framing
    "ped_n_series", "ped_series", "ped_age_min", "ped_age_max",
    "adult_n_series", "adult_series", "adult_age_min", "adult_age_max",
    # Bundled tables (carried for ED Fig 4 caption + SuppTable producer)
    "het_table", "series_detail", "ascertain_summary"
  ),

  # ---- Figure 2 (Section 2 — genotype vs. angioarchitecture null) --------
  # Producer: analysis/01_main_analysis/13_F2_genotype_phenotype.R
  fig2 = c(
    "cohort_n_sm_graded", "cohort_n_hr_forest",
    "sm_kw_binary_p", "sm_kw_variant_p",
    "sm_size_kw_p", "sm_drainage_fisher_p", "sm_eloquence_fisher_p",
    "sm_comp_min_p",
    # BH-FDR over the 3-component family (added 2026-05-19; prose cites
    # `sm_comp_min_q` to keep the "no SM sub-component differs by
    # genotype" claim defensible at the manuscript's FDR convention).
    "sm_size_kw_q", "sm_drainage_fisher_q", "sm_eloquence_fisher_q",
    "sm_comp_min_q",
    "hr_min_fdr", "hr_max_fdr",
    "rupt_meta_or", "rupt_meta_ci_lo", "rupt_meta_ci_hi",
    "rupt_meta_p", "rupt_meta_i2_pct", "rupt_meta_n_studies",
    "rupt_uni_or", "rupt_uni_ci_lo", "rupt_uni_ci_hi", "rupt_uni_p",
    "rupt_multi_or", "rupt_multi_ci_lo", "rupt_multi_ci_hi", "rupt_multi_p",
    "clinical_min_fdr", "clinical_max_fdr",
    "par_kras_rupt_pct", "par_nonkras_rupt_pct", "par_neg_rupt_pct",
    "par_kras_n",
    "par_rupt_interaction_beta", "par_rupt_interaction_p",
    "n_lobe_tests", "n_loc_int_outcomes", "n_loc_int_tests",
    "hr_forest", "location_distribution", "location_interactions"
  ),

  # ---- Figure 3 (Section 3 — age at presentation / KM / VAF × age) -------
  # Producer: analysis/01_main_analysis/09_F1_km_age.R
  # Scope: panels a (KM presentation), b (age density binary), c (age density
  # by variant), d (KM rupture), e (VAF × age scatter). The registry reorders
  # d/e vs the producer's old 4a/4b/4c/4d letters by first-citation order.
  fig3 = c(
    # Panel a — KM time to any clinical presentation
    "km_pres_n", "km_pres_n_g12d", "km_pres_n_g12v", "km_pres_n_neg",
    "km_pres_median_g12d", "km_pres_median_g12v", "km_pres_median_neg",
    "km_pres_ci_lo_g12d", "km_pres_ci_hi_g12d",
    "km_pres_ci_lo_g12v", "km_pres_ci_hi_g12v",
    "km_pres_ci_lo_neg",  "km_pres_ci_hi_neg",
    "km_pres_logrank_p",
    "km_pres_pw_g12d_neg", "km_pres_pw_g12v_neg", "km_pres_pw_g12d_g12v",
    # Panel b — age density, Mut+ vs Neg
    "age_bin_n", "age_bin_median_mut", "age_bin_median_neg", "age_bin_kw_p",
    # Pediatric (<18) Fisher from the same age-binary cohort, pulled in
    # here because the Results prose cites it alongside Panel b (see
    # `@fig[age_density_binary,age_density_variant]`).
    "ped_n_mut",   "ped_n_mut_u18",
    "ped_n_neg",   "ped_n_neg_u18",
    "ped_or", "ped_ci_lo", "ped_ci_hi", "ped_p",
    # Panel c — age density by variant (G12D / G12V / Neg)
    "age_var_n", "age_var_median_g12d", "age_var_median_g12v",
    "age_var_median_neg", "age_var_kw_p",
    # Panel e — VAF × age scatter
    "vaf_age_n", "vaf_age_n_g12d", "vaf_age_n_g12v",
    "vaf_age_slope_g12d",  "vaf_age_slope_g12d_p",
    "vaf_age_slope_g12v",  "vaf_age_slope_g12v_p",
    "vaf_age_interaction_p",
    "vaf_age_spearman_rho", "vaf_age_spearman_p",
    "vaf_age_spearman_rho_g12d", "vaf_age_spearman_p_g12d",
    "vaf_age_spearman_rho_g12v", "vaf_age_spearman_p_g12v",
    "vaf_age_g12d_cook_n_above4n",
    # Panel-d companion (moved from the producer's Fig4d → §4 rupture age):
    # median age at rupture per stratum for the prose in Results §4.
    "rupt_n", "rupt_events_total",
    "rupt_median_g12d", "rupt_median_g12v", "rupt_median_neg",
    "rupt_ci_lo_g12d",  "rupt_ci_hi_g12d",
    "rupt_ci_lo_g12v",  "rupt_ci_hi_g12v",
    "rupt_ci_lo_neg",   "rupt_ci_hi_neg",
    "rupt_logrank_p",
    # Unadjusted Cox — HEADLINE HR (Results §4 primary citation, pairs
    # one-for-one with the log-rank). 2026-05-19 reframe: see Methods
    # for why study is not in the primary multivariate model.
    "rupt_cox_unadj_hr_neg", "rupt_cox_unadj_ci_lo_neg",
    "rupt_cox_unadj_ci_hi_neg", "rupt_cox_unadj_p_neg",
    # Minimal-adjustment sensitivity (sex + sample_type; no study)
    "rupt_cox_minadj_hr_neg", "rupt_cox_minadj_ci_lo_neg",
    "rupt_cox_minadj_ci_hi_neg", "rupt_cox_minadj_p_neg",
    # Full-adjustment sensitivity (covariate-form study adjustment)
    "rupt_cox_hr_neg", "rupt_cox_ci_lo_neg", "rupt_cox_ci_hi_neg",
    "rupt_cox_p_neg",
    "rupt_cox_zph_global_p",
    # Angioarchitecture-adjusted Cox — cited in Results §4 main Cox sentence
    # (geno + nidus_size + deep_drainage + high_risk_count). Named _angio_
    # to distinguish from unadj / minadj / full / strata sensitivities.
    "rupt_cox_angio_n",
    "rupt_cox_angio_hr_neg", "rupt_cox_angio_ci_lo_neg",
    "rupt_cox_angio_ci_hi_neg", "rupt_cox_angio_p_neg",
    # G12D-vs-Negative direction (reciprocal of the *_neg keys); cited in §1
    "rupt_cox_angio_hr_g12d", "rupt_cox_angio_ci_lo_g12d",
    "rupt_cox_angio_ci_hi_g12d", "rupt_cox_angio_p_g12d",
    # Strata(study_clean) sensitivity Cox
    "rupt_cox_strata_hr_neg", "rupt_cox_strata_ci_lo_neg",
    "rupt_cox_strata_ci_hi_neg", "rupt_cox_strata_p_neg",
    "rupt_cox_strata_zph_global_p",
    # Full tables (for downstream consumers / ED replication)
    "km_pres_summary_table", "km_pres_pw_matrix",
    "rupt_summary_table",    "rupt_cox_coef"
  ),

  # ---- Figure 3 score panels (panels f + g) ------------------------------
  # Producer: analysis/01_main_analysis/21_F3_rupture_score_panels.R
  # Scope: Cox PH bedside score derived from the Option-2a strict cohort
  # used in panel d. Two main-text panels (km_by_score, rupture_lookup_
  # heatmap); the equation card and bootstrap validation live in the
  # ed_rupture_score_card section. Sized so a prose author can quote any
  # number cited in the integrated paragraph without running R.
  fig3_score = c(
    # cohort
    "n_eligible", "n_events", "n_train", "n_test",
    # selected anchor features (T1 free RFE, refit on all eligible)
    "anchor_features",
    # unrounded Cox HRs
    "hr_g12d", "hr_g12d_lo", "hr_g12d_hi", "hr_g12d_p",
    "hr_drainage", "hr_drainage_lo", "hr_drainage_hi", "hr_drainage_p",
    "hr_size", "hr_size_lo", "hr_size_hi", "hr_size_p",
    "zph_global_p", "zph_g12d_p", "zph_drainage_p", "zph_size_p",
    # integer card definition
    "scale_factor", "pts_g12d", "pts_drainage", "pts_size",
    "score_min", "score_max",
    # integer-score Cox refit
    "score1_hr", "score1_hr_lo", "score1_hr_hi", "score1_p",
    "score2_hr", "score2_hr_lo", "score2_hr_hi", "score2_p",
    # holdout discrimination
    "c_holdout_lp", "c_holdout_int",
    "c_holdout_delta",
    # bootstrap validation (Harrell, B = 500)
    "boot_B", "c_apparent",
    "c_corrected", "optimism_mean",
    "optimism_ci_lo", "optimism_ci_hi",
    # lookup grid (point estimate + bootstrap CIs at named ages)
    # rendered as a wide table for the heatmap panel + caption
    "lookup_grid_point", "lookup_grid_ci",
    # AJK #C38: youngest REF_AGES age (years) at which the score 0 vs 2
    # bootstrap 95% CIs become non-overlapping
    "score_separation_age_y",
    # convenience scalars for the prose paragraph (score 0 vs 2 by age 20)
    "risk_age20_score0_pct", "risk_age20_score0_lo", "risk_age20_score0_hi",
    "risk_age20_score2_pct", "risk_age20_score2_lo", "risk_age20_score2_hi",
    "risk_age10_score0_pct", "risk_age10_score2_pct",
    "risk_age30_score0_pct", "risk_age30_score2_pct",
    # patient counts per integer score
    "n_score0", "n_score1", "n_score2",
    "events_score0", "events_score1", "events_score2"
  ),

  # ed_rupture_score_card section retired 2026-04-26: producer 30 (the
  # ed_rupture_card script) was kept for the SuppTable_score_card.xlsx
  # output but no longer emits a manifest fragment — every value the
  # §4 prose cites is duplicated in `fig3_score` written by
  # 21_F3_rupture_score_panels.R, so there's no consumer for a separate
  # ed_rupture_score_card section.

  # ---- Extended Data Fig. 16 (time to rupture KM) ------------------------
  # Producer: analysis/01_main_analysis/09_F1_km_age.R
  # Keys largely overlap with fig3 rupture-age block; edfig16 is the
  # supporting view so we duplicate the scalars rather than reach into
  # `fig3` from the ED caption.
  edfig16 = c(
    "n", "n_events", "n_censored",
    "median_g12d", "median_g12v", "median_neg",
    "logrank_p",
    "pw_g12d_neg", "pw_g12v_neg", "pw_g12d_g12v",
    "cox_hr_neg", "cox_ci_lo_neg", "cox_ci_hi_neg", "cox_p_neg",
    "summary_table", "cox_coef"
  ),

  # ---- Extended Data Fig. 17 (rare-variant KM) ---------------------------
  # Producer: analysis/01_main_analysis/09_F1_km_age.R
  edfig17 = c(
    "n", "n_other_kras", "n_braf_v600e", "n_neg",
    "median_other_kras", "median_braf_v600e", "median_neg",
    "logrank_p",
    "pw_other_kras_neg", "pw_braf_v600e_neg", "pw_other_kras_braf",
    "summary_table"
  ),

  # ---- Figure 4 (venous stenosis waffle; scRNA panels retired 2026-05-19) --
  # Producer: analysis/01_main_analysis/19_ED10_venous_stenosis.R
  # Scope: stenosis waffle + Firth-OR ancillary numbers. scRNA panels and
  # their stats sections were retired with the rest of the scRNA story
  # (now Fig 1 of a forthcoming scRNA-only paper). Many sc_* keys remain
  # listed below for backwards-compat with older fragments; the producer
  # no longer writes them. Phase 3 of the restructure merges this group
  # into gxp_associations / new Fig 2 alongside null_phenotype panels.
  fig4 = c(
    # Panel a + b — stenosis 2x2 + forest
    "stenosis_n_total", "stenosis_n_mutpos", "stenosis_n_neg",
    "stenosis_n_mutpos_stenosis", "stenosis_n_neg_stenosis",
    "or_simple", "or_simple_lo", "or_simple_hi",
    "or_firth",  "or_firth_lo",  "or_firth_hi",  "or_firth_p",
    "firth_n_adjusted",
    # Per-study ascertainment (cited in §5 caveat)
    "stenosis_n_bch", "stenosis_n_chop", "stenosis_n_uab",
    "stenosis_total_positive"
    # sc_* keys removed 2026-05-23: single-cell analysis dropped categorically.
  ),

  # Fig 4d endovascular liquid-biopsy VAF (TERT-referenced), produced by
  # 25_F4_biopsy_vaf.R from this repo's raw QIAcuity export + cohort tissue.
  biopsy = c(
    "cfdna_patient", "cfdna_g12d_partitions", "cfdna_tert_partitions",
    "cfdna_vaf_pct", "tissue_sample_uid", "tissue_vaf_pct"
  ),

  # ---- Extended Data Fig 22 (sex-stratified KM) -------------------------
  # Producer: analysis/01_main_analysis/12_ED6_km_sex_stratified.R
  # §4 prose cites "genotype × sex interaction P = 0.712" for presentation.
  edfig22 = c(
    "n_pres", "n_pres_female", "n_pres_male",
    "logrank_p_pres_female", "logrank_p_pres_male",
    "cox_lrt_chisq_pres", "cox_lrt_df_pres", "cox_lrt_p_pres",
    "n_rupt", "n_rupt_events", "n_rupt_female", "n_rupt_male",
    "events_rupt_female", "events_rupt_male",
    "logrank_p_rupt_female", "logrank_p_rupt_male",
    "cox_lrt_chisq_rupt", "cox_lrt_df_rupt", "cox_lrt_p_rupt",
    "cox_main_coef_pres", "cox_int_coef_pres",
    "cox_main_coef_rupt", "cox_int_coef_rupt"
  ),

  # ---- Extended Data Fig 08 (anatomy) -----------------------------------
  # Producer: analysis/01_main_analysis/17_ED8_anatomy.R
  # §3 prose cites "all BH-adjusted P >= 0.83" across the 9 lobes plus
  # "27 location-by-genotype interaction tests" (parietal × rupture = sole
  # nominal-sig cell).
  edfig08_anatomy = c(
    "n_lobes",
    "panelA_min_fisher_p", "panelA_min_p_fdr", "panelA_max_p_fdr",
    "n_interaction_tests", "n_outcome_families", "n_tests_nominal_sig",
    "panelC_min_p_fdr", "panelC_min_p_bonf",
    "per_lobe_prevalence", "per_variant_anatomy", "interaction_tests"
  ),

  # ---- Extended Data Fig 4 (VAF × phenotype correlations) ---------------
  # Producer: analysis/01_main_analysis/14_ED7_vaf_phenotype.R (renamed from
  # the legacy 04_fig3.R on 2026-04-25; script feeds the resolver-assigned
  # ED4 group `ed_vaf_phenotype`). Key-space covers the full Spearman-rho +
  # Wilcoxon panel cited in Results §2: "VAF did not scale with SM grade,
  # rupture prevalence, deep venous drainage, eloquence, or high-risk
  # feature count (all Spearman ρ ∈ [−0.07, 0.04]; all P ≥ 0.13)".
  edfig_vaf_phenotype = c(
    # Primary correlations (Spearman for ordinal; Wilcoxon for binary)
    "rho_sm_total",       "p_sm_total",       "n_sm_total",
    "rho_sm_size",        "p_sm_size",        "n_sm_size",
    "rho_high_risk",      "p_high_risk",      "n_high_risk",
    "rho_age",            "p_age",            "n_age",
    "wilcox_drainage",    "p_drainage",       "n_drainage",
    "wilcox_eloquence",   "p_eloquence",      "n_eloquence",
    "wilcox_rupture",     "p_rupture",        "n_rupture",
    # Headline aggregates for §2 prose
    "rho_min", "rho_max", "p_min", "p_max",
    "q_fdr_min", "q_fdr_max",
    # Rupture ~ VAF adjusted (multivariable) — headline row
    "rupt_vaf_or", "rupt_vaf_ci_lo", "rupt_vaf_ci_hi", "rupt_vaf_p",
    # Full per-outcome table
    "correlation_table",
    "rupture_adjusted_coefs"
  ),

  # ---- Extended Data Fig. 14 (stenosis ascertainment bar) ---------------
  # Produced alongside Fig 4 by 19_ED10_venous_stenosis.R. Kept as its own
  # section so the ED caption reads from a focused key space.
  edfig14_stenosis_asc = c(
    "n_total", "n_positive",
    "n_bch", "n_bch_positive",
    "n_chop", "n_chop_positive",
    "n_uab", "n_uab_positive",
    "study_table"
  ),

  # ---- ED Fig: VAF outlier sensitivity / dose-response diagnostics ------
  # Producer: analysis/01_main_analysis/16_ST_vaf_outlier_sensitivity.R. Surfaces the numeric
  # handles the caption needs (Shapiro-Wilk W + p, n, % outliers) so caption
  # text replaces the "annotated on-panel" claim with actual values pulled
  # from this fragment.
  edfig_vaf_outlier = c(
    "n_vaf", "vaf_mean_pct", "vaf_sd_pct",
    "shapiro_W_raw", "shapiro_p_raw",
    "shapiro_W_log", "shapiro_p_log",
    "n_outliers", "pct_outliers",
    "n_full_cohort", "n_trim_cohort"
  )

  # edfig19 (VAF x age slope sensitivity) retired 2026-04-26: the figure was
  # removed and the sensitivity content now lives only in Supplementary Table
  # `vaf_age_sensitivity` (sheets 1-3 written by 20_ST_supp_tables_8_9_10.R, plus
  # auxiliary sheets appended by 16_ST_vaf_outlier_sensitivity.R). No stats keys are
  # consumed by the manuscript prose (they were caption-only).

  # Minimum-detectable-effect (MDE) summary derived from the per-outcome MDE
  # table (results/stats/mde_table.rds). Used by §4 prose to avoid hardcoded
  # "ruling out odds ratios >= 2-4" thresholds.
  ,
  mde = c("mde_or_binary_min", "mde_or_binary_max",
          "mde_d_continuous_round1")
)

# Which script is expected to produce each section. Used by
# check_stats_manifest.R to point the user at the right re-run target if a
# fragment is stale.
stats_producers <- list(
  cohort  = "analysis/01_main_analysis/03_prep_cohort_counts.R",
  fig1    = "analysis/01_main_analysis/08_F1_cohort_variants.R",
  fig2    = "analysis/01_main_analysis/13_F2_genotype_phenotype.R",
  fig3    = "analysis/01_main_analysis/09_F1_km_age.R",
  fig4    = "analysis/01_main_analysis/19_ED10_venous_stenosis.R",
  biopsy  = "analysis/01_main_analysis/25_F4_biopsy_vaf.R",
  edfig01 = "analysis/01_main_analysis/02_prep_analysis_dataset.R",
  edfig02 = "analysis/01_main_analysis/08_F1_cohort_variants.R",
  ed_cohort_heterogeneity = "analysis/01_main_analysis/11_ED4_power_forest_meta.R",
  mde                     = "analysis/01_main_analysis/11_ED4_power_forest_meta.R",
  edfig_vaf_phenotype  = "analysis/01_main_analysis/14_ED7_vaf_phenotype.R",
  edfig14_stenosis_asc = "analysis/01_main_analysis/19_ED10_venous_stenosis.R",
  edfig08_anatomy                 = "analysis/01_main_analysis/17_ED8_anatomy.R",
  edfig22                         = "analysis/01_main_analysis/12_ED6_km_sex_stratified.R",
  edfig16 = "analysis/01_main_analysis/09_F1_km_age.R",
  edfig17 = "analysis/01_main_analysis/09_F1_km_age.R",
  edfig_vaf_outlier      = "analysis/01_main_analysis/16_ST_vaf_outlier_sensitivity.R",
  fig3_score             = "analysis/01_main_analysis/21_F3_rupture_score_panels.R",
  vaf_rupture            = "analysis/01_main_analysis/15_ED7_vaf_age_at_rupture.R",
  demographics           = "analysis/01_main_analysis/06_ST_residual_demographics.R"
)
