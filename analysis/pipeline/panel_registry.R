# =============================================================================
# panel_registry.R — source of truth for figure/table tokens
# -----------------------------------------------------------------------------
# Hand-edited catalog that the prose resolver (analysis/pipeline/helpers/resolve_panels.R)
# consumes. Every panel and table referenced from the manuscript narrative must
# have an entry here.
#
# Model:
#   figures$<group>$track   — "Fig" or "EDFig"
#   figures$<group>$panels  — ordered character vector of panel tokens;
#                             order here is only used as a default composite
#                             layout hint. Panel LETTERS are assigned by
#                             order of first appearance in prose, so moving
#                             text around reletters the composite.
#   tables$<token>          — track string ("Table" or "SuppTable"); the token
#                             is the citation handle. Table NUMBERS are
#                             assigned by order of first appearance in prose
#                             per track.
#
# Citation syntax in prose (scanned by resolve_panels.R):
#   @fig[token]               -> **Fig. Nx**
#   @fig[tok1,tok2]           -> **Fig. Na,b**   (2 panels, comma)
#   @fig[tok1,tok2,tok3]      -> **Fig. Na-c**   (>=3 contiguous, hyphen)
#   @fig[tok1,tok3]           -> **Fig. Na,c**   (non-contiguous, comma)
#   @edfig[...]               -> **Extended Data Fig. ...**
#   @tab[token]               -> **Table N**
#   @supptab[tok1,tok2,tok3]  -> **Supplementary Tables N-M**  (range rules
#                                identical to panels)
#
# Tokens are globally unique across the figures namespace so the resolver can
# identify which group a panel belongs to from the token alone. Tables live
# in a flat token namespace.
# =============================================================================

panel_registry <- list(

  # ══ FIGURES (main + Extended Data) ══════════════════════════════════════
  # Group keys are author-facing slugs. The manuscript-facing figure NUMBER
  # is derived from order of first appearance in prose, per track.

  figures = list(

    # ---- Main Figures (post 2026-05-19 restructure: 6 -> 3 main figs) --
    # The original 6-group main-figure layout (cohort_variants, age,
    # null_phenotype, venous_stenosis, liquid_biopsy, precision_framework)
    # collapsed into 3 merged groups after the scRNA retirement. Each
    # merged group has author-facing slug naming the unified narrative.
    # First-citation order in the prose determines Fig 1/2/3 assignment.

    # Fig 1: cohort + genotype-stratified natural history
    # Combines old cohort_variants (CONSORT + variant landscape + VAF
    # distribution) with old age (KM curves, age densities, VAF x age
    # dose response, integer rupture-risk score). 10 panels total;
    # layout decisions on which to drop are deferred to the design pass.
    cohort_natural_history = list(
      track  = "Fig",
      # 2026-05-19 layout decision: 9 panels in 3x3 grid; dropped
      # age_density_binary (visually redundant with age_density_variant
      # which already includes panel-negative). Producer 09_F1_km_age.R
      # still saves the binary panel to disk; it just isn't registered
      # or composed into the main figure.
      # 2026-05-23: km_by_score + rupture_lookup_heatmap moved to
      # standalone rupture_score group (new Fig 3). 7 panels remain.
      panels = c("consort_flow", "variant_landscape", "vaf_by_variant",
                 "km_presentation", "age_density_variant",
                 "km_rupture", "vaf_age_scatter")
    ),

    # Fig 2: genotype-phenotype associations + venous outflow stenosis
    # Combines old null_phenotype (parietal x KRAS x rupture, SM grade,
    # SM components, HR features, clinical history) with the surviving
    # stenosis_waffle from old venous_stenosis. 6 panels.
    gxp_associations = list(
      track  = "Fig",
      # 2026-05-20 v3: stenosis_waffle dropped from Fig 2; the
      # venous-outflow-stenosis finding still lives as one of the four
      # features in the hr_features_OR forest, so the waffle was
      # redundant. Five panels (letters = first-citation order in prose;
      # visual reading order a→e is sequential left-to-right top-to-bottom):
      #   a sm_grade_dist               stacked bar (top-left, 3 cols)
      #   b sm_components_dumbbell      SM size / drainage / eloquence (top-right, 2 cols)
      #   c clinical_history_dumbbell   seizure / radiation / embolization (mid-left)
      #   d parietal_kras_rupture       lollipop (mid-right)
      #   e hr_features_OR              high-risk feature forest (full-width bottom)
      panels = c("sm_grade_dist", "sm_components_dumbbell",
                 "clinical_history_dumbbell", "parietal_kras_rupture",
                 "hr_features_OR")
    ),

    # Fig 3: genotype-informed rupture-risk score
    # Two panels moved from cohort_natural_history 2026-05-23 to give the
    # score section its own standalone main figure between gxp_associations
    # and precision_medicine.
    rupture_score = list(
      track  = "Fig",
      panels = c("km_by_score", "rupture_lookup_heatmap")
    ),

    # Fig 4: precision medicine (liquid biopsy + framework)
    # Combines old liquid_biopsy (biopsy schematic, angiograms,
    # TapeStation, dPCR waterfall) with the framework diagram. 5 panels.
    precision_medicine = list(
      track  = "Fig",
      panels = c("biopsy_schematic", "biopsy_angiograms",
                 "biopsy_tapestation", "biopsy_dpcr_waterfall",
                 "framework_diagram")
    ),

    # ---- Extended Data Figures -----------------------------------------
    ed_prisma = list(
      track  = "EDFig",
      panels = "prisma_flow"
    ),

    ed_completeness = list(
      track  = "EDFig",
      panels = "missingness_heatmap"
    ),

    # Multiplex dPCR assay validation. Six per-variant scatter plots
    # demonstrating fluorescence-amplitude clustering for each probe channel,
    # including the TERT internal positive control. Each panel: positive
    # control (left) vs non-template control / NTC (right).
    # Producer: analysis/01_main_analysis/23_ED3_dpcr_validation.R reads the
    # de-identified REAL QIAcuity extract at
    # ../avm-variant-detection/data/dpcr_real/edfig3_dpcr_real.csv, built from
    # the "UAB Priority Samples Run 3" exports by
    # ../avm-variant-detection/scripts/09_figure_generation/edfig3_dpcr_validation/extract_real_dpcr.py
    # Channel->variant: G=G12D, Y=G12V, O=BRAF V600E, R=G12C, C=G12A, Fr=TERT.
    ed_dpcr_validation = list(
      track  = "EDFig",
      panels = c("dpcr_scatter_g12d", "dpcr_scatter_g12v",
                 "dpcr_scatter_g12c", "dpcr_scatter_g12a",
                 "dpcr_scatter_v600e", "dpcr_scatter_tert")
    ),

    # 2026-05-20 (Phase 3): ed_cohort_heterogeneity split into TWO ED
    # figures because 6 panels at Nature double-col (7.20 in) made every
    # hand-rolled forest text element collide. New layout — both at
    # 3 panels in a 3-column or 1-col-x-3-row grid:
    #   ed_cohort_heterogeneity (this group) — per-series spread + pooled
    #     mutation rate (the "cohorts are heterogeneous" story).
    #   ed_age_meta_variants    (new group)  — phenotype-pooled meta
    #     forests (rupture, age) + rare-variant KM (the "but the genotype
    #     signal converges across cohorts" story). Cited from the §1
    #     rupture / §2 age paragraphs.
    ed_cohort_heterogeneity = list(
      track  = "EDFig",
      panels = c("per_series_rate", "per_series_vaf", "pooled_rate_meta")
    ),

    # v6.22 (2026-05-20): rare_variants moved out of ed_age_meta_variants
    # into ed_km_diagnostics (the natural home for any KM panel). This
    # group now carries only the two phenotype meta-forests.
    ed_age_meta_variants = list(
      track  = "EDFig",
      panels = c("rupture_meta", "age_meta")
    ),

    # 2026-05-19 (Phase 2 / Iteration 3): ed_vaf_phenotype +
    # ed_vaf_outlier + ed_vaf_age_rupture merged into ed_vaf_deep_dive,
    # bringing the ED count to the Plan-B target of 9.
    # v6.29 (2026-05-20): vaf_outlier_combined dropped; 7 panels.
    # 2026-06-04: panel order corrected so the VAF->age-at-rupture
    # dose-response scatter is FIRST (panel a), then the six VAF x
    # phenotype null comparisons (b-g):
    #   a    vaf_age_rupture_scatter   (VAF -> age-at-rupture dose-response)
    #   b-g  vaf_sm_total / vaf_sm_size / vaf_drainage / vaf_eloquence
    #        / vaf_rupture / vaf_highrisk        (vaf x phenotype nulls)
    # Rationale: both prose citations are group-level (@edfig[ed_vaf_deep_dive]),
    # so the resolver assigns letters in this registry order. The figure is
    # first invoked in §2 for the dose-response (alongside the rupture claim),
    # then re-cited in §3 for the null comparisons — so first-appearance order
    # puts the scatter at panel a. Composer 49 already lays the scatter out as
    # panel a, so figure and caption now agree. (Supersedes the earlier
    # "v6.53 moved scatter to last" decision, which misread the prose order.)
    ed_vaf_deep_dive = list(
      track  = "EDFig",
      panels = c("vaf_age_rupture_scatter",
                 "vaf_sm_total", "vaf_sm_size", "vaf_drainage",
                 "vaf_eloquence", "vaf_rupture", "vaf_highrisk")
    ),

    ed_per_variant_pheno = list(
      track  = "EDFig",
      panels = c("sm_ordinal", "sm_comp_variant",
                 "hr_variant", "clinical_variant")
    ),

    # 2026-04-26 merge: ed_mde_audit, ed_age_adj_vaf_pheno, and
    # ed_stenosis_ascertainment all attach to the same null-phenotype
    # paragraph in §3 (MDE coverage of the null comparisons, age-adjusted
    # VAF×phenotype sensitivity, and per-series stenosis ascertainment
    # caveat). Combined into a single 3-panel ED figure ed_null_audit so
    # the reader sees all three caveats in one place.
    # v6.45 (2026-05-21): age_adj_correlation split into two panels
    # (binary vs continuous outcomes) so each gets its own 3.30-in cell
    # in the ED10 2x2 grid. Total panels: 4.
    ed_null_audit = list(
      track  = "EDFig",
      panels = c("mde_all", "stenosis_ascertainment",
                 "age_adj_binary", "age_adj_continuous")
    ),

    # 2026-05-19 (Phase 2 / Iteration 1): ed_anatomy + ed_parietal merged
    # into ed_anatomic_localization to bring the ED count under Nature's
    # cap (13 -> 9 across three merges). Panels render in 2x3 layout with
    # the three anatomy panels (a/b/c) on top and the three parietal
    # deep-dive panels (d/e/f) below.
    #
    #   a — per_lobe_prevalence    9-lobe prevalence heatmap (mut+ vs neg)
    #   b — per_variant_anatomy    per-variant lobe composition
    #   c — interaction_cleveland  27-cell location × outcome cleveland
    #   d — rupture_categories     3-group rupture-category mosaic
    #                              (Negative / KRAS-NonParietal / KRAS-Parietal)
    #   e — parietal_forest        adjusted-OR forest, parietal subgroup
    #   f — rupture_meta_forest    within-study log-OR meta-analysis
    # v6.51 (2026-05-21): rupture_meta_forest dropped — rupture prevalence
    # is no longer claimed in the §3 narrative (the within-study log-OR
    # meta-analysis is implicit in the parietal-by-rupture finding and
    # the cohort_heterogeneity composite). Four panels remain.
    ed_anatomic_localization = list(
      track  = "EDFig",
      panels = c("per_lobe_prevalence", "per_variant_anatomy",
                 "rupture_categories", "parietal_forest")
    ),

    # ed_rupture_score_card group retired 2026-04-26: the score-derivation
    # numbers (anchor-model hazard ratios + Harrell bootstrap validation)
    # are tiny tabular results — three HR rows and three C-index numbers —
    # that read more cleanly as a Supplementary Table than as a
    # forest + lollipop figure. They now live in
    # SuppTable_score_card.xlsx (registered below as `score_card`).
    # The clinician-facing equation card (score_equation) was retired the
    # same day; the equation is fully described in §4 prose.

    # 2026-05-19 (Phase 2 / Iteration 2): ed_sex_stratified_km + ed_km_atrisk
    # merged into ed_km_diagnostics. Seven panels:
    #   a/b  km_pres_sex_F / km_pres_sex_M  (sex-stratified presentation KM)
    #   c/d  km_rupt_sex_F / km_rupt_sex_M  (sex-stratified rupture KM)
    #   e    km_pres_atrisk                 (presentation curve + at-risk table)
    #   f    km_rupt_atrisk                 (rupture curve + at-risk table)
    #   g    km_score_atrisk                (by-score curve + at-risk table)
    # The at-risk panels reuse the Fig 2 km_panel pointer RDS files; the
    # globally-unique-token constraint holds because the at-risk tokens
    # carry the `_atrisk` suffix distinct from the Fig 2 producers.
    # v6.22 (2026-05-20): the 3 *_atrisk compound panels (curve + at-risk
    # table) dropped — the tables ate the visual budget without adding
    # quantitative info that isn't already in the rare_variants n=*
    # annotations. rare_variants moved in from the (split) ed_age_meta_
    # variants group so all KM panels live in one place.
    ed_km_diagnostics = list(
      track  = "EDFig",
      panels = c("rare_variants",
                 "km_pres_sex_F", "km_pres_sex_M",
                 "km_rupt_sex_F", "km_rupt_sex_M")
    )

    # ed_vaf_age_sensitivity figure retired 2026-04-26: the LOO / outlier-trim /
    # institutional-restriction analyses were folded into Supplementary Table
    # `vaf_age_sensitivity` (sheets 1-3 in SuppTable10_vaf_age_sensitivity.xlsx,
    # populated by 20_ST_supp_tables_8_9_10.R). Keeping the supp-table token
    # (registered below) means the prose @supptab[vaf_age_sensitivity]
    # citation continues to work.

    # ed_vaf_outlier + ed_vaf_age_rupture entries retired in the
    # 2026-05-19 Phase-2 merge into ed_vaf_deep_dive (above). Producers
    # 16_ST_vaf_outlier_sensitivity.R and 15_ED7_vaf_age_at_rupture.R now write directly
    # into ed_vaf_deep_dive/ alongside the six vaf_phenotype panels.

    # ---- scRNA-seq Extended Data (retired 2026-05-19) ----------------
    # Six ED groups (ed_scrna_methods, ed_scrna_phred_qc,
    # ed_scrna_per_patient_sorted, ed_scrna_endothelial_restriction,
    # ed_scrna_subtype_vaf, ed_scrna_subtype_sensitivity) were retired
    # along with the rest of the scRNA story. The avm-variant-detection
    # sibling repo continues to produce these panels for a forthcoming
    # scRNA-only paper; remove from this manuscript's panel registry
    # only. Producers moved to analysis/_retired_/scrna_2026-05-19/.

    # ---- Sample-identity verification (registered, producer pending) ----
    # Verifies sc / bulk pairing claim that ED-scrna_phred_qc's scbulk
    # per-patient panel rests on. Producer at
    #   ../avm-variant-detection/scripts/01_alignment/m_somalier_extract_relate.lsf
    # has already run on the 7 scRNA BAMs (5 sorted_EC + 2 unsorted) and
    # produced an interim sc-only matrix at
    #   ../avm-variant-detection/results/01_sample_identity/somalier/somalier_heatmap.{png,pdf,svg}
    # confirming AVM4 + AVM5 paired-prep identity (relatedness 0.96 / 0.98)
    # and 19 unrelated cross-patient pairs (relatedness ≤ 0.005).
    #
    # Bulk RNA-seq BAMs are 0-byte placeholders (purged from scratch + project
    # clones); FASTQ source not yet identified — see MINERVA_TODO M-B11. Once
    # bulk is re-aligned, the somalier matrix will be re-rendered to include
    # all 14 BAMs (7 bulk + 5 sorted_EC + 2 unsorted) and this panel will
    # cite the full sc/bulk identity proof. Until then the token is registered
    # so manuscript prose can already reference it via @edfig[scrna_sample_identity].
    #
    # Pending producer: analysis/01_main_analysis/24_ed_scrna_sample_identity.R
    # (TODO — wire after bulk re-alignment lands).
    #
    # 2026-05-12: TEMPORARILY COMMENTED OUT so the strict registry-vs-disk
    # validator passes. No prose currently references this token (verified
    # with grep across sections/, captions.R, and ExtendedData.md — only
    # the caption stub at captions.R::ed_scrna_sample_identity references
    # it, and that stub is unused without a registered token). Restore
    # this block as a one-line uncomment when the producer lands.
    # ed_scrna_sample_identity = list(
    #   track  = "EDFig",
    #   panels = "scrna_sample_identity_heatmap"
    # )
  ),

  # ══ TABLES (main + Supplementary) ═══════════════════════════════════════
  # Flat token -> track namespace. Table NUMBERS are derived from order of
  # first appearance in prose, per track.

  tables = list(
    # Main tables
    demographics           = "Table",

    # Supplementary tables
    mutation_rate_by_study = "SuppTable",
    variant_freq           = "SuppTable",
    cohort_heterogeneity   = "SuppTable",
    hierarchical_gxp       = "SuppTable",
    high_risk_OR           = "SuppTable",
    vaf_phenotype_corr     = "SuppTable",
    mde_audit              = "SuppTable",
    parietal_kras          = "SuppTable",
    location_interactions  = "SuppTable",
    coarse_groupings       = "SuppTable",
    age_phenotype_interact = "SuppTable",
    vaf_age_sensitivity    = "SuppTable",
    vaf_age_rupture        = "SuppTable",
    # Bedside rupture-age score (Cox anchor-model HRs + Harrell bootstrap
    # validation). Replaced the retired ed_rupture_score_card ED figure
    # 2026-04-26; numbers now live as a single supplementary table.
    score_card             = "SuppTable"

    # scRNA variant-detection supplements: ALL retired 2026-05-19 along
    # with the rest of the scRNA story. The avm-variant-detection sibling
    # repo still owns these tables for a forthcoming scRNA-only paper.
  )
)
