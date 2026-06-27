# Script index — `analysis/01_main_analysis/`

Scripts are named `NN_<TARGET>_<desc>.R` (`F1`–`F4`, `ED1`–`ED10`, `T1` = Table 1,
`ST` = Supplementary Tables, `prep` = setup) and grouped so each **main figure sits
with its companion Extended Data figures**. `run_all.R` runs them in number order.

**Flow:** data prep → producers (`02`–`25`) → composers (`26`–`37`). Producers are
independent (each reads only the prepared dataset); composers run after all producers.
Several producers are multi-figure (e.g. `08`→Fig 1+ED 4, `13`→Fig 2+ED 9,
`11`→ED 4/5/10) — the **Produces** column is the authoritative map.

## Data prep — `analysis/00_data_prep/` (run first)

| # | script | produces |
|---|--------|----------|
| 00 | `00_mask_phi.R` | masks identifiers in the raw master → `phi_safe/` |
| 01 | `01_clean_master.R` | column-by-column cleaning → analysis-ready dataset |

## Producers — `02`–`25`

### Setup & cohort description
| # | script | produces |
|---|--------|----------|
| 02 | `02_prep_analysis_dataset.R` | analysis-ready dataset + cohort subsetting; **ED 2** (completeness) |
| 03 | `03_prep_cohort_counts.R` | main-text cohort counts |
| 04 | `04_T1_table1.R` | **Table 1** |
| 05 | `05_T1_table1_detailed.R` | **Table 1** (detailed, per-variant) |
| 06 | `06_ST_residual_demographics.R` | **Table 1** residuals + Supp Tables |

### Figure 1 + companions (ED 1, 4, 5, 6)
| # | script | produces |
|---|--------|----------|
| 07 | `07_F1_consort_flow.R` | **Fig 1A** (CONSORT flow) |
| 08 | `08_F1_cohort_variants.R` | **Fig 1B,C** + **ED 4A,B** + Supp Tables |
| 09 | `09_F1_km_age.R` | **Fig 1D–G** + **ED 6A** |
| 10 | `10_ED1_prisma_flow.R` | **ED 1** (PRISMA flow) |
| 11 | `11_ED4_power_forest_meta.R` | **ED 4C** + **ED 5A,B** + **ED 10A** + Supp Tables |
| 12 | `12_ED6_km_sex_stratified.R` | **ED 6B–E** (sex-stratified KM) |

### Figure 2 + companions (ED 7, 8, 9, 10)
| # | script | produces |
|---|--------|----------|
| 13 | `13_F2_genotype_phenotype.R` | **Fig 2A–E** + **ED 9A–D** + Supp Tables |
| 14 | `14_ED7_vaf_phenotype.R` | **ED 7B–G** + Supp Tables |
| 15 | `15_ED7_vaf_age_at_rupture.R` | **ED 7A** + Supp Tables |
| 16 | `16_ST_vaf_outlier_sensitivity.R` | Supp Table + VAF-outlier sensitivity stats |
| 17 | `17_ED8_anatomy.R` | **ED 8A,B** + Supp Tables |
| 18 | `18_ED10_age_adj_vaf_pheno.R` | **ED 10C,D** + Supp Tables |
| 19 | `19_ED10_venous_stenosis.R` | **ED 10B** + Fig 4 venous-stenosis stats |
| 20 | `20_ST_supp_tables_8_9_10.R` | **Supp Tables 8–10** |

### Figure 3
| # | script | produces |
|---|--------|----------|
| 21 | `21_F3_rupture_score_panels.R` | **Fig 3A,B** (bedside rupture-risk score) |
| 22 | `22_ST_rupture_score_card.R` | **Supp Table** (rupture-age score card) |

### Figure 4 + companion (ED 3)
| # | script | produces |
|---|--------|----------|
| 23 | `23_ED3_dpcr_validation.R` | **ED 3A–F** (multiplex dPCR validation) |
| 24 | `24_F4_dpcr_waterfall.R` | **Fig 4D** (endovascular liquid-biopsy dPCR) |
| 25 | `25_F4_biopsy_vaf.R` | **Fig 4** biopsy VAF (stats feeding panel D) |

## Composers — `26`–`37` (one per figure; assemble producer panels)

| # | script | assembles |
|---|--------|-----------|
| 26 | `26_F1_assemble.R` | **Figure 1** |
| 27 | `27_ED4_assemble.R` | **Extended Data Fig 4** |
| 28 | `28_ED5_assemble.R` | **Extended Data Fig 5** |
| 29 | `29_ED6_assemble.R` | **Extended Data Fig 6** |
| 30 | `30_F2_assemble.R` | **Figure 2** |
| 31 | `31_ED7_assemble.R` | **Extended Data Fig 7** |
| 32 | `32_ED8_assemble.R` | **Extended Data Fig 8** |
| 33 | `33_ED9_assemble.R` | **Extended Data Fig 9** |
| 34 | `34_ED10_assemble.R` | **Extended Data Fig 10** |
| 35 | `35_F3_assemble.R` | **Figure 3** |
| 36 | `36_F4_assemble.R` | **Figure 4** |
| 37 | `37_ED3_assemble.R` | **Extended Data Fig 3** |

Figure numbers and panel letters come from the frozen
`results/stats/panel_assignments.rds`. Extended Data Fig 2 (completeness) has no
dedicated script — it is emitted by `02_prep_analysis_dataset.R`.
