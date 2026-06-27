<div align="center">
  <img src="docs/kundishora-lab-logo.jpeg" alt="Kundishora Lab" width="480"/>
</div>

# Precision diagnosis and genotype-stratified natural history of brain arteriovenous malformations

Analysis code and figure/table outputs for the study, from the **Kundishora Lab**
(Hale & Keller et al.).

Brain arteriovenous malformations (bAVMs) are high-flow cerebrovascular lesions
associated with a lifelong risk of intracranial hemorrhage, yet no biologically
anchored framework exists to predict their clinical trajectory or guide medical
therapy. Although activating somatic variants in *KRAS* and *BRAF* drive sporadic
bAVMs, the relationship between genotype, disease severity, and clinical behavior has
remained poorly defined, and molecular diagnosis has largely been restricted to
surgically resected tissue.

This study assembles the largest integrated genomic–phenomic cohort of sporadic bAVMs
to date (**n = 475**) — harmonizing somatic genotypes, clinical phenotypes,
angioarchitectural features, and neuroanatomical localization — to define a
**genotype-stratified natural history** of the disease. Pathogenic somatic variants
were identified in 67.4% of lesions, dominated by *KRAS*<sup>G12D</sup>; variant-positive lesions
presented and ruptured roughly two decades earlier than panel-negative bAVMs, variant
allele fraction (VAF) tracked inversely with age at rupture as a dose-dependent modifier
of severity, and endovascular-enabled liquid biopsy detected *KRAS*<sup>G12D</sup> in perinidal
cell-free DNA — proof-of-concept for in situ molecular diagnosis without surgical
resection. Together these results outline a framework for precision diagnosis,
biologically anchored prognosis, and genotype-directed therapeutic trials in sporadic
bAVMs.

This repository holds the analysis code and the figure/table/statistics outputs
underlying those findings.

## Repository layout

```
analysis/
  00_data_prep/      00 mask PHI · 01 clean the cohort master → analysis-ready dataset
  01_main_analysis/  producers + composers, named NN_<TARGET>_<desc>.R and grouped so
                     each main figure sits with its companion Extended Data figures
                     (F1–F4, ED1–ED10, T1 = Table 1, ST = Supplementary Tables)
  pipeline/          run_all.R orchestrator + figure-assembly / stats-manifest helpers
  helper_scripts/    shared R utilities (palettes, themes, plotting, stats)
results/
  figures/           manuscript-facing deliverable, organized by figure + panel:
                       Figure 1/ … Figure 4/        Fig N (composite) + Fig NA, NB, …
                       Extended Data Fig 1/ … 10/   ED Fig N (composite) + ED Fig NA, …
                       Tables/                      Table 1 + Supplementary Tables
  stats/             text readouts of the reported statistics
  _sessionInfo.txt   R session the outputs were generated under
```

No data are committed to this repository — see **Data availability** below.

`analysis/01_main_analysis/SCRIPT_INDEX.md` maps every script to the figure(s),
table(s), and statistics it produces.

## Reproducing the figures

Everything runs as one sequence from `00` to the end:

```r
# R 4.5.1
source("analysis/pipeline/builders/run_all.R")
```

`run_all.R` (1) reprocesses the cohort data, (2) runs every producer in
`analysis/01_main_analysis/` in number order — each writing its panels, tables, and a
statistics fragment, (3) aggregates the fragments into
`results/stats/manuscript_stats.rds`, (4) validates the manifest and panel-token
uniqueness, and (5) organizes the panels into the labeled `results/figures/` tree.
Figure numbers and panel letters come from the frozen
`results/stats/panel_assignments.rds`.

## Data availability

No data are committed to this repository. The cleaned cohort dataset that drives the
analysis is available from the corresponding author on reasonable request, subject to
the governing IRB and data-use agreements. With that dataset placed under `data/`,
`analysis/00_data_prep/01_clean_master.R` and the full producer chain regenerate every
figure, table, and statistic.

One extended-data panel (dPCR validation) draws its source measurements from the
companion variant-detection dataset; its rendered panel is included here.

## Citation

If you use this code, please cite the associated manuscript (Hale & Keller et al.).
