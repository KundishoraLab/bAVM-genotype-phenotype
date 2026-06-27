# 17_ED8_anatomy.R — Extended Data Figure 8 (anatomical location).
#
# Three panels backing §3 of the Results (location + most demographic nulls):
#   a) Per-lobe prevalence in mut+ vs panel-negative lesions across nine
#      anatomical territories, with Fisher exact P and BH-FDR across the 9
#      comparisons.
#   b) Anatomical composition per variant class (KRAS G12D / KRAS G12V /
#      Other KRAS / BRAF / Negative), shown as stacked proportions per lobe.
#   c) Cleveland dot plot of 27 location-by-genotype interaction tests
#      (nine lobes × three phenotypic outcomes: rupture, age, SM size).
#      The parietal × KRAS × rupture cell is highlighted as the sole
#      nominally significant test (β = 1.82, P = 0.0019).
#
# Inputs:  data/processed/bAVM_analysis_ready.rds
# Output:  results/ExtendedData/ed_anatomy/panel_A/*.{png,pdf}
#          results/ExtendedData/ed_anatomy/panel_B/*.{png,pdf}
#          results/ExtendedData/ed_anatomy/panel_C/*.{png,pdf}
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(here); library(patchwork); library(scales)
})

source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(mutation_positive))

# Nine anatomical territories (aligned with 13_F2_genotype_phenotype.R §location interactions)
loc_vars   <- c("loc_frontal", "loc_temporal", "loc_parietal", "loc_occipital",
                "loc_cerebellar", "loc_basal_ganglia", "loc_thalamus",
                "loc_brainstem", "loc_insular")
loc_labels <- c("Frontal", "Temporal", "Parietal", "Occipital",
                "Cerebellar", "Basal Ganglia", "Thalamus", "Brainstem",
                "Insular")

df <- df %>%
  mutate(
    geno_binary = factor(
      ifelse(mutation_positive, "Variant-positive", "Panel-negative"),
      levels = c("Panel-negative", "Variant-positive")),
    variant_group = case_when(
      !mutation_positive              ~ "Negative",
      mutation == "KRAS G12D"         ~ "KRAS G12D",
      mutation == "KRAS G12V"         ~ "KRAS G12V",
      mutation_gene == "KRAS"          ~ "Other KRAS",
      mutation_gene == "BRAF"          ~ "BRAF",
      TRUE                             ~ "Negative"
    ),
    variant_group = factor(variant_group,
      levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF", "Negative"))
  )

# 2026-05-19 (Phase 2 / Iteration 1): panels now route through
# slot_dir(token) so the registry merger of ed_anatomy + ed_parietal into
# ed_anatomic_localization auto-redirects the output directories without
# producer-side edits. RDS saving is now ON so 32_ED8_assemble.R can rebuild the merged composite from the per-panel
# ggplot objects rather than re-running this script's analysis.
source(here("analysis", "pipeline", "helpers", "save_panel.R"))
# Token-routed save: ignore the legacy `_subdir` argument and resolve the
# canonical slot directory from the panel token (the panel-assignments
# resolver maps token -> group -> panel_<LETTER> dir, picking up the
# merged ed_anatomic_localization group automatically).
save_panel <- function(.subdir_ignored, token, plot, w, h) {
  d <- panel_slot_dir(token)
  if (is.null(d)) {
    stop(sprintf("[17_ED8_anatomy.R] panel_slot_dir() could not resolve token '%s'.",
                 token), call. = FALSE)
  }
  save_panel_impl(d, token, plot, w, h, device = "cairo", save_rds = TRUE)
}
out_root <- here("results", "ExtendedData", "ed_anatomic_localization")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ═════════════════════════════════════════════════════════════════════════════
# Panel A — Per-lobe prevalence (mut+ vs neg) with Fisher P + BH-FDR
# ═════════════════════════════════════════════════════════════════════════════
cat("\n── Panel A: per-lobe prevalence mut+ vs neg ──\n")

panelA_df <- map2_dfr(loc_vars, loc_labels, function(v, lab) {
  s <- df %>% filter(!is.na(.data[[v]]))
  tbl <- table(s$geno_binary, s[[v]])
  mut_pct <- 100 * mean(s[[v]][s$geno_binary == "Variant-positive"] == 1)
  neg_pct <- 100 * mean(s[[v]][s$geno_binary == "Panel-negative"] == 1)
  p <- tryCatch(fisher.test(tbl)$p.value, error = function(e) NA_real_)
  tibble(lobe = lab, mut_pct = mut_pct, neg_pct = neg_pct,
         fisher_p = p, n = nrow(s))
}) %>%
  mutate(p_fdr = p.adjust(fisher_p, method = "BH"))

print(panelA_df)

panelA_long <- panelA_df %>%
  pivot_longer(c(mut_pct, neg_pct), names_to = "grp", values_to = "pct") %>%
  mutate(grp = recode(grp,
                      "mut_pct" = "Variant-positive",
                      "neg_pct" = "Panel-negative"),
         grp = factor(grp,
                      levels = c("Variant-positive", "Panel-negative")),
         lobe = factor(lobe, levels = rev(loc_labels)))

pA <- ggplot(panelA_long, aes(x = grp, y = lobe, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            size = NM$text$body_mm) +
  # Anatomy-locked Carbon Brown gradient (PAL_ANATOMY) replaces the prior
  # blue heatmap ramp so all anatomy-themed panels share one hue family
  # distinct from the blue/green/teal/purple used elsewhere.
  scale_fill_gradient(low = PAL_ANATOMY$LOW, high = PAL_ANATOMY$HIGH,
                      limits = c(0, max(panelA_long$pct, na.rm = TRUE)),
                      name = "Prevalence (%)") +
  labs(x = NULL, y = NULL) +
  theme_nature_panel() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

save_panel("panel_A", "per_lobe_prevalence", pA, 3.30, 1.40)

# ═════════════════════════════════════════════════════════════════════════════
# Panel B — Per-variant anatomical composition (stacked bars)
# ═════════════════════════════════════════════════════════════════════════════
cat("\n── Panel B: per-variant anatomical composition ──\n")

panelB_df <- df %>%
  select(variant_group, all_of(loc_vars)) %>%
  pivot_longer(cols = all_of(loc_vars), names_to = "loc", values_to = "hit") %>%
  filter(!is.na(hit), !is.na(variant_group)) %>%
  mutate(loc_label = loc_labels[match(loc, loc_vars)],
         loc_label = factor(loc_label, levels = rev(loc_labels))) %>%
  group_by(variant_group, loc_label) %>%
  summarise(pct = 100 * mean(hit == 1), n = n(), .groups = "drop")

# v6.36 (2026-05-21): drop the "\n(n = N)" annotation — the 2-line
# tick labels overflow the narrow Nature-spec cell and the per-stratum
# n is already disclosed in the figure caption + SuppTable. Display
# only the variant name; "Negative" gets the manuscript-wide
# "Panel-negative" relabel.
variant_n <- df %>% count(variant_group) %>%
  mutate(col_label = dplyr::recode(
    as.character(variant_group),
    "Negative" = "Panel-negative"
  ))

panelB_plot_df <- panelB_df %>%
  left_join(variant_n, by = "variant_group") %>%
  mutate(col_label = factor(col_label,
                            levels = variant_n$col_label[
                              order(variant_n$variant_group)]))

pB <- ggplot(panelB_plot_df,
             aes(x = col_label, y = loc_label, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            size = NM$text$body_mm) +
  # Anatomy-locked Carbon Brown gradient (PAL_ANATOMY).
  scale_fill_gradient(low = PAL_ANATOMY$LOW, high = PAL_ANATOMY$HIGH,
                      limits = c(0, max(panelB_plot_df$pct, na.rm = TRUE)),
                      name = "Prevalence (%)") +
  labs(x = NULL, y = NULL) +
  theme_nature_panel() +
  # Angle the x-axis tick labels so the multi-line "<variant>\n(n = NN)"
  # labels stop colliding when rendered side-by-side at composite width.
  # angle=30 + hjust=1 gives clean diagonal placement without overflowing
  # the panel footprint.
  theme(legend.position = "right",
        # v6.37 (2026-05-21): angle 30 -> 45 — at narrow Nature-spec
        # cell width the 30-degree labels still overflowed/collided.
        axis.text.x = element_text(angle = 45, hjust = 1, lineheight = 0.9))

save_panel("panel_B", "per_variant_anatomy", pB, 3.30, 1.40)

# ═════════════════════════════════════════════════════════════════════════════
# Panel C — Cleveland dot of 27 interaction tests
# ═════════════════════════════════════════════════════════════════════════════
cat("\n── Panel C: 27 location × genotype interaction tests ──\n")

int_outcomes <- tribble(
  ~var,               ~label,            ~type,
  "ever_ruptured_num", "Rupture",         "binary",
  "age",               "Age",             "continuous",
  "sm_size_num",       "Spetzler\u2013Martin size", "continuous"
)

# Fit glm(outcome ~ geno_binary * lobe) for each (lobe, outcome) cell.
# AUDIT 2026-05-12 (F11 deep, canonical helper): the per-cell GLM is now
# fit by `fit_loc_geno_interaction()` from `analysis/helper_scripts/utils.R`,
# shared with `13_F2_genotype_phenotype.R` (SuppTable08) so the raw interaction P
# and beta sign at any overlapping (lobe, outcome) cell match exactly. The
# helper filters internally on !is.na(outcome) & !is.na(lobe) & !is.na(
# geno_binary), and relevels geno_binary with Panel-negative as the
# explicit reference category so the interaction beta sign is convention-
# stable across producers.
panelC_df <- map2_dfr(loc_vars, loc_labels, function(v, lab) {
  map_dfr(seq_len(nrow(int_outcomes)), function(i) {
    o   <- int_outcomes[i, ]
    fit <- fit_loc_geno_interaction(df, o$var, v, o$type)
    tibble(lobe    = lab,
           outcome = o$label,
           beta    = fit$beta,
           se      = fit$se,
           p       = fit$p,
           n       = fit$n)
  })
}) %>%
  mutate(p_bonf = pmin(1, p * n()),
         p_fdr  = p.adjust(p, method = "BH"),
         log10p = -log10(pmax(p, .Machine$double.eps)),
         is_parietal_rupture = (lobe == "Parietal" & outcome == "Rupture"),
         outcome = factor(outcome, levels = int_outcomes$label),
         lobe = factor(lobe, levels = rev(loc_labels)))

print(panelC_df %>% select(lobe, outcome, beta, p, p_bonf))

alpha_bonf <- 0.05 / nrow(panelC_df)

pC <- ggplot(panelC_df, aes(x = log10p, y = lobe, color = is_parietal_rupture)) +
  geom_segment(aes(xend = 0, yend = lobe), linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2) +
  ref_vline(-log10(0.05),       kind = "threshold") +
  ref_vline(-log10(alpha_bonf), kind = "bonf") +
  # Anatomy-locked palette (PAL_ANATOMY) — dots use the brown ACCENT so the
  # Cleveland panel reads as part of the anatomy figure family. The
  # parietal × rupture call-out uses PAL_ANATOMY$HIGHLIGHT (purple) for
  # contrast that sits outside Tier A (no KRAS/BRAF colour collision).
  scale_color_manual(values = c(`TRUE`  = PAL_ANATOMY$HIGHLIGHT,
                                `FALSE` = PAL_ANATOMY$ACCENT),
                     guide = "none") +
  facet_wrap(~ outcome, nrow = 1) +
  labs(x = expression(-log[10](italic(P))~"(interaction)"), y = NULL) +
  theme_nature_panel() +
  theme(panel.spacing = unit(1.2, "lines"))

# interaction_cleveland dropped from figure in v6.36 (2026-05-21); pC
# analysis kept for the panelC_df stats fragment and CSV export below.

# Composite assembly moved to 32_ED8_assemble.R as of
# the 2026-05-19 Phase-2 ed_anatomy + ed_parietal merge. This producer now
# only writes the three anatomy panel RDS / PNG / PDF triples; the
# six-panel merged composite is built downstream from the RDS files.

# Dump numeric back-up as CSV for reviewer audit
write.csv(panelA_df,
          file.path(out_root, "per_lobe_prevalence.csv"),
          row.names = FALSE)
write.csv(panelB_df,
          file.path(out_root, "per_variant_anatomy.csv"),
          row.names = FALSE)
write.csv(panelC_df %>% select(-is_parietal_rupture, -log10p),
          file.path(out_root, "interaction_tests.csv"),
          row.names = FALSE)

# ═════════════════════════════════════════════════════════════════════════════
# Stats manifest fragment — edfig08_anatomy
# -----------------------------------------------------------------------------
# §3 prose cites "lobe-level prevalence did not differ by genotype at any of
# the nine anatomical territories examined (all BH-adjusted P ≥ 0.83)" —
# that's panelA_df. Also cites the interaction test count (27 = 9 lobes ×
# 3 outcomes, with parietal × rupture the sole nominally significant cell)
# — that's panelC_df. Emit both scalars + tables.
# ═════════════════════════════════════════════════════════════════════════════
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

.panelC_nominal <- panelC_df[!is.na(panelC_df$p) & panelC_df$p < 0.05, , drop = FALSE]

edfig08_anatomy_fragment <- list(
  # Panel A — per-lobe Fisher + BH-FDR aggregate
  n_lobes              = nrow(panelA_df),
  panelA_min_fisher_p  = min(panelA_df$fisher_p, na.rm = TRUE),
  panelA_min_p_fdr     = min(panelA_df$p_fdr,    na.rm = TRUE),
  panelA_max_p_fdr     = max(panelA_df$p_fdr,    na.rm = TRUE),
  # Panel C — interaction-test aggregates
  n_interaction_tests  = nrow(panelC_df),
  # Min BH-FDR + Bonferroni q across the 27-cell location × outcome grid
  # (used by §4 prose to honestly cite the corrected significance of the
  # parietal × rupture cell rather than the panel-A prevalence FDR).
  panelC_min_p_fdr     = min(panelC_df$p_fdr,  na.rm = TRUE),
  panelC_min_p_bonf    = min(panelC_df$p_bonf, na.rm = TRUE),
  n_outcome_families   = length(unique(panelC_df$outcome)),
  n_tests_nominal_sig  = nrow(.panelC_nominal),
  # Full tables for ED caption / SI spreadsheet use
  per_lobe_prevalence  = as.data.frame(panelA_df),
  per_variant_anatomy  = as.data.frame(panelB_df),
  interaction_tests    = as.data.frame(panelC_df %>% dplyr::select(-dplyr::any_of(c("is_parietal_rupture", "log10p"))))
)
write_stats_section(section = "edfig08_anatomy",
                    stats   = edfig08_anatomy_fragment)

cat("\n══ 17_ED8_anatomy.R complete ══\n")
cat(sprintf("   Panel A:     %s\n", file.path(out_root, "panel_A")))
cat(sprintf("   Panel B:     %s\n", file.path(out_root, "panel_B")))
cat(sprintf("   Panel C:     %s\n", file.path(out_root, "panel_C")))
cat(sprintf("   Composite:   %s\n",
            file.path(here("results", "ExtendedData", "ed_anatomic_localization"))))
