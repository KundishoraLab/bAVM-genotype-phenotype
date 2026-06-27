# 14_ED7_vaf_phenotype.R — VAF × phenotype null panels feed.
#
# Registry tokens written here: ed_vaf_phenotype (6 panels: vaf_sm_total,
# vaf_sm_size, vaf_drainage, vaf_eloquence, vaf_rupture, vaf_highrisk) plus
# panel feeds for ed_per_variant_pheno. Manuscript figure number is
# resolver-driven; do not hard-code an ED-N label here or in the caption.
#
# Input:  data/processed/bAVM_analysis_ready.rds
# Output: results/ExtendedData/{ed_vaf_phenotype, ed_per_variant_pheno, ...}/
#           per-panel pdf + png + rds (rds added 2026-04-26 so the
#           ed_vaf_phenotype composer can read panels via compose_figure)
#         + manifest fragment edfig_vaf_phenotype (see stats_schema.R)
#
# VAF analyses restricted to variant-positive patients with VAF data (n ≈ 144).
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(broom)
library(ordinal)  # proportional odds logistic regression for ordinal outcomes

source(here("analysis", "helper_scripts", "utils.R"))

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
# Force cairo device — Arial doesn't embed via the default grDevices::pdf
# device (see ED02 missingness-heatmap fix in 02_prep_analysis_dataset.R).
save_panel <- function(dir, name, plot, w, h)
  save_panel_impl(dir, name, plot, w, h, save_rds = FALSE, device = "cairo")
# NB: line ~155 of this file redefines save_panel() with a different
# signature (dir_path, base_name, pooled, ...) for the pooled-fit panels.
# That redefinition is preserved as-is — its semantics differ from the
# canonical helper (NULL-pooled short-circuit) and the redefine kicks in
# at the right point in the script's flow.

# ── 1. Load data ─────────────────────────────────────────────────────────────

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))
palettes <- readRDS(here("results", "palettes.rds"))
output_dir <- here("results")
efig_dir <- file.path(output_dir, "ExtendedData")
# 2026-05-19 (Phase 2 / Iteration 3): ed_vaf_phenotype merged into the
# umbrella ed_vaf_deep_dive group (alongside vaf_outlier_combined and
# vaf_age_rupture_scatter). edvp_dir now points at ed_vaf_deep_dive/;
# panels remain saved flat at the top level so the merged composer
# (31_ED7_assemble.R) can read them via <token>.rds.
edvp_dir <- file.path(efig_dir, "ed_vaf_deep_dive")

# VAF analyses: variant-positive with VAF.
# relabel_geno_factor() swaps any residual "Negative" factor level to the
# canonical "Panel-negative" display label so facet strips and legends in
# every panel below honour the manuscript-wide display contract — a single
# canonical fix avoids per-ggplot patches at lines 69/81/97/113/130/145/217/240.
# PAL_VARIANT is dual-keyed (both raw and display strings → same hex), so
# scale_color_manual(values = PAL_VARIANT) lookup still resolves correctly.
vaf_data <- df %>%
  filter(mutation_positive == TRUE & !is.na(vaf_prop)) %>%
  mutate(vaf_pct = vaf_prop * 100,
         geno_variant = relabel_geno_factor(geno_variant))

vaf_data <- vaf_data %>%
  mutate(gene_group = case_when(
    str_detect(geno_variant, "KRAS") ~ "KRAS",
    str_detect(geno_variant, "BRAF") ~ "BRAF",
    TRUE ~ "Other"
  ),
  # Variant identity reduced to four buckets for SHAPE encoding on every
  # ED10 jittered-dot panel: KRAS G12D / KRAS G12V / Other KRAS / BRAF.
  # Replaces the previous colour-by-gene mapping (blue=KRAS, pink=BRAF)
  # so the visual reads as a single-colour scatter and the symbol carries
  # variant identity — printable in black-and-white and consistent with
  # the rest of the manuscript's colour usage.
  variant_shape_group = case_when(
    str_detect(as.character(geno_variant), "G12D")        ~ "KRAS G12D",
    str_detect(as.character(geno_variant), "G12V")        ~ "KRAS G12V",
    str_detect(as.character(geno_variant), "^KRAS")       ~ "Other KRAS",
    str_detect(as.character(geno_variant), "BRAF")        ~ "BRAF",
    TRUE                                                   ~ "Other"
  ),
  variant_shape_group = factor(variant_shape_group,
    levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF", "Other")))

PAL_GENE <- c("KRAS" = unname(PAL_VARIANT["KRAS G12D"]), "BRAF" = unname(PAL_VARIANT["BRAF"]))

# Shape contract for VAF scatter panels — symbols are colour-blind safe and
# the open/filled circle pair encodes the two dominant KRAS variants. Other
# KRAS gets an open triangle so it stays visually distinct from G12V's open
# circle; BRAF uses the X mark to read as categorically different from any
# KRAS shape.
SHAPE_VARIANT <- c(
  "KRAS G12D"  = 16,   # filled circle
  "KRAS G12V"  = 1,    # open circle
  "Other KRAS" = 2,    # open triangle
  "BRAF"       = 4,    # X mark
  "Other"      = 8     # asterisk (rare; kept as fallback)
)
VARIANT_DOT_COLOUR <- "#1A1A1A"  # paper-black; uniform across every panel

cat(sprintf("VAF analysis cohort: %d variant-positive patients with VAF\n", nrow(vaf_data)))

# ── Helper: faceted scatter by mutation ──────────────────────────────────────
facet_scatter_mutation <- function(data, x_var, y_var, x_lab, y_lab, title) {
  sub <- data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  if (nrow(sub) < 5) return(NULL)
  ggplot(sub, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = geno_variant), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, linetype = "dashed") +
    facet_wrap(~geno_variant, scales = "free_x") +
    scale_color_manual(values = PAL_VARIANT, guide = "none") +
    labs(title = paste0(title, " (by mutation)"), x = x_lab, y = y_lab) +
    theme_nature_panel()
}

facet_scatter_gene <- function(data, x_var, y_var, x_lab, y_lab, title) {
  sub <- data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  if (nrow(sub) < 5) return(NULL)
  ggplot(sub, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = gene_group), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, linetype = "dashed") +
    facet_wrap(~gene_group) +
    scale_color_manual(values = PAL_GENE, guide = "none") +
    labs(title = paste0(title, " (by gene)"), x = x_lab, y = y_lab) +
    theme_nature_panel()
}

# ── Helper: faceted boxplot by mutation ──────────────────────────────────────
facet_box_mutation <- function(data, binary_var, var_label, title) {
  sub <- data %>% filter(!is.na(.data[[binary_var]]))
  if (nrow(sub) < 5 || length(unique(sub[[binary_var]])) < 2) return(NULL)
  ggplot(sub, aes(x = factor(.data[[binary_var]]), y = vaf_pct)) +
    geom_violin(aes(fill = factor(.data[[binary_var]])), alpha = 0.3,
      scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(color = geno_variant), width = 0.1, size = 1.5, alpha = 0.6) +
    facet_wrap(~geno_variant) +
    scale_fill_manual(values = c("0" = unname(PAL_BINARY[["Panel-negative"]]), "1" = unname(PAL_BINARY[["Variant-positive"]])), guide = "none") +
    scale_color_manual(values = PAL_VARIANT, guide = "none") +
    scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
    labs(title = paste0(title, " (by mutation)"), x = var_label, y = "VAF (%)") +
    theme_nature_panel()
}

facet_box_gene <- function(data, binary_var, var_label, title) {
  sub <- data %>% filter(!is.na(.data[[binary_var]]))
  if (nrow(sub) < 5 || length(unique(sub[[binary_var]])) < 2) return(NULL)
  ggplot(sub, aes(x = factor(.data[[binary_var]]), y = vaf_pct)) +
    geom_violin(aes(fill = factor(.data[[binary_var]])), alpha = 0.3,
      scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(color = gene_group), width = 0.1, size = 1.5, alpha = 0.6) +
    facet_wrap(~gene_group) +
    scale_fill_manual(values = c("0" = unname(PAL_BINARY[["Panel-negative"]]), "1" = unname(PAL_BINARY[["Variant-positive"]])), guide = "none") +
    scale_color_manual(values = PAL_GENE, guide = "none") +
    scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
    labs(title = paste0(title, " (by gene)"), x = var_label, y = "VAF (%)") +
    theme_nature_panel()
}

# ── Helper: faceted ordinal violin+boxplot by mutation ─────────────────────
facet_ordinal_mutation <- function(data, ordinal_var, var_label, title) {
  sub <- data %>% filter(!is.na(.data[[ordinal_var]]))
  if (nrow(sub) < 5 || length(unique(sub[[ordinal_var]])) < 2) return(NULL)
  sub$ord_factor <- factor(sub[[ordinal_var]])
  ggplot(sub, aes(x = ord_factor, y = vaf_pct)) +
    geom_violin(aes(fill = ord_factor), alpha = 0.3, scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(color = geno_variant), width = 0.1, size = 1.5, alpha = 0.6) +
    facet_wrap(~geno_variant) +
    scale_fill_brewer(palette = "Blues", guide = "none") +
    scale_color_manual(values = PAL_VARIANT, guide = "none") +
    labs(title = paste0(title, " (by mutation)"), x = var_label, y = "VAF (%)") +
    theme_nature_panel()
}

facet_ordinal_gene <- function(data, ordinal_var, var_label, title) {
  sub <- data %>% filter(!is.na(.data[[ordinal_var]]))
  if (nrow(sub) < 5 || length(unique(sub[[ordinal_var]])) < 2) return(NULL)
  sub$ord_factor <- factor(sub[[ordinal_var]])
  ggplot(sub, aes(x = ord_factor, y = vaf_pct)) +
    geom_violin(aes(fill = ord_factor), alpha = 0.3, scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(color = gene_group), width = 0.1, size = 1.5, alpha = 0.6) +
    facet_wrap(~gene_group) +
    scale_fill_brewer(palette = "Blues", guide = "none") +
    scale_color_manual(values = PAL_GENE, guide = "none") +
    labs(title = paste0(title, " (by gene)"), x = var_label, y = "VAF (%)") +
    theme_nature_panel()
}

# ── Helper: save the pooled panel ────────────────────────────────────────────
# 2026-04-26: simplified from the legacy 3-arg trio (pooled / by_mut /
# by_gene). Only the pooled panel was ever embedded — by_mut and by_gene
# were computed by the call sites and dropped silently. The call sites
# below were trimmed to pass only the pooled panel; the dead facet
# variants would need to be revived as separate ED panels if reactivated.
# Adds .rds persistence so 31_ED7_assemble.R can read each
# panel without re-fitting the analysis.
# v6.29 (2026-05-20): native ED07 cell dims (6.60 / 2 cols = 3.0 wide
# x 6.09 / 4 rows = 1.6 in tall — close to Fig 1 Panel C/G size). Also
# resolves each panel's canonical slot via panel_slot_dir() so the
# sync_panel_prefixes "registry-misplaced" sweep doesn't delete the
# group-root copies the producer used to drop here.
source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa_vaf <- load_panel_assignments()
save_panel <- function(dir_path, base_name, pooled, w = 3.0, h = 1.6) {
  if (is.null(pooled)) return(invisible(NULL))
  slot <- tryCatch(panel_slot_dir(base_name, .pa_vaf), error = function(e) dir_path)
  if (is.null(slot) || identical(slot, "")) slot <- dir_path
  dir_path <- slot
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir_path, paste0(base_name, ".pdf")), pooled,
         width = w, height = h, device = grDevices::cairo_pdf,
         family = NM$font_family)
  ggsave(file.path(dir_path, paste0(base_name, ".png")), pooled,
         width = w, height = h, dpi = 300, type = "cairo")
  saveRDS(pooled, file.path(dir_path, paste0(base_name, ".rds")))
  invisible(NULL)
}

# ── Helper: Spearman scatter with correlation ────────────────────────────────

vaf_scatter <- function(data, x_var, y_var, x_lab, y_lab, title, subtitle_extra = "") {
  sub <- data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  if (nrow(sub) < 5) {
    cat(sprintf("  Skipping %s: only %d observations\n", title, nrow(sub)))
    return(NULL)
  }
  cor_test <- cor.test(sub[[x_var]], sub[[y_var]], method = "spearman", exact = FALSE)
  p <- ggplot(sub, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.6, size = 2, color = PAL_VARIANT["KRAS G12D"]) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, linetype = "dashed") +
    # title/subtitle stripped — captions carry n, rho, p
    labs(x = x_lab, y = y_lab) +
    theme_nature_panel()
  list(plot = p, cor = cor_test, n = nrow(sub))
}

# ── Helper: violin + jitter for VAF by ordinal category ──────────────────────
# Uses proportional odds (cumulative link) model for the statistical test
vaf_ordinal <- function(data, ordinal_var, var_label, title) {
  sub <- data %>% filter(!is.na(.data[[ordinal_var]]))
  if (nrow(sub) < 10 || length(unique(sub[[ordinal_var]])) < 2) {
    cat(sprintf("  Skipping %s: n=%d, levels=%d\n", title, nrow(sub),
      length(unique(sub[[ordinal_var]]))))
    return(NULL)
  }
  sub$ord_factor <- factor(sub[[ordinal_var]])
  # Proportional odds logistic regression: ordinal outcome ~ VAF
  clm_mod <- tryCatch(
    clm(ord_factor ~ vaf_pct, data = sub),
    error = function(e) NULL
  )
  if (!is.null(clm_mod)) {
    cs <- summary(clm_mod)$coefficients
    vaf_row <- grep("vaf_pct", rownames(cs))
    clm_beta <- cs[vaf_row, 1]
    clm_p <- cs[vaf_row, 4]
    sub_text <- sprintf("n = %d | Proportional odds: \u03B2 = %.3f, p = %.3f", nrow(sub), clm_beta, clm_p)
  } else {
    clm_beta <- NA; clm_p <- NA
    sub_text <- sprintf("n = %d", nrow(sub))
  }
  # Spearman as supplementary
  sp <- cor.test(sub$vaf_pct, sub[[ordinal_var]], method = "spearman", exact = FALSE)
  p <- ggplot(sub, aes(x = ord_factor, y = vaf_pct)) +
    # v6.28 (2026-05-20): Fig 1 Panel C 1-to-1 polish — alpha 0.4 -> 0.2,
    # violin linewidth 0.4, boxplot linewidth 0.3, jitter size 1.8 -> 0.7
    # so panels render cleanly at the new 2.20 x 1.545 in cell footprint.
    geom_violin(aes(fill = ord_factor), alpha = 0.2, scale = "width", trim = FALSE, linewidth = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, linewidth = 0.3) +
    geom_jitter(aes(shape = variant_shape_group), colour = VARIANT_DOT_COLOUR,
                width = 0.1, size = 0.7, alpha = 0.5, stroke = 0.3) +
    scale_fill_brewer(palette = "Blues", guide = "none") +
    scale_shape_manual(values = SHAPE_VARIANT, name = "Variant",
                       drop = FALSE) +
    # Two-row legend keeps the shape key within panel A's column width in
    # the ED10 composite (5 variants in one row was clipping past the
    # left edge of the figure canvas).
    guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
    # title/subtitle stripped — caption carries n, POLR β + p
    labs(x = var_label, y = "VAF (%)") +
    theme_nature_panel() +
    theme(legend.justification = "left")
  list(plot = p, clm_beta = clm_beta, clm_p = clm_p, spearman = sp, n = nrow(sub))
}

# ── Helper: VAF violin+jitter for binary outcome ─────────────────────────────

vaf_boxplot <- function(data, binary_var, var_label, title) {
  sub <- data %>% filter(!is.na(.data[[binary_var]]))
  if (nrow(sub) < 5 || length(unique(sub[[binary_var]])) < 2) {
    cat(sprintf("  Skipping %s: n=%d, levels=%d\n", title, nrow(sub),
      length(unique(sub[[binary_var]]))))
    return(NULL)
  }
  wt <- wilcox.test(sub$vaf_pct ~ sub[[binary_var]])
  p <- ggplot(sub, aes(x = factor(.data[[binary_var]]), y = vaf_pct)) +
    # v6.28 (2026-05-20): Fig 1 Panel C 1-to-1 polish (see vaf_ordinal above).
    geom_violin(aes(fill = factor(.data[[binary_var]])), alpha = 0.2,
      scale = "width", trim = FALSE, linewidth = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, linewidth = 0.3) +
    geom_jitter(aes(shape = variant_shape_group), colour = VARIANT_DOT_COLOUR,
                width = 0.1, size = 0.7, alpha = 0.5, stroke = 0.3) +
    # Drop Tier B Mut+ green for the "Present" violin — Mut+ green is
    # reserved for the cohort split, not for "outcome-present" inside the
    # variant-positive sub-cohort. Use a sequential blue pair (light/dark
    # heat blues) so binary violins match the Blues fill used by the
    # ordinal violins (A / F).
    scale_fill_manual(values = c("0" = "#D1E5F0", "1" = "#4393C3"), guide = "none") +
    scale_shape_manual(values = SHAPE_VARIANT, name = "Variant",
                       drop = FALSE) +
    scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
    # title/subtitle stripped — caption carries n + Wilcoxon p
    labs(x = var_label, y = "VAF (%)") +
    theme_nature_panel()
  list(plot = p, wilcox = wt, n = nrow(sub))
}

# ── Helper: variant × phenotype-level heatmap of median VAF ──────────────
# Used by ED10 panels B (vaf_sm_size), C (vaf_drainage), D (vaf_eloquence)
# in place of the prior violin+jitter layout. Rows = phenotype levels,
# columns = variant_shape_group (KRAS G12D / G12V / Other KRAS / BRAF),
# cell fill = median VAF (%), in-cell label = "<med%>\nN=<n>".
#
# 2026-05-17 — switched from the heat-blue ramp to white → PAL_HEAT_HIGH
# (= PAL_SM[["IV"]] = #2166AC) so the heatmap fill family matches the
# SM Total Score Brewer Blues violins in panel A. Contrast-aware text
# colour (white on dark cells, black on light) + bumped + bolded in-cell
# label so the value reads through the gradient. (Previously briefly
# used the Fig 2G white → PAL_SCORE[["2"]] purple ramp; reverted to keep
# every panel that touches an SM-tier feature in the same blue family.)
vaf_heatmap <- function(data, level_var, level_label,
                        x_labels = NULL,   # optional per-level relabel
                        title    = NULL) {
  sub <- data %>%
    filter(!is.na(.data[[level_var]]),
           !is.na(variant_shape_group),
           variant_shape_group != "Other")
  if (nrow(sub) < 5) return(NULL)
  variant_levels <- c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF")
  agg <- sub %>%
    mutate(level = .data[[level_var]],
           level_str = as.character(level),
           level_str = if (!is.null(x_labels))
             dplyr::recode(level_str, !!!x_labels)
             else level_str) %>%
    group_by(variant_shape_group, level, level_str) %>%
    summarise(median_vaf = stats::median(vaf_pct, na.rm = TRUE),
              n = dplyr::n(),
              .groups = "drop") %>%
    arrange(level)
  # Ensure all four variant columns appear in every heatmap (KRAS G12D /
  # G12V / Other KRAS / BRAF) — tidyr::complete() fills any combination
  # with no patients as NA so the column-axis stays uniform across the
  # B/C/D panels even when a rare variant has zero rows at some level.
  level_str_order <- agg %>% distinct(level, level_str) %>%
    arrange(level) %>% pull(level_str)
  agg <- agg %>%
    mutate(variant_shape_group = factor(variant_shape_group,
                                        levels = variant_levels),
           level_str = factor(level_str, levels = level_str_order)) %>%
    tidyr::complete(variant_shape_group, level_str,
                    fill = list(median_vaf = NA_real_, n = 0L))
  # In-cell annotation text is always black (per author 2026-05-27) — kept
  # uniform across all ED heatmaps rather than switching to white on dark
  # cells. fill_max still bounds the fill ramp (scale_fill_gradient limits).
  fill_max <- max(agg$median_vaf, na.rm = TRUE)
  agg <- agg %>% mutate(text_col = "black")
  p <- ggplot(agg, aes(x = variant_shape_group, y = level_str,
                       fill = median_vaf)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # N= dropped 2026-05-17 — % alone reads cleanly at composite scale;
    # per-variant N is already disclosed elsewhere (panel-level captions,
    # SuppTable 4). Same edit applied to ED11's prevalence_heatmap.
    geom_text(aes(label = ifelse(n == 0L, "—",
                                 sprintf("%.1f%%", median_vaf)),
                  colour = text_col),
              # v6.31 (2026-05-20): size 4.5 -> 2.0 — mirrors Fig 1 Panel I
              # (rupture_lookup_heatmap) which uses size = NM$text$body_mm (~5.7 pt) to
              # stay inside Nature's 5-7 pt body band at native cell footprint.
              size = NM$text$body_mm, fontface = "bold") +
    scale_fill_gradient(low = "white",
                        high = PAL_HEAT_HIGH,
                        na.value = "grey92",
                        limits = c(0, fill_max),
                        name = "Median\nVAF (%)") +
    scale_colour_identity() +
    labs(x = NULL, y = level_label, title = title) +
    theme_nature_panel() +
    # Legend suppressed: the in-cell "<med%>\nN=<n>" label carries the
    # value directly, so a separate fill key would just steal width from
    # the heatmap when packed into the ED10 composite's 3-column row.
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 30, hjust = 1,
                                         lineheight = 0.9),
          panel.grid      = element_blank())
  list(plot = p, n = nrow(sub))
}

# ── Collect all stats for Supplementary Table 4 ──────────────────────────────────────────

all_stats <- list()

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3A: VAF × SM total score
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_vaf_phenotype: vaf_sm_total ──\n")
fig3a <- vaf_ordinal(vaf_data, "sm_total_num", "SM Total Score",
  "VAF by SM Total Score")
if (!is.null(fig3a)) {
  fig3a_gene <- facet_ordinal_gene(vaf_data, "sm_total_num", "SM Total Score",
    "VAF vs SM Total Score")
  save_panel(edvp_dir, "vaf_sm_total", fig3a$plot)
  all_stats[["VAF_vs_SM_total"]] <- list(rho = fig3a$spearman$estimate,
    p = fig3a$spearman$p.value, clm_p = fig3a$clm_p, n = fig3a$n)
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3B–D: VAF × SM components
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_vaf_phenotype: vaf_sm_size ──\n")
# Panels B/C/D rendered as median-VAF heatmaps (variant × phenotype level).
# Stat (POLR β / Wilcoxon p, n) is still computed from the existing
# ordinal/boxplot helper and pushed into `all_stats` so SuppTable 4 stays
# unchanged; only the visual representation switches to heatmap.
fig3b <- vaf_ordinal(vaf_data, "sm_size_num", "SM Size Score",
  "VAF by SM Size")
fig3b_hm <- vaf_heatmap(vaf_data, "sm_size_num", "SM Size Score")
if (!is.null(fig3b) && !is.null(fig3b_hm)) {
  fig3b_gene <- facet_ordinal_gene(vaf_data, "sm_size_num", "SM Size Score",
    "VAF vs SM Size")
  save_panel(edvp_dir, "vaf_sm_size", fig3b_hm$plot)
  all_stats[["VAF_vs_SM_size"]] <- list(rho = fig3b$spearman$estimate,
    p = fig3b$spearman$p.value, clm_p = fig3b$clm_p, n = fig3b$n)
}

cat("\n── ed_vaf_phenotype: vaf_drainage ──\n")
fig3c <- vaf_boxplot(vaf_data, "sm_drainage_num", "Deep Drainage",
  "VAF vs Deep Venous Drainage")
fig3c_hm <- vaf_heatmap(vaf_data, "sm_drainage_num", "Deep Drainage",
  x_labels = c("0" = "Absent", "1" = "Present"))
if (!is.null(fig3c) && !is.null(fig3c_hm)) {
  fig3c_mut <- facet_box_mutation(vaf_data, "sm_drainage_num", "Deep Drainage",
    "VAF vs Deep Venous Drainage")
  fig3c_gene <- facet_box_gene(vaf_data, "sm_drainage_num", "Deep Drainage",
    "VAF vs Deep Venous Drainage")
  save_panel(edvp_dir, "vaf_drainage", fig3c_hm$plot)
  all_stats[["VAF_vs_drainage"]] <- list(W = fig3c$wilcox$statistic, p = fig3c$wilcox$p.value, n = fig3c$n)
}

cat("\n── ed_vaf_phenotype: vaf_eloquence ──\n")
fig3d <- vaf_boxplot(vaf_data, "sm_eloquence_num", "Eloquent Location",
  "VAF vs Eloquent Location")
fig3d_hm <- vaf_heatmap(vaf_data, "sm_eloquence_num", "Eloquent Location",
  x_labels = c("0" = "Absent", "1" = "Present"))
if (!is.null(fig3d) && !is.null(fig3d_hm)) {
  fig3d_mut <- facet_box_mutation(vaf_data, "sm_eloquence_num", "Eloquent Location",
    "VAF vs Eloquent Location")
  fig3d_gene <- facet_box_gene(vaf_data, "sm_eloquence_num", "Eloquent Location",
    "VAF vs Eloquent Location")
  save_panel(edvp_dir, "vaf_eloquence", fig3d_hm$plot)
  all_stats[["VAF_vs_eloquence"]] <- list(W = fig3d$wilcox$statistic, p = fig3d$wilcox$p.value, n = fig3d$n)
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3E: VAF × composite high-risk count
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_vaf_phenotype: vaf_highrisk ──\n")
fig3e <- vaf_ordinal(vaf_data, "n_high_risk_num", "Number of High-Risk Features",
  "VAF by High-Risk Feature Count")
if (!is.null(fig3e)) {
  fig3e_gene <- facet_ordinal_gene(vaf_data, "n_high_risk_num", "Number of High-Risk Features",
    "VAF vs High-Risk Feature Count")
  # Drop the redundant shape legend on this panel — panel A (vaf_sm_total)
  # carries the single canonical Variant legend in the ED10 composite.
  save_panel(edvp_dir, "vaf_highrisk",
             fig3e$plot + theme(legend.position = "none"))
  all_stats[["VAF_vs_high_risk"]] <- list(rho = fig3e$spearman$estimate,
    p = fig3e$spearman$p.value, clm_p = fig3e$clm_p, n = fig3e$n)
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3G: VAF × rupture
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── ed_vaf_phenotype: vaf_rupture ──\n")
fig3g <- vaf_boxplot(vaf_data, "ever_ruptured_num", "Ever Ruptured",
  "VAF vs Rupture Status")
if (!is.null(fig3g)) {
  fig3g_mut <- facet_box_mutation(vaf_data, "ever_ruptured_num", "Ever Ruptured",
    "VAF vs Rupture Status")
  fig3g_gene <- facet_box_gene(vaf_data, "ever_ruptured_num", "Ever Ruptured",
    "VAF vs Rupture Status")
  # See vaf_highrisk note above — single canonical legend on panel A.
  save_panel(edvp_dir, "vaf_rupture",
             fig3g$plot + theme(legend.position = "none"))
  # Logistic regression: rupture ~ VAF + SM size + drainage + age + sample_type
  rupt_vaf_data <- vaf_data %>%
    filter(!is.na(ever_ruptured_num) & !is.na(sm_size_num) & !is.na(sm_drainage_num) &
           !is.na(age) & !is.na(sample_type_clean))
  if (nrow(rupt_vaf_data) >= 20) {
    rupt_vaf_mod <- glm(ever_ruptured_num ~ vaf_pct + sm_size_num + sm_drainage_num +
      age + sample_type_clean,
      data = rupt_vaf_data, family = binomial)
    cat("Rupture ~ VAF (adjusted incl. sample type) model:\n")
    print(summary(rupt_vaf_mod)$coefficients)
    all_stats[["rupture_vaf_adjusted"]] <- tidy(rupt_vaf_mod, conf.int = TRUE, exponentiate = TRUE)
  }
  all_stats[["VAF_vs_rupture"]] <- list(W = fig3g$wilcox$statistic, p = fig3g$wilcox$p.value, n = fig3g$n)
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3H: VAF × growth (PLACEHOLDER — only 2 genotyped patients have data)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 3H: VAF × growth (placeholder) ──\n")
growth_vaf <- vaf_data %>% filter(!is.na(growing_num))
cat(sprintf("  Growth data available: %d patients\n", nrow(growth_vaf)))

fig3h_placeholder <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
    label = paste0("VAF vs Growth\n\n",
      sprintf("Only %d genotyped patients have growth data.\n", nrow(growth_vaf)),
      "145 additional patients with growth data\n",
      "are pending genotyping (BCH/CHOP).\n\n",
      "This panel will be populated when\npending genotypes become available."),
    size = NM$text$body_mm, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(20, 20, 20, 20))

# NB: Fig 3H VAF × growth was a placeholder pending genotype completion;
# never reached the rendered manuscript. Panel computation kept above for
# documentary purposes; disk write disabled until pending genotypes
# unlock the analysis.

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3I: VAF × recurrence (PLACEHOLDER — all 13 cases pending genotype)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Fig 3I: VAF × recurrence (placeholder) ──\n")
recur_vaf <- vaf_data %>% filter(recurrence_num == 1)
cat(sprintf("  Recurrence among VAF patients: %d\n", nrow(recur_vaf)))

fig3i_placeholder <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
    label = paste0("VAF vs Recurrence\n\n",
      "All 13 patients with documented recurrence\n",
      "are pending genotyping (CHOP).\n\n",
      "0 variant-positive patients currently\nhave recurrence data.\n\n",
      "This panel will be populated when\npending genotypes become available."),
    size = NM$text$body_mm, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(20, 20, 20, 20))

# NB: Fig 3I VAF × recurrence was a placeholder pending genotype data;
# never reached the rendered manuscript. Panel computation kept above for
# documentary purposes; disk write disabled until pending genotypes
# unlock the analysis.

# ══════════════════════════════════════════════════════════════════════════════
# Extended Data Fig. 10: VAF × individual high-risk features
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Extended Data Fig. 10: VAF × individual high-risk features ──\n")

hr_features <- c("intranidal_aneurysm_num", "venous_varix_num",
                  "venous_outflow_stenosis_num", "flow_related_aneurysm_num")
hr_labels <- c("Intranidal Aneurysm", "Venous Varix",
               "Venous Outflow Stenosis", "Flow-Related Aneurysm")

hr_vaf_plots <- map2(hr_features, hr_labels, function(feat, lab) {
  vaf_boxplot(vaf_data, feat, lab, paste0("VAF vs ", lab))
})

# Combine non-null plots
hr_vaf_valid <- compact(map(hr_vaf_plots, "plot"))
if (length(hr_vaf_valid) > 0) {
  library(patchwork)
  efig10 <- wrap_plots(hr_vaf_valid, ncol = 2) +
    plot_annotation(
      # title/subtitle in caption
    )
  # Per-feature HR sub-panels (vaf_highrisk_<feature>) were dropped from the
  # supplementary scope; ed_vaf_phenotype now ships only the combined
  # vaf_highrisk panel saved above.
}

# ══════════════════════════════════════════════════════════════════════════════
# Extended Data Fig. 11: VAF × age
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Extended Data Fig. 6A: VAF × age ──\n")
efig6a <- vaf_scatter(vaf_data, "vaf_pct", "age",
  "VAF (%)", "Age at Diagnosis (years)",
  "VAF vs Age at Diagnosis",
  " | Tests clonal expansion over time")
if (!is.null(efig6a)) {
  # vaf_age panel itself isn't shipped; the canonical VAF × age scatter is
  # Fig 2e (vaf_age_scatter, producer 09_F1_km_age.R). We compute efig6a here
  # only for the all_stats hand-off (rho/p/n) downstream.
  all_stats[["VAF_vs_age"]] <- list(rho = efig6a$cor$estimate, p = efig6a$cor$p.value, n = efig6a$n)
}

# ══════════════════════════════════════════════════════════════════════════════
# Extended Data Fig. 6B: VAF by anatomical location (boxplots per lobe)
#
# Hale v2 G3: replace the prior "VAF vs number of lobes involved" panel
# (an ordinal count with no biological interpretation) with a VAF
# distribution *per anatomical region*. Each patient contributes to every
# lobe their AVM involves, so a single patient can appear in multiple
# boxes. Box order is by median VAF within lobe (descending). Only lobes
# with >= 5 involved patients are shown so boxplot summaries are stable.
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Extended Data Fig. 6B: VAF by anatomical location (boxplots per lobe) ──\n")

# Lobe columns + display labels (mirrors 13_F2_genotype_phenotype.R).
loc_cols <- c("loc_frontal", "loc_temporal", "loc_insular", "loc_basal_ganglia",
              "loc_thalamus", "loc_periventricular", "loc_parietal",
              "loc_occipital", "loc_cerebellar", "loc_brainstem",
              "loc_corpus_callosum", "loc_cingulate", "loc_sylvian_fissure")
loc_labels_vec <- c("Frontal", "Temporal", "Insular", "Basal Ganglia",
                    "Thalamus", "Periventricular", "Parietal",
                    "Occipital", "Cerebellar", "Brainstem",
                    "Corpus Callosum", "Cingulate", "Sylvian Fissure")

vaf_loc_long <- vaf_data %>%
  select(vaf_pct, all_of(loc_cols)) %>%
  pivot_longer(cols = all_of(loc_cols), names_to = "location",
               values_to = "involved") %>%
  filter(!is.na(involved), involved == 1, !is.na(vaf_pct)) %>%
  mutate(loc_label = loc_labels_vec[match(location, loc_cols)])

# Keep lobes with >= 5 patients; sort by median VAF desc.
lobe_n <- vaf_loc_long %>%
  count(loc_label, name = "n") %>%
  filter(n >= 5)
lobe_medians <- vaf_loc_long %>%
  filter(loc_label %in% lobe_n$loc_label) %>%
  group_by(loc_label) %>%
  summarise(med = median(vaf_pct), .groups = "drop") %>%
  arrange(desc(med))

vaf_loc_plot <- vaf_loc_long %>%
  filter(loc_label %in% lobe_medians$loc_label) %>%
  left_join(lobe_n, by = "loc_label") %>%
  mutate(loc_label = factor(loc_label, levels = lobe_medians$loc_label),
         x_label   = sprintf("%s\n(n = %d)", loc_label, n),
         x_label   = factor(x_label,
                            levels = unique(x_label[order(loc_label)])))

# Kruskal-Wallis across lobes (treats each patient-lobe row as an
# independent observation; acknowledged caveat: patients appear more
# than once when the AVM is multi-lobar. Reported as descriptive only).
#
# AUDIT 2026-05-12 (F10): this statistic is currently NOT surfaced in
# any manifest fragment, prose, or caption (see edfig_vaf_phenotype
# fragment construction at the bottom of this file — VAF_vs_location
# is in all_stats but not exported). If a future revision wants to
# cite this number, it MUST disclose the non-independence (one
# patient contributes multiple rows when the AVM is multi-lobar) and
# label it explicitly as descriptive. Do NOT cite as an inferential P.
kw_loc <- kruskal.test(vaf_pct ~ loc_label, data = vaf_loc_plot)

p_ef6b <- ggplot(vaf_loc_plot, aes(x = x_label, y = vaf_pct)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey92",
               color = "grey40") +
  geom_jitter(width = 0.18, size = 1.4, alpha = 0.55, color = unname(PAL_VARIANT[["KRAS G12D"]])) +
  # Statistical values (Kruskal-Wallis p, n) live in the figure caption +
  # companion stats CSV — never on the panel.
  labs(x = NULL,
       y = "Variant allele frequency (%)") +
  theme_nature_panel() +
  theme(axis.text.x = element_text(lineheight = 0.95))

all_stats[["VAF_vs_location"]] <- list(
  # 'n' is the number of patient-lobe observations (patients appear once
  # per involved lobe); downstream stats table expects a scalar 'n'.
  n = nrow(vaf_loc_plot),
  p = kw_loc$p.value,
  kw_p = kw_loc$p.value,
  n_patients_unique = length(unique(vaf_data$patient_id)),
  lobes = as.character(lobe_medians$loc_label)
)

# ══════════════════════════════════════════════════════════════════════════════
# Extended Data Fig. 13: VAF × deep venous drainage specifically
# ══════════════════════════════════════════════════════════════════════════════

cat("\n── Extended Data Fig. 6C: VAF × deep drainage (SM's strongest component) ──\n")
efig6c <- vaf_boxplot(vaf_data, "sm_drainage_num", "Deep Venous Drainage",
  "VAF vs Deep Venous Drainage")
if (!is.null(efig6c)) {
  # vaf_drainage_extended panel was dropped from the supplementary scope;
  # vaf_drainage (the canonical pooled boxplot above) is the surviving
  # panel. Stats hand-off retained for downstream stats fragments.
  all_stats[["VAF_vs_deep_drainage"]] <- list(W = efig6c$wilcox$statistic, p = efig6c$wilcox$p.value, n = efig6c$n)
}

# ══════════════════════════════════════════════════════════════════════════════
# HIERARCHICAL VAF ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  HIERARCHICAL VAF-PHENOTYPE ANALYSIS\n")
cat("══════════════════════════════════════════════════════════════\n")

# ── Full genotyped cohort for dose models ────────────────────────────────────
genotyped <- df %>% filter(!is.na(mutation_positive))

# ── A. VAF → All Outcomes (within mut+ only) ────────────────────────────────

vaf_outcomes <- tribble(
  ~var,                       ~label,            ~type,
  "age",                      "Age",             "continuous",
  "ever_ruptured_num",        "Rupture",         "binary",
  "sm_size_num",              "SM Size",         "continuous",
  "sm_drainage_num",          "Drainage",        "binary",
  "sm_eloquence_num",         "Eloquence",       "binary",
  "n_high_risk_num",          "High-Risk",       "continuous",
  "intranidal_aneurysm_num",  "Intranidal",      "binary",
  "compact_nidus_num",        "Compact Nidus",   "binary",
  "prior_seizure_num",        "Seizure",         "binary"
)

cat("\n── A. VAF → Outcomes (all mut+ with VAF) ──\n")
vaf_all_mutp <- map_dfr(seq_len(nrow(vaf_outcomes)), function(i) {
  o <- vaf_outcomes[i, ]
  s <- vaf_data %>% filter(!is.na(.data[[o$var]]))
  if (nrow(s) < 15) return(tibble())
  fam <- if (o$type == "binary") binomial() else gaussian()
  m <- tryCatch(glm(as.formula(paste(o$var, "~ vaf_prop")), data = s, family = fam), error = function(e) NULL)
  if (is.null(m)) return(tibble())
  cs <- summary(m)$coefficients; vr <- grep("vaf", rownames(cs))
  if (length(vr) == 0) return(tibble())
  sp_rho <- NA_real_; sp_p <- NA_real_
  if (o$type == "continuous") {
    sp <- suppressWarnings(cor.test(s$vaf_prop, s[[o$var]], method = "spearman"))
    sp_rho <- sp$estimate; sp_p <- sp$p.value
  }
  tibble(outcome = o$label, beta = cs[vr, 1], p = cs[vr, 4], rho = sp_rho, rho_p = sp_p, n = nrow(s))
})

vaf_all_mutp %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.)))
    cat(sprintf("  %-16s β=%-8.1f rho=%-7s p=%.4f %s (n=%d)\n",
      .$outcome[j], .$beta[j],
      ifelse(is.na(.$rho[j]), "—", sprintf("%.3f", .$rho[j])),
      .$p[j], .$sig[j], .$n[j]))}

# ── B. VAF → Outcomes by specific variant (G12D, G12V) ──────────────────────

cat("\n── B. VAF → Outcomes (KRAS G12D only) ──\n")
g12d_vaf <- vaf_data %>% filter(geno_variant == "KRAS G12D")
cat(sprintf("  G12D with VAF: n=%d\n", nrow(g12d_vaf)))
vaf_g12d <- map_dfr(seq_len(nrow(vaf_outcomes)), function(i) {
  o <- vaf_outcomes[i, ]
  s <- g12d_vaf %>% filter(!is.na(.data[[o$var]]))
  if (nrow(s) < 10) return(tibble())
  fam <- if (o$type == "binary") binomial() else gaussian()
  m <- tryCatch(glm(as.formula(paste(o$var, "~ vaf_prop")), data = s, family = fam), error = function(e) NULL)
  if (is.null(m)) return(tibble())
  cs <- summary(m)$coefficients; vr <- grep("vaf", rownames(cs))
  if (length(vr) == 0) return(tibble())
  sp_rho <- NA_real_
  if (o$type == "continuous") {
    sp <- suppressWarnings(cor.test(s$vaf_prop, s[[o$var]], method = "spearman"))
    sp_rho <- sp$estimate
  }
  tibble(outcome = o$label, beta = cs[vr, 1], p = cs[vr, 4], rho = sp_rho, n = nrow(s))
})

vaf_g12d %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.)))
    cat(sprintf("  %-16s β=%-8.1f rho=%-7s p=%.4f %s (n=%d)\n",
      .$outcome[j], .$beta[j],
      ifelse(is.na(.$rho[j]), "—", sprintf("%.3f", .$rho[j])),
      .$p[j], .$sig[j], .$n[j]))}

cat("\n── B. VAF → Outcomes (KRAS G12V only) ──\n")
g12v_vaf <- vaf_data %>% filter(geno_variant == "KRAS G12V")
cat(sprintf("  G12V with VAF: n=%d\n", nrow(g12v_vaf)))
vaf_g12v <- map_dfr(seq_len(nrow(vaf_outcomes)), function(i) {
  o <- vaf_outcomes[i, ]
  s <- g12v_vaf %>% filter(!is.na(.data[[o$var]]))
  if (nrow(s) < 8) return(tibble())
  fam <- if (o$type == "binary") binomial() else gaussian()
  m <- tryCatch(glm(as.formula(paste(o$var, "~ vaf_prop")), data = s, family = fam), error = function(e) NULL)
  if (is.null(m)) return(tibble())
  cs <- summary(m)$coefficients; vr <- grep("vaf", rownames(cs))
  if (length(vr) == 0) return(tibble())
  sp_rho <- NA_real_
  if (o$type == "continuous") {
    sp <- suppressWarnings(cor.test(s$vaf_prop, s[[o$var]], method = "spearman"))
    sp_rho <- sp$estimate
  }
  tibble(outcome = o$label, beta = cs[vr, 1], p = cs[vr, 4], rho = sp_rho, n = nrow(s))
})

vaf_g12v %>% mutate(sig = ifelse(p < 0.05, "***", ifelse(p < 0.1, ".", ""))) %>%
  {for (j in seq_len(nrow(.)))
    cat(sprintf("  %-16s β=%-8.1f rho=%-7s p=%.4f %s (n=%d)\n",
      .$outcome[j], .$beta[j],
      ifelse(is.na(.$rho[j]), "—", sprintf("%.3f", .$rho[j])),
      .$p[j], .$sig[j], .$n[j]))}

# ── C. VAF Dose Model (all patients: VAF=0 for neg) ─────────────────────────

cat("\n── C. VAF Dose Model: outcome ~ mutation_positive + vaf_dose (all patients) ──\n")
genotyped_dose <- genotyped %>%
  mutate(vaf_dose = ifelse(mutation_positive & !is.na(vaf_prop), vaf_prop,
    ifelse(!mutation_positive, 0, NA)))

vaf_dose_results <- map_dfr(seq_len(nrow(vaf_outcomes)), function(i) {
  o <- vaf_outcomes[i, ]
  s <- genotyped_dose %>% filter(!is.na(.data[[o$var]]) & !is.na(vaf_dose))
  if (nrow(s) < 30) return(tibble())
  fam <- if (o$type == "binary") binomial() else gaussian()
  m <- tryCatch(glm(as.formula(paste(o$var, "~ mutation_positive + vaf_dose")), data = s, family = fam),
    error = function(e) NULL)
  if (is.null(m)) return(tibble())
  cs <- summary(m)$coefficients
  mp_r <- grep("mutation_positiveTRUE", rownames(cs))
  vd_r <- grep("vaf_dose", rownames(cs))
  if (length(mp_r) == 0 | length(vd_r) == 0) return(tibble())
  tibble(
    outcome = o$label,
    mut_beta = cs[mp_r, 1], mut_p = cs[mp_r, 4],
    vaf_beta = cs[vd_r, 1], vaf_p = cs[vd_r, 4],
    n = nrow(s)
  )
})

vaf_dose_results %>% {for (j in seq_len(nrow(.))) {
  mut_sig <- ifelse(.$mut_p[j] < 0.05, "***", ifelse(.$mut_p[j] < 0.1, ".", ""))
  vaf_sig <- ifelse(.$vaf_p[j] < 0.05, "***", ifelse(.$vaf_p[j] < 0.1, ".", ""))
  cat(sprintf("  %-16s mut+: β=%.2f p=%.3f %s | VAF_dose: β=%.1f p=%.4f %s (n=%d)\n",
    .$outcome[j], .$mut_beta[j], .$mut_p[j], mut_sig,
    .$vaf_beta[j], .$vaf_p[j], vaf_sig, .$n[j]))
}}

# ── D. VAF × Location Interactions (within mut+ only) ───────────────────────

cat("\n── D. VAF × Location → Outcomes (within mut+ only) ──\n")
loc_vars <- c("loc_frontal", "loc_temporal", "loc_parietal", "loc_occipital",
  "loc_cerebellar", "loc_basal_ganglia", "loc_insular")
loc_labs <- c("Frontal", "Temporal", "Parietal", "Occipital",
  "Cerebellar", "Basal Ganglia", "Insular")

vaf_loc_int <- list()
for (li in seq_along(loc_vars)) {
  loc <- loc_vars[li]; loc_lab <- loc_labs[li]
  for (i in seq_len(nrow(vaf_outcomes))) {
    o <- vaf_outcomes[i, ]
    s <- vaf_data %>% filter(!is.na(.data[[o$var]]) & !is.na(.data[[loc]]))
    if (nrow(s) < 25 || n_distinct(s[[loc]]) < 2) next
    fam <- if (o$type == "binary") binomial() else gaussian()
    m <- tryCatch(glm(as.formula(paste(o$var, "~ vaf_prop *", loc)), data = s, family = fam),
      error = function(e) NULL)
    if (is.null(m)) next
    cs <- summary(m)$coefficients; ir <- grep(":", rownames(cs))
    if (length(ir) == 0) next
    p_int <- cs[ir, 4]
    if (p_int < 0.15) {
      sig <- ifelse(p_int < 0.05, "***", ifelse(p_int < 0.1, ".", ""))
      cat(sprintf("  VAF × %-12s → %-12s β=%.1f p=%.4f %s (n=%d)\n",
        loc_lab, o$label, cs[ir, 1], p_int, sig, nrow(s)))
      vaf_loc_int[[paste(loc, o$var)]] <- tibble(
        location = loc_lab, outcome = o$label,
        beta = cs[ir, 1], p = p_int, n = nrow(s))
    }
  }
}

vaf_loc_int_df <- bind_rows(vaf_loc_int)

# ══════════════════════════════════════════════════════════════════════════════
# Supplementary Table 4: All VAF correlation stats
# ══════════════════════════════════════════════════════════════════════════════

# Compile summary (exclude regression table entries that have different format)
simple_stats <- all_stats[!names(all_stats) %in% "rupture_vaf_adjusted"]
etable4 <- tibble(
  Analysis = names(simple_stats),
  n = map_int(simple_stats, "n"),
  Statistic = map_chr(simple_stats, function(s) {
    if (!is.null(s$rho)) sprintf("rho=%.3f", s$rho)
    else if (!is.null(s$W)) sprintf("W=%.0f", s$W)
    else "—"
  }),
  p_value = map_dbl(simple_stats, "p")
)
etable4$p_BH_FDR <- p.adjust(etable4$p_value, method = "BH")

cat("\n── Supplementary Table 4: VAF correlation summary ──\n")
print(etable4)

all_stats[["vaf_all_mutp"]] <- vaf_all_mutp
all_stats[["vaf_g12d"]] <- vaf_g12d
all_stats[["vaf_g12v"]] <- vaf_g12v
all_stats[["vaf_dose_model"]] <- vaf_dose_results
all_stats[["vaf_location_interactions"]] <- vaf_loc_int_df

saveRDS(all_stats, file.path(output_dir, "stats", "fig3_stats.rds"))

library(writexl)
etable_dir <- file.path(output_dir, "SupplementaryTables")
dir.create(etable_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(etable_dir, "_supporting"), recursive = TRUE, showWarnings = FALSE)

write_xlsx(etable4, file.path(etable_dir, "_supporting", "SupplementaryTable_vaf_correlations.xlsx"))

# ─────────────────────────────────────────────────────────────────────────────
# Manifest fragment — edfig_vaf_phenotype
# -----------------------------------------------------------------------------
# Feeds Extended Data Figs 05 (VAF vs individual high-risk features) and 06
# (VAF × age / VAF by location / VAF × drainage). Flat-scalar key space
# lets Results §2 cite the aggregate Spearman ρ range + min p without
# reaching into the nested list structure of `all_stats`.
# ─────────────────────────────────────────────────────────────────────────────
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

.stat_or <- function(lst, key, default = NA_real_) {
  v <- lst[[key]]
  if (is.null(v)) return(default)
  unname(v)
}

.maybe_num <- function(x) if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x)

# Consolidated correlation table (one row per outcome, filled with whatever
# combination of rho/p/W/p was computed — Spearman for ordinal, Wilcoxon
# for binary). NA rows left intact so downstream ED captions can detect
# missing measurements.
corr_rows <- list(
  c("sm_total",   "spearman", .maybe_num(all_stats$VAF_vs_SM_total$rho),
    .maybe_num(all_stats$VAF_vs_SM_total$p),   .maybe_num(all_stats$VAF_vs_SM_total$n)),
  c("sm_size",    "spearman", .maybe_num(all_stats$VAF_vs_SM_size$rho),
    .maybe_num(all_stats$VAF_vs_SM_size$p),    .maybe_num(all_stats$VAF_vs_SM_size$n)),
  c("high_risk",  "spearman", .maybe_num(all_stats$VAF_vs_high_risk$rho),
    .maybe_num(all_stats$VAF_vs_high_risk$p),  .maybe_num(all_stats$VAF_vs_high_risk$n)),
  c("age",        "spearman", .maybe_num(all_stats$VAF_vs_age$rho),
    .maybe_num(all_stats$VAF_vs_age$p),        .maybe_num(all_stats$VAF_vs_age$n)),
  c("drainage",   "wilcoxon", .maybe_num(all_stats$VAF_vs_drainage$W),
    .maybe_num(all_stats$VAF_vs_drainage$p),   .maybe_num(all_stats$VAF_vs_drainage$n)),
  c("eloquence",  "wilcoxon", .maybe_num(all_stats$VAF_vs_eloquence$W),
    .maybe_num(all_stats$VAF_vs_eloquence$p),  .maybe_num(all_stats$VAF_vs_eloquence$n)),
  c("rupture",    "wilcoxon", .maybe_num(all_stats$VAF_vs_rupture$W),
    .maybe_num(all_stats$VAF_vs_rupture$p),    .maybe_num(all_stats$VAF_vs_rupture$n))
)
corr_df <- do.call(rbind, lapply(corr_rows, function(r) {
  data.frame(outcome = r[1], test = r[2],
             estimate = as.numeric(r[3]),
             p        = as.numeric(r[4]),
             n        = as.integer(as.numeric(r[5])),
             stringsAsFactors = FALSE)
}))

rho_vals <- corr_df$estimate[corr_df$test == "spearman"]
p_vals   <- corr_df$p

# Null-phenotype subset: excludes VAF × age (the dose-response, which is the
# *non-null* correlation reported in Fig 3E / §3). Prose in §2 cites the
# range of the FIVE null-phenotype correlations only (sm_total, sm_size,
# high_risk rhos plus drainage / eloquence / rupture Wilcoxons).
null_rho_vals <- corr_df$estimate[corr_df$test == "spearman" &
                                   corr_df$outcome != "age"]
null_p_vals   <- corr_df$p[corr_df$outcome != "age"]

rupt_adj <- all_stats$rupture_vaf_adjusted
rupt_adj_vaf_row <- if (!is.null(rupt_adj) && "term" %in% names(rupt_adj)) {
  rupt_adj[rupt_adj$term == "vaf_pct", , drop = FALSE]
} else NULL

edfig_vaf_phenotype_fragment <- list(
  rho_sm_total      = .maybe_num(all_stats$VAF_vs_SM_total$rho),
  p_sm_total        = .maybe_num(all_stats$VAF_vs_SM_total$p),
  n_sm_total        = .maybe_num(all_stats$VAF_vs_SM_total$n),
  rho_sm_size       = .maybe_num(all_stats$VAF_vs_SM_size$rho),
  p_sm_size         = .maybe_num(all_stats$VAF_vs_SM_size$p),
  n_sm_size         = .maybe_num(all_stats$VAF_vs_SM_size$n),
  rho_high_risk     = .maybe_num(all_stats$VAF_vs_high_risk$rho),
  p_high_risk       = .maybe_num(all_stats$VAF_vs_high_risk$p),
  n_high_risk       = .maybe_num(all_stats$VAF_vs_high_risk$n),
  rho_age           = .maybe_num(all_stats$VAF_vs_age$rho),
  p_age             = .maybe_num(all_stats$VAF_vs_age$p),
  n_age             = .maybe_num(all_stats$VAF_vs_age$n),
  wilcox_drainage   = .maybe_num(all_stats$VAF_vs_drainage$W),
  p_drainage        = .maybe_num(all_stats$VAF_vs_drainage$p),
  n_drainage        = .maybe_num(all_stats$VAF_vs_drainage$n),
  wilcox_eloquence  = .maybe_num(all_stats$VAF_vs_eloquence$W),
  p_eloquence       = .maybe_num(all_stats$VAF_vs_eloquence$p),
  n_eloquence       = .maybe_num(all_stats$VAF_vs_eloquence$n),
  wilcox_rupture    = .maybe_num(all_stats$VAF_vs_rupture$W),
  p_rupture         = .maybe_num(all_stats$VAF_vs_rupture$p),
  n_rupture         = .maybe_num(all_stats$VAF_vs_rupture$n),
  rho_min           = if (length(null_rho_vals)) min(null_rho_vals, na.rm = TRUE) else NA_real_,
  rho_max           = if (length(null_rho_vals)) max(null_rho_vals, na.rm = TRUE) else NA_real_,
  p_min             = if (length(null_p_vals))   min(null_p_vals,   na.rm = TRUE) else NA_real_,
  p_max             = if (length(null_p_vals))   max(null_p_vals,   na.rm = TRUE) else NA_real_,
  # BH-FDR across the 6-test null family (audit B4: surface to prose so
  # §4 ¶11 can cite "all q ≥ X" consistent with the adjacent paragraph's
  # BH-adjusted P framing).
  q_fdr_min         = if (length(null_p_vals))
                        min(p.adjust(null_p_vals, method = "BH"), na.rm = TRUE)
                      else NA_real_,
  q_fdr_max         = if (length(null_p_vals))
                        max(p.adjust(null_p_vals, method = "BH"), na.rm = TRUE)
                      else NA_real_,
  rupt_vaf_or       = if (!is.null(rupt_adj_vaf_row) && nrow(rupt_adj_vaf_row) == 1L)
                        rupt_adj_vaf_row$estimate else NA_real_,
  rupt_vaf_ci_lo    = if (!is.null(rupt_adj_vaf_row) && nrow(rupt_adj_vaf_row) == 1L)
                        rupt_adj_vaf_row$conf.low else NA_real_,
  rupt_vaf_ci_hi    = if (!is.null(rupt_adj_vaf_row) && nrow(rupt_adj_vaf_row) == 1L)
                        rupt_adj_vaf_row$conf.high else NA_real_,
  rupt_vaf_p        = if (!is.null(rupt_adj_vaf_row) && nrow(rupt_adj_vaf_row) == 1L)
                        rupt_adj_vaf_row$p.value else NA_real_,
  correlation_table      = corr_df,
  rupture_adjusted_coefs = if (!is.null(rupt_adj)) as.data.frame(rupt_adj) else NA
)
write_stats_section(section = "edfig_vaf_phenotype",
                    stats   = edfig_vaf_phenotype_fragment)

cat("\n══ 14_ED7_vaf_phenotype.R complete ══\n")
cat("Stats: VAF all mut+, G12D, G12V, dose model, location interactions\n")
