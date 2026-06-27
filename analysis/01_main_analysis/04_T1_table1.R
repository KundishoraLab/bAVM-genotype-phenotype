# 04_T1_table1.R — Cohort characterization: Table 1
#
# Input:  data/processed/bAVM_analysis_ready.rds
# Output: Table 1 (Tier 1: binary), Table 1 alt (Tier 2: per-variant),
#         Supplementary Table 1 (extended), Supplementary Table 2 (full cohort including pending)
#
# All continuous variables: median (IQR), Kruskal-Wallis test
# All categorical variables: n (%), Fisher's exact test
# BH-FDR correction applied across all comparisons within each table
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(writexl)

source(here("analysis", "helper_scripts", "utils.R"))

# ── 1. Load analysis-ready data ─────────────────────────────────────────────

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))
genotyped <- df %>% filter(!is.na(mutation_positive))

cat(sprintf("Analysis-ready dataset: %d total, %d genotyped\n", nrow(df), nrow(genotyped)))

# ── 2. Helper functions ─────────────────────────────────────────────────────

# Format median (IQR) for a numeric vector
fmt_median_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  sprintf("%.1f (%.1f–%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

# Format n (%) for a binary (0/1) or logical variable
fmt_n_pct <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  n_pos <- sum(x == 1 | x == TRUE)
  sprintf("%d (%.1f%%)", n_pos, 100 * n_pos / length(x))
}

# Format n (%) for categorical levels
fmt_cat_levels <- function(x, levels_to_show = NULL) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(setNames("-", "all"))
  if (is.null(levels_to_show)) levels_to_show <- sort(unique(x))
  sapply(levels_to_show, function(lev) {
    n_lev <- sum(x == lev)
    sprintf("%d (%.1f%%)", n_lev, 100 * n_lev / length(x))
  })
}

# Get n available for a variable in a subset
n_avail <- function(x) sum(!is.na(x))

# ── 3. Build Table 1 function ───────────────────────────────────────────────

build_table1 <- function(data, group_var, group_levels, table_label = "Table 1") {

  # Initialize results
  rows <- list()
  p_values <- c()

  # Helper: add a row to the table
  add_row <- function(label, variable, type = "binary", levels_to_show = NULL, suffix = "") {
    vals <- list()
    all_vals_by_group <- list()

    for (g in group_levels) {
      sub <- data %>% filter(.data[[group_var]] == g)
      x <- sub[[variable]]
      n <- n_avail(x)
      if (type == "continuous") {
        vals[[g]] <- paste0(fmt_median_iqr(x), "\n(n=", n, ")")
      } else if (type == "binary") {
        vals[[g]] <- paste0(fmt_n_pct(x), "\n(n=", n, ")")
      } else if (type == "categorical") {
        level_strs <- fmt_cat_levels(x, levels_to_show)
        vals[[g]] <- paste(paste0("  ", names(level_strs), ": ", level_strs), collapse = "\n")
        vals[[g]] <- paste0(vals[[g]], "\n(n=", n, ")")
      }
      all_vals_by_group[[g]] <- x[!is.na(x)]
    }

    # Total column
    x_total <- data[[variable]]
    n_total <- n_avail(x_total)
    if (type == "continuous") {
      vals[["Total"]] <- paste0(fmt_median_iqr(x_total), "\n(n=", n_total, ")")
    } else if (type == "binary") {
      vals[["Total"]] <- paste0(fmt_n_pct(x_total), "\n(n=", n_total, ")")
    } else if (type == "categorical") {
      level_strs <- fmt_cat_levels(x_total, levels_to_show)
      vals[["Total"]] <- paste(paste0("  ", names(level_strs), ": ", level_strs), collapse = "\n")
      vals[["Total"]] <- paste0(vals[["Total"]], "\n(n=", n_total, ")")
    }

    # Statistical test
    test_data <- data %>%
      filter(!is.na(.data[[variable]]) & !is.na(.data[[group_var]])) %>%
      select(group = all_of(group_var), value = all_of(variable))

    p_val <- NA_real_
    test_name <- ""
    tryCatch({
      if (type == "continuous") {
        if (length(unique(test_data$group)) >= 2 && nrow(test_data) >= 3) {
          kw <- kruskal.test(value ~ group, data = test_data)
          p_val <- kw$p.value
          test_name <- "Kruskal-Wallis"
        }
      } else {
        if (length(unique(test_data$group)) >= 2 && nrow(test_data) >= 3) {
          tbl <- table(test_data$group, test_data$value)
          if (ncol(tbl) >= 2 && nrow(tbl) >= 2) {
            # 2x2: use exact Fisher (deterministic, no simulation noise).
            # Larger tables: simulated P with a fixed seed so the
            # reported P is bit-reproducible across runs (audit D6).
            if (nrow(tbl) == 2L && ncol(tbl) == 2L) {
              ft <- fisher.test(tbl)
              test_name <- "Fisher's exact"
            } else {
              set.seed(MASTER_SEED)   # audit F13: canonical seed from utils.R
              ft <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)
              test_name <- "Fisher's exact (simulated)"
            }
            p_val <- ft$p.value
          }
        }
      }
    }, error = function(e) {
      p_val <<- NA_real_
      test_name <<- paste0("Error: ", e$message)
    })

    list(
      label = paste0(label, suffix),
      values = vals,
      p_raw = p_val,
      test = test_name
    )
  }

  # ── Build rows ──

  # N per group
  n_row <- list(label = "n", values = list(), p_raw = NA_real_, test = "")
  for (g in group_levels) {
    n_row$values[[g]] <- as.character(sum(data[[group_var]] == g, na.rm = TRUE))
  }
  n_row$values[["Total"]] <- as.character(nrow(data))
  rows[[length(rows) + 1]] <- n_row

  # Demographics
  rows[[length(rows) + 1]] <- add_row("Sex (% male)", "sex", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Age at presentation, years", "age", type = "continuous")
  rows[[length(rows) + 1]] <- add_row("Age at surgery, years", "age_surgery_numeric", type = "continuous")
  rows[[length(rows) + 1]] <- add_row("Pediatric (<18)", "age_group", type = "binary")

  # SM components
  rows[[length(rows) + 1]] <- add_row("SM size score", "sm_size_num", type = "continuous")
  rows[[length(rows) + 1]] <- add_row("SM drainage (deep)", "sm_drainage_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("SM eloquence", "sm_eloquence_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("SM total score", "sm_total_num", type = "continuous")
  rows[[length(rows) + 1]] <- add_row("SM grade (composite)", "sm_grade", type = "categorical",
    levels_to_show = c("I", "II", "III", "IV", "V"))

  # Angioarchitecture
  rows[[length(rows) + 1]] <- add_row("Rupture (ever)", "ever_ruptured_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("High-risk feature count", "n_high_risk_num", type = "continuous")
  rows[[length(rows) + 1]] <- add_row("Intranidal aneurysm", "intranidal_aneurysm_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Venous varix", "venous_varix_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Venous outflow stenosis", "venous_outflow_stenosis_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Flow-related aneurysm", "flow_related_aneurysm_num", type = "binary")

  # Clinical history
  rows[[length(rows) + 1]] <- add_row("Seizure history", "prior_seizure_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Prior radiation", "prior_radiation_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Prior embolization", "prior_embolization_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Growth", "growing_num", type = "binary")
  rows[[length(rows) + 1]] <- add_row("Compact nidus", "compact_nidus_num", type = "binary")

  # VAF (variant-positive only — report overall for the mut+ column)
  rows[[length(rows) + 1]] <- add_row("VAF (%, mut+ only)", "vaf_prop", type = "continuous")

  # ── Compile into data frame ──
  p_raw_vec <- sapply(rows, function(r) r$p_raw)
  # BH-FDR across testable rows
  testable <- !is.na(p_raw_vec)
  p_adj_vec <- rep(NA_real_, length(p_raw_vec))
  if (sum(testable) > 0) {
    p_adj_vec[testable] <- p.adjust(p_raw_vec[testable], method = "BH")
  }

  # Format p-values
  fmt_p <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("<0.001")
    sprintf("%.3f", p)
  }

  # Build output data frame
  out <- tibble(Variable = sapply(rows, function(r) r$label))
  for (g in c(group_levels, "Total")) {
    out[[g]] <- sapply(rows, function(r) r$values[[g]] %||% "")
  }
  out$p_raw <- sapply(p_raw_vec, fmt_p)
  out$p_BH_FDR <- sapply(p_adj_vec, fmt_p)
  out$test <- sapply(rows, function(r) r$test)

  cat(sprintf("\n── %s built: %d rows, %d groups ──\n", table_label, nrow(out), length(group_levels)))

  return(out)
}

# ── 4. Table 1 — Tier 1 (Binary: Mut+ vs Neg) ──────────────────────────────

cat("\n════════════════════════════════════════\n")
cat("Building Table 1 — Tier 1 (Binary)\n")
cat("════════════════════════════════════════\n")

table1_binary <- build_table1(
  data = genotyped,
  group_var = "geno_binary",
  group_levels = c("Variant-positive", "Panel-negative"),
  table_label = "Table 1 (Binary)"
)

# ── 5. Table 1 alt — Tier 2 (Per-variant) ───────────────────────────────────

cat("\n════════════════════════════════════════\n")
cat("Building Table 1 alt — Tier 2 (Per-variant)\n")
cat("════════════════════════════════════════\n")

# Relabel raw "Negative" → canonical "Panel-negative" before build_table1
# so the Excel sheet column header matches the manuscript convention
# (the Tier-1 binary call above already uses "Panel-negative" — this
# brings the Tier-2 per-variant table into alignment).
table1_variant <- build_table1(
  data = genotyped %>%
    mutate(geno_variant = relabel_geno_factor(geno_variant)),
  group_var = "geno_variant",
  group_levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF",
                   "Panel-negative"),
  table_label = "Table 1 (Per-variant)"
)

# ── 6. Supplementary Table 2 — Full cohort including pending ─────────────────────────────

cat("\n════════════════════════════════════════\n")
cat("Building Supplementary Table 2 — Full cohort (including pending)\n")
cat("════════════════════════════════════════\n")

etable2_full <- build_table1(
  data = df,
  group_var = "geno_status",
  group_levels = c("Variant-positive", "Panel-negative", "Pending"),
  table_label = "Supplementary Table 2 (Full cohort)"
)

# ── 7. Export ────────────────────────────────────────────────────────────────

output_dir <- here("results")

# Clean up newlines for Excel export (replace \n with " | ")
clean_for_excel <- function(tbl) {
  tbl %>% mutate(across(where(is.character), ~ str_replace_all(., "\n", " | ")))
}

etable_dir <- file.path(output_dir, "SupplementaryTables")
dir.create(etable_dir, recursive = TRUE, showWarnings = FALSE)
write_xlsx(
  list(
    "Table1_Binary" = clean_for_excel(table1_binary),
    "Table1_PerVariant" = clean_for_excel(table1_variant),
    "SupplementaryTable2_FullCohort" = clean_for_excel(etable2_full)
  ),
  path = file.path(etable_dir, "Table1_all_versions.xlsx")
)

# Also save as RDS for downstream use
saveRDS(
  list(
    table1_binary = table1_binary,
    table1_variant = table1_variant,
    etable2_full = etable2_full
  ),
  file.path(output_dir, "stats", "table1_results.rds")
)

cat("\n── Tables exported to output/ ──\n")
cat("  Table1_all_versions.xlsx (3 sheets)\n")
cat("  table1_results.rds\n")

# ── 8. Print summary to console ─────────────────────────────────────────────

cat("\n══════════════════════════════════════════════\n")
cat("Table 1 Binary — Significant results (FDR < 0.05)\n")
cat("══════════════════════════════════════════════\n")
sig_binary <- table1_binary %>%
  filter(p_BH_FDR != "" & p_BH_FDR != "<0.001") %>%
  filter(as.numeric(p_BH_FDR) < 0.05)
sig_binary_strict <- table1_binary %>%
  filter(p_BH_FDR == "<0.001")
bind_rows(sig_binary, sig_binary_strict) %>%
  select(Variable, p_raw, p_BH_FDR, test) %>%
  print(n = 50)

cat("\nTable 1 Per-variant — Significant results (FDR < 0.05)\n")
sig_variant <- table1_variant %>%
  filter(p_BH_FDR != "" & p_BH_FDR != "<0.001") %>%
  filter(as.numeric(p_BH_FDR) < 0.05)
sig_variant_strict <- table1_variant %>%
  filter(p_BH_FDR == "<0.001")
bind_rows(sig_variant, sig_variant_strict) %>%
  select(Variable, p_raw, p_BH_FDR, test) %>%
  print(n = 50)

cat("\n══ 04_T1_table1.R complete ══\n")
