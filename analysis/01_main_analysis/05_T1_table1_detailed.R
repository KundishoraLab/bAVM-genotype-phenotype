# 05_T1_table1_detailed.R — Manuscript-format Table 1 (fine-grained per-variant)
# =============================================================================
# Produces the manuscript main-text Table 1 at the per-variant resolution
# specified in the Apr 21 author direction:
#
#   Columns: KRAS G12D | KRAS G12V | KRAS G12C | KRAS G12A | BRAF V600E |
#            Panel-negative tissue-tested | Total cohort
#
#   Rows: n; sex (% male); age at diagnosis, median (IQR); age at surgery,
#         median (IQR); race, n (%); SM total score, median (IQR); SM grade,
#         n (%); eloquent location, n (%); deep venous drainage, n (%);
#         intranidal aneurysm, n (%); venous varix, n (%); flow-related
#         aneurysm, n (%); venous outflow stenosis, n (%); composite
#         high-risk feature count, median (IQR); rupture at time of surgery,
#         n (%); prior rupture history, n (%); seizure history, n (%); bAVM
#         growth prompting surgery, n (%); prior radiation, n (%); radiation
#         response >50%, n (%); pre-operative embolization, n (%); VAF among
#         variant-positive cases, median % (IQR).
#
# Statistical tests:
#   * Continuous variables: Kruskal-Wallis across the six genotype strata
#   * Categorical variables: Fisher's exact test with simulated P (B = 10000)
#   * Benjamini-Hochberg FDR correction applied across all comparisons
#
# Output:
#   results/Table1/Table1_detailed.xlsx     — machine-readable source
#   results/Table1/Table1_detailed.docx     — rendered flextable (insertable)
#   results/Table1/Table1_detailed_ft.rds   — flextable object
#
# The flextable object renders Table 1 in the standalone Word export.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(writexl)
  library(flextable)
  library(officer)
})

source(here("analysis", "helper_scripts", "utils.R"))

# ── 1. Load analysis-ready data ──────────────────────────────────────────────

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))

# Define the manuscript Table 1 strata. `geno_variant == "Negative"` is the
# single source of truth for the panel-negative-tissue-tested column (155
# lesions that went through the assay and were WT at every interrogated
# locus); see 04_T1_table1.R and stats_schema.R for the mapping. We use the
# fine-grained `mutation` column for the five variant-positive columns so
# KRAS G12C (n = 16), G12A (n = 1), and BRAF V600E (n = 4) are broken out
# individually rather than folded into "Other KRAS"/"BRAF" buckets.
df <- df %>%
  mutate(
    table1_stratum = case_when(
      mutation == "KRAS G12D" ~ "KRAS G12D",
      mutation == "KRAS G12V" ~ "KRAS G12V",
      mutation == "KRAS G12C" ~ "KRAS G12C",
      mutation == "KRAS G12A" ~ "KRAS G12A",
      mutation == "BRAF V600E" ~ "BRAF V600E",
      geno_variant == "Negative" ~ "Panel-negative",
      TRUE ~ NA_character_
    ),
    table1_stratum = factor(
      table1_stratum,
      levels = c("KRAS G12D", "KRAS G12V", "KRAS G12C", "KRAS G12A",
                 "BRAF V600E", "Panel-negative")
    )
  )

strata   <- levels(df$table1_stratum)
n_strata <- length(strata)

# Rare variant-positive variants outside the six displayed columns (KRAS
# Q61H, KRAS c.191_196dup, BRAF Q636X) are contained within the Total column
# but do not contribute to any individual stratum. Tissue-pending lesions
# (mutation == NA) are also in Total but not in the strata.
cat(sprintf(
  "Stratum counts: %s | Total cohort n = %d | in-strata n = %d\n",
  paste(sprintf("%s=%d", strata,
                tabulate(df$table1_stratum, nbins = n_strata)),
        collapse = ", "),
  nrow(df),
  sum(!is.na(df$table1_stratum))
))

# ── 2. Formatters ────────────────────────────────────────────────────────────

fmt_n_only  <- function(x) as.character(sum(!is.na(x)))
fmt_med_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  sprintf(paste0("%.", digits, "f (%.", digits, "f\u2013%.", digits, "f)"),
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}
fmt_n_pct_of_group <- function(x, positive_values = c(1, TRUE)) {
  xnm <- x[!is.na(x)]
  if (length(xnm) == 0) return("-")
  n_pos <- sum(xnm %in% positive_values)
  # n/N (%) — N is the non-missing denominator for this cell, so the
  # percentage is reproducible from the cell itself (e.g. BRAF V600E sex =
  # 2/3 (67%), not 2/4). Missingness is therefore visible per cell rather
  # than needing a separate column.
  sprintf("%d/%d (%.0f%%)", n_pos, length(xnm), 100 * n_pos / length(xnm))
}
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Column values for a single continuous variable across all strata + Total.
values_continuous <- function(var, d = df, digits = 1) {
  out <- vapply(strata, function(s) {
    fmt_med_iqr(d[[var]][d$table1_stratum == s & !is.na(d$table1_stratum)],
                digits = digits)
  }, character(1))
  c(out, Total = fmt_med_iqr(d[[var]], digits = digits))
}

# Column values for a single binary / indicator variable.
values_binary <- function(var, d = df) {
  out <- vapply(strata, function(s) {
    fmt_n_pct_of_group(d[[var]][d$table1_stratum == s & !is.na(d$table1_stratum)])
  }, character(1))
  c(out, Total = fmt_n_pct_of_group(d[[var]]))
}

# Column values for a categorical variable returning one row per level.
values_categorical <- function(var, levels_to_show, d = df) {
  map_dfr(levels_to_show, function(lev) {
    per_stratum <- vapply(strata, function(s) {
      x <- d[[var]][d$table1_stratum == s & !is.na(d$table1_stratum)]
      x <- x[!is.na(x)]
      if (length(x) == 0) return("-")
      n_lev <- sum(x == lev)
      sprintf("%d/%d (%.0f%%)", n_lev, length(x), 100 * n_lev / length(x))
    }, character(1))
    x_total <- d[[var]][!is.na(d[[var]])]
    total_cell <- if (length(x_total) == 0) "-" else
      sprintf("%d/%d (%.0f%%)", sum(x_total == lev), length(x_total),
              100 * sum(x_total == lev) / length(x_total))
    tibble(
      Variable = sprintf("  %s", as.character(lev)),
      !!!setNames(as.list(per_stratum), strata),
      Total = total_cell
    )
  })
}

# ── 3. Omnibus tests across the six strata (for BH-FDR) ──────────────────────

kw_p <- function(var, d = df) {
  sub <- d %>% filter(!is.na(.data[[var]]), !is.na(table1_stratum))
  if (length(unique(sub$table1_stratum)) < 2 || nrow(sub) < 3) return(NA_real_)
  tryCatch(kruskal.test(sub[[var]] ~ sub$table1_stratum)$p.value,
           error = function(e) NA_real_)
}

fisher_p <- function(var, d = df) {
  sub <- d %>% filter(!is.na(.data[[var]]), !is.na(table1_stratum))
  if (length(unique(sub$table1_stratum)) < 2) return(NA_real_)
  tbl <- table(sub$table1_stratum, sub[[var]])
  if (ncol(tbl) < 2 || nrow(tbl) < 2) return(NA_real_)
  tryCatch({
    set.seed(MASTER_SEED)   # audit F13: pin Monte-Carlo Fisher P
    fisher.test(tbl, simulate.p.value = TRUE, B = 10000)$p.value
  }, error = function(e) NA_real_)
}

# ── 4. Assemble the table row-by-row ─────────────────────────────────────────

# Each entry in `row_spec` captures (Variable label, values-function, p-value
# function). The ordering is the manuscript row order. For categorical
# variables rendered as one-row-per-level (race, SM grade) we use a single
# omnibus p on the parent row and emit the level rows as sub-rows with
# blank p cells — matching standard Table 1 rendering conventions.

sex_is_male <- as.integer(df$sex_f == "Male")  # 1 = male, 0 = female
# previous_rupture is coded 0/1 as binary already; ruptured_at_surgery_num
# likewise. radiation_shrink_num has levels 0 (no), 1 (yes), 3 (N/A no
# radiation) — only `1` is "response >50%" per data dictionary.

rows <- list()

rows[[length(rows) + 1]] <- list(
  Variable = "n", ValuesFn = function() {
    out <- vapply(strata, function(s) as.character(sum(df$table1_stratum == s, na.rm = TRUE)),
                  character(1))
    c(out, Total = as.character(nrow(df)))
  }, PFn = function() NA_real_, type = "simple"
)

rows[[length(rows) + 1]] <- list(
  Variable = "Sex, % male",
  ValuesFn = function() {
    d <- df %>% mutate(.male = as.integer(sex_f == "Male"))
    values_binary(".male", d)
  },
  PFn = function() {
    d <- df %>% mutate(.male = as.integer(sex_f == "Male"))
    fisher_p(".male", d)
  },
  type = "simple"
)

rows[[length(rows) + 1]] <- list(
  Variable = "Age at diagnosis, median (IQR), years",
  ValuesFn = function() values_continuous("age_dx_numeric"),
  PFn      = function() kw_p("age_dx_numeric"),
  type = "simple"
)

rows[[length(rows) + 1]] <- list(
  Variable = "Age at surgery, median (IQR), years",
  ValuesFn = function() values_continuous("age_surgery_numeric"),
  PFn      = function() kw_p("age_surgery_numeric"),
  type = "simple"
)

# Race — omnibus p on the parent row, per-level cells rendered as sub-rows.
race_levels <- c("White", "Black/African American", "Hispanic/Latino",
                 "Asian", "Other", "Unable to collect")
rows[[length(rows) + 1]] <- list(
  Variable = "Race, n (%)",
  ValuesFn = function() {
    blank <- setNames(rep("", n_strata + 1), c(strata, "Total"))
    list(parent = blank,
         children = values_categorical("race_f", race_levels))
  },
  PFn = function() fisher_p("race_f"),
  type = "categorical"
)

rows[[length(rows) + 1]] <- list(
  Variable = "Spetzler-Martin total score, median (IQR)",
  ValuesFn = function() values_continuous("sm_total_num", digits = 0),
  PFn      = function() kw_p("sm_total_num"),
  type = "simple"
)

# SM grade — also rendered as one row per grade level under an omnibus p.
sm_levels <- c("I", "II", "III", "IV", "V")
rows[[length(rows) + 1]] <- list(
  Variable = "Spetzler-Martin grade, n (%)",
  ValuesFn = function() {
    blank <- setNames(rep("", n_strata + 1), c(strata, "Total"))
    list(parent = blank,
         children = values_categorical("sm_grade", sm_levels))
  },
  PFn = function() fisher_p("sm_grade"),
  type = "categorical"
)

rows[[length(rows) + 1]] <- list(
  Variable = "Eloquent location, n (%)",
  ValuesFn = function() values_binary("sm_eloquence_num"),
  PFn      = function() fisher_p("sm_eloquence_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Deep venous drainage, n (%)",
  ValuesFn = function() values_binary("sm_drainage_num"),
  PFn      = function() fisher_p("sm_drainage_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Intranidal aneurysm, n (%)",
  ValuesFn = function() values_binary("intranidal_aneurysm_num"),
  PFn      = function() fisher_p("intranidal_aneurysm_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Venous varix, n (%)",
  ValuesFn = function() values_binary("venous_varix_num"),
  PFn      = function() fisher_p("venous_varix_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Flow-related aneurysm, n (%)",
  ValuesFn = function() values_binary("flow_related_aneurysm_num"),
  PFn      = function() fisher_p("flow_related_aneurysm_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Venous outflow stenosis, n (%)",
  ValuesFn = function() values_binary("venous_outflow_stenosis_num"),
  PFn      = function() fisher_p("venous_outflow_stenosis_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "High-risk feature count, median (IQR)",
  ValuesFn = function() values_continuous("n_high_risk_num", digits = 0),
  PFn      = function() kw_p("n_high_risk_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Rupture at time of surgery, n (%)",
  ValuesFn = function() values_binary("ruptured_at_surgery_num"),
  PFn      = function() fisher_p("ruptured_at_surgery_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Prior rupture history, n (%)",
  ValuesFn = function() values_binary("previous_rupture"),
  PFn      = function() fisher_p("previous_rupture"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Seizure history, n (%)",
  ValuesFn = function() values_binary("prior_seizure_num"),
  PFn      = function() fisher_p("prior_seizure_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "bAVM growth prompting surgery, n (%)",
  ValuesFn = function() values_binary("growing_num"),
  PFn      = function() fisher_p("growing_num"),
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Prior radiation, n (%)",
  ValuesFn = function() values_binary("prior_radiation_num"),
  PFn      = function() fisher_p("prior_radiation_num"),
  type = "simple"
)

# Radiation response >50% — `radiation_shrink_num` is coded 0/1 among those
# who received radiation (3 = N/A no radiation); restrict denominator to
# irradiated lesions and use `1` as the positive indicator.
rows[[length(rows) + 1]] <- list(
  Variable = "Radiation response >50%, n (%) of irradiated",
  ValuesFn = function() {
    d <- df %>% mutate(
      .resp = ifelse(radiation_shrink_num %in% c(0, 1),
                     as.integer(radiation_shrink_num == 1), NA_integer_)
    )
    values_binary(".resp", d)
  },
  PFn = function() {
    d <- df %>% mutate(
      .resp = ifelse(radiation_shrink_num %in% c(0, 1),
                     as.integer(radiation_shrink_num == 1), NA_integer_)
    )
    fisher_p(".resp", d)
  },
  type = "simple"
)
rows[[length(rows) + 1]] <- list(
  Variable = "Pre-operative embolization, n (%)",
  ValuesFn = function() values_binary("prior_embolization_num"),
  PFn      = function() fisher_p("prior_embolization_num"),
  type = "simple"
)

# VAF is defined only for variant-positive cases, so the panel-negative
# column is "-" by construction. Kruskal-Wallis across the five mut+ strata
# tests whether per-variant VAF distributions differ.
rows[[length(rows) + 1]] <- list(
  Variable = "VAF among variant-positive, median % (IQR)",
  ValuesFn = function() {
    d <- df
    out <- vapply(strata, function(s) {
      x <- d$vaf_pct[d$table1_stratum == s & !is.na(d$table1_stratum)]
      fmt_med_iqr(x, digits = 2)
    }, character(1))
    out["Panel-negative"] <- "-"
    c(out, Total = fmt_med_iqr(
      d$vaf_pct[d$table1_stratum != "Panel-negative"], digits = 2))
  },
  PFn = function() {
    d <- df %>% filter(table1_stratum != "Panel-negative",
                       !is.na(table1_stratum), !is.na(vaf_pct))
    if (length(unique(d$table1_stratum)) < 2 || nrow(d) < 3) return(NA_real_)
    tryCatch(kruskal.test(d$vaf_pct ~ droplevels(d$table1_stratum))$p.value,
             error = function(e) NA_real_)
  },
  type = "simple"
)

# ── 5. Materialise rows into a tibble with raw + BH-FDR-adjusted p ───────────

# Flatten: each entry becomes one or more rows, with omnibus p attached to
# the parent row (sub-rows have blank p cells).

flat_rows <- list()
raw_p     <- c()  # one p per parent row, in order

for (r in rows) {
  if (r$type == "simple") {
    vals <- r$ValuesFn()
    flat_rows[[length(flat_rows) + 1]] <- tibble(
      Variable = r$Variable,
      !!!setNames(as.list(vals[strata]), strata),
      Total    = unname(vals["Total"]),
      .row_has_p = TRUE
    )
    raw_p <- c(raw_p, r$PFn())
  } else if (r$type == "categorical") {
    nested <- r$ValuesFn()
    flat_rows[[length(flat_rows) + 1]] <- tibble(
      Variable = r$Variable,
      !!!setNames(as.list(nested$parent[strata]), strata),
      Total    = unname(nested$parent["Total"]),
      .row_has_p = TRUE
    )
    raw_p <- c(raw_p, r$PFn())
    children <- nested$children %>% mutate(.row_has_p = FALSE)
    flat_rows[[length(flat_rows) + 1]] <- children
  }
}

table1 <- bind_rows(flat_rows)

# Attach per-row p values (raw + BH across the family of omnibus tests).
parent_ix <- which(table1$.row_has_p)
stopifnot(length(parent_ix) == length(raw_p))
testable  <- !is.na(raw_p)
adj_p     <- rep(NA_real_, length(raw_p))
if (any(testable)) adj_p[testable] <- p.adjust(raw_p[testable], method = "BH")

p_col     <- rep("", nrow(table1))
p_adj_col <- rep("", nrow(table1))
p_col[parent_ix]     <- vapply(raw_p, fmt_p, character(1))
p_adj_col[parent_ix] <- vapply(adj_p, fmt_p, character(1))

table1$p_raw    <- p_col
table1$p_BH_FDR <- p_adj_col
table1$.row_has_p <- NULL

# ── 6. Header with n-per-column injected ─────────────────────────────────────

n_by_col <- c(
  vapply(strata, function(s) sum(df$table1_stratum == s, na.rm = TRUE),
         integer(1)),
  Total = nrow(df)
)

col_headers <- c(
  Variable = "Variable",
  setNames(
    sprintf("%s\n(n = %d)", c(strata, "Total cohort"), n_by_col),
    c(strata, "Total")
  ),
  p_raw    = "Raw P",
  p_BH_FDR = "BH-FDR P"
)

# ── 7. Flextable rendering ───────────────────────────────────────────────────

ft <- flextable(table1) %>%
  set_header_labels(values = as.list(col_headers)) %>%
  fontsize(size = 8, part = "all") %>%
  padding(padding.top = 2, padding.bottom = 2, part = "all") %>%
  align(align = "left", part = "all") %>%
  align(j = 2:(n_strata + 2), align = "center", part = "all") %>%
  align(j = c("p_raw", "p_BH_FDR"), align = "right", part = "all") %>%
  bold(part = "header") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1)) %>%
  width(j = 1, width = 1.70) %>%
  width(j = 2:(n_strata + 2), width = 0.60) %>%
  width(j = (n_strata + 3):(n_strata + 4), width = 0.475) %>%
  set_table_properties(layout = "fixed")

# ── 8. Output: Excel (data), Word (standalone), RDS ───────

out_dir <- here("results", "Table1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Excel: swap column names to the human-readable headers for machine use too.
excel_table <- table1
names(excel_table) <- c(col_headers[names(table1)])
write_xlsx(
  list("Table1_detailed" = excel_table),
  path = file.path(out_dir, "Table1_detailed.xlsx")
)

# Standalone Word doc (useful to paste into coauthor drafts).
tbl_doc <- read_docx() %>%
  body_add_par("Table 1. Clinical and demographic characteristics stratified by somatic mutation.",
               style = "heading 2") %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par(paste0(
    "Continuous variables compared by Kruskal-Wallis test across the six ",
    "genotype strata; categorical variables compared by Fisher's exact test ",
    "(simulated P, B = 10000). Benjamini-Hochberg false-discovery-rate ",
    "correction applied across all omnibus comparisons. Categorical cells ",
    "are shown as n/N (%), where N is the number of non-missing observations ",
    "in that stratum; the percentage is therefore of non-missing observations."
  ), style = "Normal")
print(tbl_doc, target = file.path(out_dir, "Table1_detailed.docx"))

# RDS.
saveRDS(
  list(tibble = table1, flextable = ft, headers = col_headers,
       n_by_col = n_by_col),
  file.path(out_dir, "Table1_detailed_ft.rds")
)

cat(sprintf(
  "\n\u2713 Table 1 detailed built: %d rows across %d strata + Total.\n",
  nrow(table1), n_strata))
cat(sprintf("  \u2192 %s\n", file.path(out_dir, "Table1_detailed.xlsx")))
cat(sprintf("  \u2192 %s\n", file.path(out_dir, "Table1_detailed.docx")))
cat(sprintf("  \u2192 %s\n", file.path(out_dir, "Table1_detailed_ft.rds")))
