# =============================================================================
# formatters.R — unified formatters for manuscript stats
# -----------------------------------------------------------------------------
# Every number that appears in the Quarto .qmd files should pass through a
# formatter here. Guarantees: consistent Unicode glyphs (−, ×, –, ρ, β),
# consistent sig-fig rules, consistent P-value rendering, consistent OR/CI
# format. All formatters return plain character strings suitable for inline
# use in Quarto markdown.
#
# Contract: these helpers are PURE — no side effects, no I/O. They must be
# deterministic given their inputs so renders are reproducible.
#
# Typical usage in a .qmd:
#   `r fmt_p(stats$fig4a$logrank_p)`
#   `r fmt_or_ci(stats$fig5$firth_or, stats$fig5$firth_ci)`
#   `r fmt_median_ci(stats$fig4a$median_g12d, stats$fig4a$ci_g12d)`
# =============================================================================

# ---- internal helpers -------------------------------------------------------

.unicode_minus <- function(x) gsub("-", "\u2212", x, fixed = TRUE)

.sigfig <- function(x, digits = 3) {
  if (is.na(x)) return(NA_character_)
  s <- format(signif(x, digits), scientific = FALSE, trim = TRUE,
              drop0trailing = TRUE)
  sub("\\.$", "", s)
}

# ---- P-values ---------------------------------------------------------------

#' Format a P-value for manuscript prose.
#'
#' @param p numeric P-value in [0, 1].
#' @param threshold below this, switch to scientific notation (default 1e-4).
#' @param digits significant figures for mantissa (default 2).
#' @return character string, e.g. "0.044", "1.8 \u00D7 10^-6^", "< 0.001".
#'
#' Scientific-notation exponents are emitted as Pandoc superscript syntax
#' (`^N^`). Pandoc converts these to true Word `<w:vertAlign w:val=
#' "superscript"/>` runs on docx output, rather than the Unicode-glyph
#' approximations (\u00B2, \u207B) which are plain text and don't
#' scale, copy-paste, or accessibility-read correctly. The leading sign
#' uses a U+2212 minus for proper typesetting; spaces are not allowed
#' inside `^...^` so the mantissa stays outside the superscript span.
# Format a p-value with the appropriate operator already attached:
#   fmt_p_op(0.0008)  ->  "< 0.001"
#   fmt_p_op(0.043)   ->  "= 0.043"
#   fmt_p_op(2e-7)    ->  "= 2 × 10^-7^"
# Use as: `*P* `r fmt_p_op(p)`` — avoids the "P = < 0.001" double-symbol bug.
fmt_p_op <- function(p, threshold = 1e-4, digits = 2) {
  s <- fmt_p(p, threshold = threshold, digits = digits)
  if (length(s) == 0L) return("NA")
  if (grepl("^<", s) || identical(s, "NA")) return(s)
  paste0("= ", s)
}

fmt_p <- function(p, threshold = 1e-4, digits = 2) {
  if (length(p) == 0L || is.na(p)) return("NA")
  if (p < 1e-300) return("< 10^\u2212300^")
  if (p >= threshold) {
    # conventional decimal
    rounded <- round(p, 3)
    if (rounded < 0.001) return("< 0.001")
    return(formatC(rounded, format = "f", digits = 3))
  }
  # scientific: "mantissa \u00D7 10^-exp^"
  exp_part <- floor(log10(p))
  mantissa <- p / 10 ^ exp_part
  mantissa_str <- formatC(signif(mantissa, digits), format = "fg")
  sign_char <- if (exp_part < 0) "\u2212" else ""
  exp_str <- paste0("^", sign_char, abs(exp_part), "^")
  paste0(mantissa_str, " \u00D7 10", exp_str)
}

# ---- Effect sizes + CIs -----------------------------------------------------

#' Format an odds ratio with 95% CI.
#'
#' @param or numeric point estimate.
#' @param ci length-2 numeric vector (lower, upper).
#' @param digits significant figures (default 2).
#' @return e.g. "OR = 8.05 (95% CI, 1.00–1043)".
# Harmonized 2026-05-27: fixed decimal places (was significant figures via
# .sigfig). HR/OR values render uniformly to `digits` dp (2 by convention),
# so a CI that straddles 1.0 no longer mixes precisions (e.g. 0.62\u20131.14, not
# 0.616\u20131.14). All rendered HR/OR values lie in ~0.6\u20136.9, so 2 dp is exact.
fmt_or_ci <- function(or, ci, digits = 2, prefix = "OR") {
  if (length(or) == 0L || length(ci) != 2L) return("NA")
  sprintf("%s = %s (95%% CI, %s\u2013%s)",
          prefix,
          fmt_num(or, digits),
          fmt_num(ci[1], digits),
          fmt_num(ci[2], digits))
}

#' Format a median with 95% CI or IQR.
#'
#' `suffix` (e.g., " years") is placed IMMEDIATELY AFTER the median value
#' and BEFORE the parenthetical CI / IQR \u2014 matches the form authors
#' typically expect ("16 years (95% CI, 14.6\u201317.9)" rather than
#' "16 (95% CI, 14.6\u201317.9) years"). Same convention for fmt_median_iqr.
# Harmonized 2026-05-27: fixed decimal places (was significant figures via
# .sigfig). Survival ages now render uniformly to `digits` decimal places and
# match the figure-legend captions, which use fmt_num(., 1). Callers pass
# digits = 1; this also removes the sigfig-vs-fmt_num boundary discrepancy
# (e.g. a 31.15 CI bound that read 31.2 in prose but 31.1 in the legend).
fmt_median_ci <- function(med, ci, digits = 1, suffix = "") {
  if (length(med) == 0L || length(ci) != 2L) return("NA")
  sprintf("%s%s (95%% CI, %s\u2013%s)",
          fmt_num(med, digits),
          suffix,
          fmt_num(ci[1], digits),
          fmt_num(ci[2], digits))
}

fmt_median_iqr <- function(med, q1, q3, digits = 3, suffix = "") {
  if (length(med) == 0L || length(q1) == 0L || length(q3) == 0L) return("NA")
  sprintf("%s%s (%s\u2013%s)",
          .sigfig(med, digits),
          suffix,
          .sigfig(q1, digits),
          .sigfig(q3, digits))
}

# ---- Counts and percentages -------------------------------------------------

#' Format count (percent) — n (XX.X%).
fmt_n_pct <- function(n, total, digits = 1) {
  stopifnot(length(n) == length(total))
  pct <- 100 * n / total
  sprintf("%d (%s%%)", as.integer(n),
          formatC(pct, format = "f", digits = digits))
}

#' Format an already-computed percent (e.g. 64.86 -> "64.9%").
fmt_pct <- function(x, digits = 1) {
  if (length(x) == 0L || is.na(x)) return("NA")
  sprintf("%s%%", formatC(x, format = "f", digits = digits))
}

#' Format a plain number to fixed decimal places.
fmt_num <- function(x, digits = 2) {
  if (length(x) == 0L || is.na(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

#' Format a number for a "≥ floor" claim: rounds DOWN to `digits` places so the
#' printed bound never exceeds the true value. Use for "all P/Q ≥ x" statements,
#' where ordinary (nearest) rounding can print a floor above the actual minimum
#' (e.g. min q = 0.149 must read "≥ 0.14", not "≥ 0.15", to stay true).
fmt_floor <- function(x, digits = 2) {
  if (length(x) == 0L || is.na(x)) return("NA")
  m <- 10^digits
  formatC(floor(x * m) / m, format = "f", digits = digits)
}

#' Format a number for a "≤ ceiling" claim: rounds UP to `digits` places so the
#' printed bound never falls below the true value. Mirror of fmt_floor() for
#' "all P/Q ≤ x" statements.
fmt_ceiling <- function(x, digits = 2) {
  if (length(x) == 0L || is.na(x)) return("NA")
  m <- 10^digits
  formatC(ceiling(x * m) / m, format = "f", digits = digits)
}

#' Spell out a small integer (1–12) for prose. Falls back to the digit
#' rendering for values outside that range, since longer numerals read
#' more naturally as digits anyway.
.num_word <- function(n) {
  if (length(n) == 0L || is.na(n)) return("NA")
  words <- c("zero", "one", "two", "three", "four", "five",
             "six", "seven", "eight", "nine", "ten", "eleven", "twelve")
  i <- as.integer(round(n))
  if (i >= 0 && i <= 12) words[i + 1L] else format(i)
}

#' Format an IQR as "q1–q3" with a specified number of decimals.
fmt_iqr <- function(q1, q3, digits = 2) {
  sprintf("%s\u2013%s",
          formatC(q1, format = "f", digits = digits),
          formatC(q3, format = "f", digits = digits))
}

# ---- Per-study lookups ------------------------------------------------------

#' Format a single study's mut+ rate as "n_mut/n_total (XX.X%)".
#'
#' @param stats the loaded manifest
#' @param study_name character matching a row of stats$edfig02$mut_rate_by_study$study_clean
#' @return e.g. "95/110 (86.4%)"
fmt_study_rate <- function(stats, study_name) {
  tbl <- stats$edfig02$mut_rate_by_study
  if (is.null(tbl)) stop("fmt_study_rate: stats$edfig02$mut_rate_by_study missing")
  row <- tbl[tbl$study_clean == study_name, , drop = FALSE]
  if (nrow(row) == 0L) {
    stop(sprintf("fmt_study_rate: study '%s' not found", study_name))
  }
  sprintf("%d/%d (%s%%)",
          as.integer(row$n_mut[1]),
          as.integer(row$n[1]),
          formatC(row$mut_rate[1], format = "f", digits = 1))
}

# ---- Correlation + slope ----------------------------------------------------

#' Format a Spearman or Pearson correlation with P-value.
fmt_rho_p <- function(rho, p, digits = 3, symbol = "\u03C1") {
  rho_str <- .unicode_minus(.sigfig(rho, digits))
  sprintf("%s = %s, %s = %s", symbol, rho_str, "*P*", fmt_p(p))
}

#' Format a regression slope with 95% CI.
fmt_slope_ci <- function(beta, ci, digits = 3, unit = "", prefix = "\u03B2") {
  beta_str <- .unicode_minus(.sigfig(beta, digits))
  ci_lo <- .unicode_minus(.sigfig(ci[1], digits))
  ci_hi <- .unicode_minus(.sigfig(ci[2], digits))
  unit_suffix <- if (nzchar(unit)) paste0(" ", unit) else ""
  sprintf("%s = %s%s (95%% CI, %s to %s)",
          prefix, beta_str, unit_suffix, ci_lo, ci_hi)
}

# ---- Convenience: load from manifest ---------------------------------------

#' Safely extract a stat from the manifest with an informative error.
stat <- function(manifest, section, key) {
  if (!section %in% names(manifest)) {
    stop(sprintf("stat(): section '%s' missing from manifest. Schema mismatch?",
                 section), call. = FALSE)
  }
  if (!key %in% names(manifest[[section]])) {
    stop(sprintf("stat(): key '%s$%s' missing from manifest.",
                 section, key), call. = FALSE)
  }
  manifest[[section]][[key]]
}
