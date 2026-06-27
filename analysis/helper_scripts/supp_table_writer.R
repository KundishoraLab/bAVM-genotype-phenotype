# =============================================================================
# supp_table_writer.R — Validated, Nature-Medicine-styled writer for the
# manuscript's supplementary tables.
# -----------------------------------------------------------------------------
# Why this exists
# -----------------------------------------------------------------------------
# Two failure modes have shipped from the existing writexl::write_xlsx() calls:
#   (1) silently-NA rows that look like "model didn't fit" but were really
#       an extractor bug (see ST04 logistf row, 2026-06-14 audit);
#   (2) raw R-style headers ("n_mut", "interaction_p_fdr") leaking into the
#       reviewer-facing supplement.
#
# write_supp_table() prevents both:
#   - Hard-errors on NA / NaN / Inf / empty strings / wrong type. Caller must
#     convert intentional "not estimable" cells to a sentinel string column
#     (e.g. "—" or "n<30") and declare that column as type = "character".
#   - Renders Nature-Medicine-style headers: bold + thin underline, italic on
#     stat-symbol header cells (N, P, β, χ², HR, OR, etc.). Italic is
#     opt-in per column via italic_header = TRUE.
#
# API
# -----------------------------------------------------------------------------
# write_supp_table(
#   data    = <data.frame>,
#   path    = "results/SupplementaryTables/SuppTableNN_topic.xlsx",
#   sheet   = "Mutation rate by study",
#   columns = list(
#     col(name = "study_clean", label = "Study"),
#     col_int("n",     label = "N",          italic = TRUE),
#     col_int("n_mut", label = "N (mutant)", italic = TRUE),
#     col_num("mut_rate", label = "Mutation rate (%)", digits = 1)
#   )
# )
#
# For multi-sheet workbooks, call write_supp_table_workbook() with a named list
# of (data, columns) pairs.
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
})

# ── Column-spec constructors ──────────────────────────────────────────────────
# Keep the producer call sites short. Each returns a named list the writer
# reads. type is one of "character", "integer", "numeric"; digits applies only
# to numeric (Excel display format, NOT a re-rounding of the underlying value).
.supp_col <- function(name, label = name, type = "character",
                      italic = FALSE, digits = NULL, width = NULL,
                      align = NULL) {
  list(name = name, label = label, type = type,
       italic = isTRUE(italic), digits = digits, width = width,
       align = align)
}
col     <- function(name, label = name, italic = FALSE, width = NULL,
                    align = NULL)
  .supp_col(name, label, "character", italic, NULL, width, align)
col_int <- function(name, label = name, italic = FALSE, width = NULL,
                    align = "right")
  .supp_col(name, label, "integer", italic, NULL, width, align)
col_num <- function(name, label = name, italic = FALSE, digits = 2,
                    width = NULL, align = "right")
  .supp_col(name, label, "numeric", italic, digits, width, align)

# ── Validator: hard-fail on any unexpected cell value ─────────────────────────
.validate_supp_data <- function(data, columns, sheet_label = "") {
  if (!is.data.frame(data) || nrow(data) == 0L)
    stop("write_supp_table(): `data` must be a non-empty data.frame ",
         "(sheet = '", sheet_label, "').", call. = FALSE)

  spec_names <- vapply(columns, `[[`, character(1), "name")
  missing_cols <- setdiff(spec_names, names(data))
  if (length(missing_cols))
    stop("write_supp_table(): columns missing from `data` (sheet = '",
         sheet_label, "'): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  for (spec in columns) {
    x <- data[[spec$name]]
    nm <- spec$name

    if (any(is.na(x))) {
      bad <- which(is.na(x))
      stop(sprintf(
        "write_supp_table(): NA values in column '%s' (sheet = '%s'). ",
        nm, sheet_label),
        sprintf("Row(s): %s. ",
                paste(head(bad, 20), collapse = ", ")),
        "Convert intentional NAs to a sentinel string ",
        "(e.g. '—' or 'n<30') and declare the column as character.",
        call. = FALSE)
    }

    if (is.numeric(x)) {
      if (any(is.nan(x)))
        stop(sprintf("write_supp_table(): NaN in column '%s' (sheet = '%s') ",
                     nm, sheet_label),
             "at row(s) ",
             paste(head(which(is.nan(x)), 20), collapse = ", "),
             ". Likely a 0/0 division. Fix at the source.", call. = FALSE)
      if (any(is.infinite(x)))
        stop(sprintf("write_supp_table(): Inf in column '%s' (sheet = '%s') ",
                     nm, sheet_label),
             "at row(s) ",
             paste(head(which(is.infinite(x)), 20), collapse = ", "), ".",
             call. = FALSE)
    }

    if (is.character(x)) {
      empty <- which(!nzchar(trimws(x)))
      if (length(empty))
        stop(sprintf(
          "write_supp_table(): empty / whitespace-only strings in column '%s' (sheet = '%s'). ",
          nm, sheet_label),
          "Row(s): ", paste(head(empty, 20), collapse = ", "), ".",
          call. = FALSE)
    }

    # Type check vs spec
    if (spec$type == "character") {
      # Allow factor (coerced on write); error on numeric/logical.
      if (!(is.character(x) || is.factor(x)))
        stop(sprintf(
          "write_supp_table(): column '%s' declared character, got %s.",
          nm, class(x)[1]), call. = FALSE)
    } else if (spec$type == "integer") {
      if (!is.numeric(x))
        stop(sprintf(
          "write_supp_table(): column '%s' declared integer, got %s.",
          nm, class(x)[1]), call. = FALSE)
      # Non-integer numeric → error (would be misleading under integer display).
      if (any(abs(x - round(x)) > 1e-8))
        stop(sprintf(
          "write_supp_table(): column '%s' declared integer but contains non-integer values ",
          nm),
          "(e.g. ", paste(head(x[abs(x - round(x)) > 1e-8], 3),
                          collapse = ", "), "). ",
          "Use col_num() with digits = N instead.", call. = FALSE)
    } else if (spec$type == "numeric") {
      if (!is.numeric(x))
        stop(sprintf(
          "write_supp_table(): column '%s' declared numeric, got %s.",
          nm, class(x)[1]), call. = FALSE)
    } else {
      stop("write_supp_table(): unknown spec type '", spec$type, "' for column '",
           nm, "'.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ── Single-sheet writer ───────────────────────────────────────────────────────
write_supp_table <- function(data, path, columns, sheet = "Sheet1",
                             title = NULL, footnote = NULL,
                             freeze_header = TRUE) {
  wb <- createWorkbook()
  .add_supp_sheet(wb, data = data, columns = columns, sheet = sheet,
                  title = title, footnote = footnote,
                  freeze_header = freeze_header)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

# ── Multi-sheet writer (for ST3, ST8, ST9, ST10, ST11 which carry > 1 sheet)
write_supp_table_workbook <- function(sheets, path) {
  if (!is.list(sheets) || is.null(names(sheets)) || any(!nzchar(names(sheets))))
    stop("write_supp_table_workbook(): `sheets` must be a named list. Names become sheet tabs.",
         call. = FALSE)
  wb <- createWorkbook()
  for (sh_name in names(sheets)) {
    sh <- sheets[[sh_name]]
    if (!all(c("data", "columns") %in% names(sh)))
      stop("write_supp_table_workbook(): sheet '", sh_name,
           "' must have `data` and `columns`.", call. = FALSE)
    .add_supp_sheet(wb,
                    data = sh$data, columns = sh$columns,
                    sheet = sh_name,
                    title = sh$title, footnote = sh$footnote,
                    freeze_header = sh$freeze_header %||% TRUE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Worker: render one sheet ──────────────────────────────────────────────────
.add_supp_sheet <- function(wb, data, columns, sheet,
                            title = NULL, footnote = NULL,
                            freeze_header = TRUE) {
  .validate_supp_data(data, columns, sheet_label = sheet)

  # Reorder + rename data per spec.
  out <- data[, vapply(columns, `[[`, character(1), "name"), drop = FALSE]
  names(out) <- vapply(columns, `[[`, character(1), "label")

  # Excel tab name: trim to 31 chars (Excel hard limit) and replace illegal chars.
  tab <- gsub("[\\\\/?*\\[\\]:]", "_", sheet)
  tab <- substr(tab, 1, 31)
  addWorksheet(wb, tab, gridLines = FALSE)

  # Optional title row above the table.
  start_row <- 1L
  if (!is.null(title)) {
    writeData(wb, tab, title, startRow = start_row, startCol = 1)
    addStyle(wb, tab,
             createStyle(textDecoration = "bold", fontSize = 12),
             rows = start_row, cols = 1)
    start_row <- start_row + 2L  # blank gap row
  }
  header_row <- start_row
  data_start <- header_row + 1L

  # Write data block (incl. header row).
  writeData(wb, tab, out, startRow = header_row, startCol = 1,
            headerStyle = NULL)

  # Header styling: bold + thin bottom border + center, italic where requested.
  header_base <- createStyle(textDecoration = "bold",
                             border = "bottom", borderStyle = "thin",
                             halign = "center", valign = "center",
                             wrapText = TRUE)
  header_italic <- createStyle(textDecoration = c("bold", "italic"),
                               border = "bottom", borderStyle = "thin",
                               halign = "center", valign = "center",
                               wrapText = TRUE)
  for (j in seq_along(columns)) {
    sty <- if (isTRUE(columns[[j]]$italic)) header_italic else header_base
    addStyle(wb, tab, sty, rows = header_row, cols = j,
             gridExpand = FALSE, stack = FALSE)
  }
  # Thin top rule across the header so the table has a JAMA-/Nature-style frame.
  addStyle(wb, tab,
           createStyle(border = "top", borderStyle = "thin"),
           rows = header_row, cols = seq_along(columns),
           gridExpand = TRUE, stack = TRUE)

  # Body alignment + numeric display formats.
  for (j in seq_along(columns)) {
    spec <- columns[[j]]
    align <- spec$align %||% (if (spec$type == "character") "left" else "right")
    body_sty <- createStyle(halign = align, valign = "center")
    if (spec$type == "integer") {
      body_sty <- createStyle(halign = align, valign = "center",
                              numFmt = "0")
    } else if (spec$type == "numeric" && !is.null(spec$digits)) {
      fmt <- if (spec$digits == 0L) "0" else
        paste0("0.", strrep("0", spec$digits))
      body_sty <- createStyle(halign = align, valign = "center",
                              numFmt = fmt)
    }
    addStyle(wb, tab, body_sty,
             rows = data_start:(data_start + nrow(out) - 1L),
             cols = j, gridExpand = FALSE, stack = FALSE)
  }

  # Thin bottom rule under the last data row to close the frame.
  addStyle(wb, tab,
           createStyle(border = "bottom", borderStyle = "thin"),
           rows = data_start + nrow(out) - 1L,
           cols = seq_along(columns),
           gridExpand = TRUE, stack = TRUE)

  # Column widths: explicit > auto.
  for (j in seq_along(columns)) {
    w <- columns[[j]]$width
    if (!is.null(w)) {
      setColWidths(wb, tab, cols = j, widths = w)
    } else {
      label_w <- nchar(columns[[j]]$label)
      data_w  <- if (columns[[j]]$type == "character") {
        max(nchar(as.character(out[[j]])), na.rm = TRUE)
      } else {
        max(nchar(format(out[[j]])), na.rm = TRUE)
      }
      setColWidths(wb, tab, cols = j,
                   widths = max(label_w, data_w) + 3L)
    }
  }

  # Footnotes intentionally suppressed (2026-06-15): Nature MOESM3/MOESM10
  # templates carry zero in-sheet footnotes — all methodology lives in the
  # manuscript caption (captions.R). Producers still pass `footnote = ...`
  # for documentation, but it's ignored at write time.
  invisible(footnote)

  if (isTRUE(freeze_header))
    freezePane(wb, tab, firstActiveRow = data_start, firstActiveCol = 1)

  invisible(TRUE)
}
