#!/usr/bin/env Rscript
# 00_mask_phi.R
# -----------------------------------------------------------------------------
# PHI-masking ingestion step. Reads the real-MRN / real-accession master
# spreadsheets from
#   data/raw/phi_unsafe/   (gitignored; never committed)
# masks all MRNs and pathology accessions in place, and writes identical-schema
# copies to
#   data/raw/phi_safe/     (committed; safe to share)
# It also rewrites data/flags/review_ids.csv so its MRN-based audit flags
# match the masked identifiers in the safe spreadsheets.
#
# Masking:
#   MRN → "#####XXXX"            (9-char uniform; first 5 digits preserved)
#   Pathology accession → "NP-YY-##XXX"
#                               (year + first 2 digits of serial; rest masked)
#                               If two accessions collide under this rule (e.g.
#                               two 4-digit serials sharing first 2 digits),
#                               a trailing letter ('a','b',...) is appended,
#                               sorted by the original string, to guarantee
#                               a bijection from original → masked.
#
# The masked MRN column loses the last 2-5 digits of precision; the masked
# accession loses the last 1-3 digits of the serial (plus letter disambig
# where needed). Everything else in the spreadsheet passes through unchanged.
#
# Downstream: analysis/00_data_prep/01_clean_master.R reads only from phi_safe/.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl); library(writexl); library(dplyr); library(here)
})

unsafe_dir <- here("data", "raw", "phi_unsafe")
safe_dir   <- here("data", "raw", "phi_safe")
dir.create(safe_dir, recursive = TRUE, showWarnings = FALSE)

# ── MRN mask ─────────────────────────────────────────────────────────────────
mask_mrn <- function(v) {
  v <- as.character(v)
  d <- gsub("[^0-9]", "", v)
  ifelse(is.na(v) | nchar(d) == 0, v, paste0(substr(d, 1, 5), "XXXX"))
}

# ── Accession mask ───────────────────────────────────────────────────────────
# Build a deterministic collision-free mapping over the FULL set of accessions
# found across all phi_unsafe spreadsheets, then reuse that mapping for every
# (sheet, column) substitution. Using a single global mapping means the masked
# form of any given accession is identical everywhere it appears.
build_acc_map <- function(accs) {
  accs <- unique(accs[nchar(accs) > 0])
  masked <- gsub("(NP-[0-9]{2}-[0-9]{2})[0-9]+", "\\1XXX", accs, perl = TRUE)
  # Detect collisions and disambiguate with a sorted-by-original letter suffix.
  dup_masks <- unique(masked[duplicated(masked)])
  for (m in dup_masks) {
    ix  <- which(masked == m)
    ord <- order(accs[ix])
    masked[ix[ord]] <- paste0(m, letters[seq_along(ord)])
  }
  stopifnot(!anyDuplicated(masked))
  setNames(masked, accs)
}

apply_acc_mask <- function(text, acc_map) {
  if (is.null(acc_map) || !length(acc_map)) return(text)
  # Longest originals first so shorter aren't clobbered as substrings.
  ord <- order(-nchar(names(acc_map)))
  for (orig in names(acc_map)[ord]) {
    text <- gsub(orig, acc_map[[orig]], text, fixed = TRUE)
  }
  text
}

# Collect all accessions across every phi_unsafe xlsx up front (so the global
# mapping is stable even if one spreadsheet has a subset).
collect_accs <- function(dir) {
  files <- list.files(dir, pattern = "\\.xlsx$", full.names = TRUE)
  all_text <- unlist(lapply(files, function(f) {
    sheets <- excel_sheets(f)
    unlist(lapply(sheets, function(sh) {
      df <- read_excel(f, sheet = sh, col_types = "text")
      unlist(df, use.names = FALSE)
    }))
  }))
  unique(unlist(regmatches(all_text, gregexpr("NP-[0-9]{2}-[0-9]+", all_text))))
}

cat("== Scanning phi_unsafe for pathology accessions ==\n")
acc_universe <- collect_accs(unsafe_dir)
acc_map <- build_acc_map(acc_universe)
cat(sprintf("  %d unique accessions, %d unique masked values\n",
            length(acc_map), length(unique(unname(acc_map)))))
ambig <- acc_map[duplicated(substr(unname(acc_map), 1, 11)) |
                 duplicated(substr(unname(acc_map), 1, 11), fromLast = TRUE)]
if (length(ambig) > 0) {
  cat("  collision-resolution letter suffixes:\n")
  for (o in names(ambig)) cat(sprintf("    %s → %s\n", o, ambig[[o]]))
}

# ── Process each xlsx ────────────────────────────────────────────────────────
check_no_mrn_collisions <- function(orig, masked, label) {
  df <- tibble(orig = orig, masked = masked) |>
    filter(!is.na(orig), !is.na(masked), nzchar(orig), nzchar(masked))
  agg <- df |> distinct(orig, masked) |>
    count(masked, name = "n_distinct_orig") |>
    filter(n_distinct_orig > 1)
  if (nrow(agg) > 0) {
    print(agg)
    stop("MRN collision in ", label, call. = FALSE)
  }
}

process_xlsx <- function(filename) {
  src <- file.path(unsafe_dir, filename)
  dst <- file.path(safe_dir,   filename)
  sheets <- excel_sheets(src)
  out <- setNames(vector("list", length(sheets)), sheets)
  for (sh in sheets) {
    # Read with default col_types so numeric phenotype columns (VAF, age,
    # codes) stay numeric for the downstream cleaner's case_when branches.
    # Identifier columns are coerced to character only where we mask.
    df <- read_excel(src, sheet = sh)
    if ("MRN" %in% names(df)) {
      orig <- as.character(df$MRN)
      df$MRN <- mask_mrn(orig)
      check_no_mrn_collisions(orig, df$MRN, sprintf("%s[%s]:MRN", filename, sh))
    }
    # Mask accessions in every character column (patient_id lives under
    # varying header spellings; Notes may also carry inline accessions).
    # Non-character columns (numeric/logical) are left untouched — accessions
    # don't appear there.
    for (col in names(df)) {
      if (is.character(df[[col]])) df[[col]] <- apply_acc_mask(df[[col]], acc_map)
    }
    out[[sh]] <- df
  }
  write_xlsx(out, dst)
  cat(sprintf("  %s  (%d sheet%s) → phi_safe/%s\n",
              filename, length(sheets), if (length(sheets)==1) "" else "s",
              filename))
}

cat("\n== Masking phi_unsafe/*.xlsx → phi_safe/*.xlsx ==\n")
for (f in list.files(unsafe_dir, pattern = "\\.xlsx$", full.names = FALSE)) {
  process_xlsx(f)
}

# ── review_ids.csv ───────────────────────────────────────────────────────────
review_path <- here("data", "flags", "review_ids.csv")
if (file.exists(review_path)) {
  ri <- read.csv(review_path, stringsAsFactors = FALSE)
  if ("flag_type" %in% names(ri) && "value" %in% names(ri)) {
    mrn_rows <- ri$flag_type == "conflict_mrn"
    if (any(mrn_rows)) ri$value[mrn_rows] <- mask_mrn(ri$value[mrn_rows])
    ri$value <- apply_acc_mask(ri$value, acc_map)
    ri$description <- apply_acc_mask(ri$description, acc_map)
    write.csv(ri, review_path, row.names = FALSE, quote = TRUE)
    cat(sprintf("\n== Masked data/flags/review_ids.csv ==\n"))
  }
}

cat("\nDone.\n")
