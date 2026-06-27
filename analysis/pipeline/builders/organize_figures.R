#!/usr/bin/env Rscript
# =============================================================================
# organize_figures.R — assemble the manuscript-facing figure/table deliverable.
#
# Reads the producer outputs scattered under results/ and writes a clean,
# manuscript-labeled tree to results/figures/ — one folder per figure, panels
# named by figure number + panel letter:
#
#   results/figures/
#     Figure 1/
#       Fig 1 (composite).{pdf,png}
#       Fig 1A.{pdf,png}              (consort_flow)
#       Fig 1B.{pdf,png}              (variant_landscape)
#       …
#     Figure 2/ … Figure 4/
#     Extended Data Fig 1/
#       ED Fig 1.{pdf,png}            (single-panel → no letter)
#     Extended Data Fig 6/
#       ED Fig 6 (composite).{pdf,png}
#       ED Fig 6A.{pdf,png} …
#     Tables/
#       Table 1/…
#       Supplementary Tables/…
#
# Figure numbers and panel letters come from the frozen assignment cache
# (results/stats/panel_assignments.rds). Idempotent: wipes results/figures/
# and rebuilds. Run as the last step of run_all.R, or standalone:
#   Rscript analysis/pipeline/builders/organize_figures.R
# =============================================================================

suppressPackageStartupMessages(library(here))
source(here::here("analysis", "pipeline", "helpers", "panel_assignments.R"))

OUT_ROOT    <- here::here("results", "figures")
SEARCH_ROOT <- here::here("results")
EXCLUDE     <- c("/figures/", "/_review/", "/_exploration/", "/_archive/",
                 "/_retired", "/_supporting/", "/_parity/",
                 "/_manifest_fragments/", "/stats/")

assignments <- load_panel_assignments()

# ---- file discovery (mirrors aggregate_review_dir.R resolution) ------------
discover_files <- function() {
  hits <- list.files(SEARCH_ROOT, pattern = "\\.(pdf|png|svg)$",
                     recursive = TRUE, full.names = TRUE)
  for (pat in EXCLUDE) hits <- hits[!grepl(pat, hits, fixed = TRUE)]
  hits
}

# files whose basename ends in `<token>.{ext}` (handles Fig1B_<token>.pdf and <token>.pdf)
find_panel_hits <- function(token, all_files) {
  pat <- sprintf("(^|[/_-])%s\\.(pdf|png|svg)$", token)
  all_files[grepl(pat, basename(all_files))]
}

find_composite_hits <- function(group, fig_num, is_ed, all_files) {
  short <- sub("^ed_", "", group)
  num <- if (is_ed) {
    c(sprintf("/ED%02d_composite\\.(pdf|png|svg)$", fig_num),
      sprintf("/ED%d_composite\\.(pdf|png|svg)$", fig_num),
      sprintf("/EDFig%02d_composite\\.(pdf|png|svg)$", fig_num))
  } else {
    c(sprintf("/Fig%d_composite\\.(pdf|png|svg)$", fig_num),
      sprintf("/Fig%02d_composite\\.(pdf|png|svg)$", fig_num),
      # Fig 4 ships its composite as Fig04_precision_medicine.*
      sprintf("/Fig%02d_%s\\.(pdf|png|svg)$", fig_num, short))
  }
  pats <- c(sprintf("/%s_composite\\.(pdf|png|svg)$", group),
            sprintf("/%s_composite\\.(pdf|png|svg)$", short), num)
  hits <- character(0)
  for (p in pats) hits <- c(hits, all_files[grepl(p, all_files)])
  unique(hits)
}

# one file per extension; newest wins on duplicates (legacy renamed producers)
pick_one_per_ext <- function(files) {
  if (!length(files)) return(character(0))
  out <- character(0)
  for (e in unique(tools::file_ext(files))) {
    c2 <- files[tools::file_ext(files) == e]
    out <- c(out, if (length(c2) == 1L) c2 else c2[which.max(file.info(c2)$mtime)])
  }
  out
}

copy_as <- function(src_files, dest_dir, stem) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in pick_one_per_ext(src_files)) {
    file.copy(f, file.path(dest_dir, paste0(stem, ".", tools::file_ext(f))),
              overwrite = TRUE)
  }
}

# ---- build -----------------------------------------------------------------
if (dir.exists(OUT_ROOT)) unlink(OUT_ROOT, recursive = TRUE)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
all_files <- discover_files()
missing <- character(0)

for (track in names(assignments$group_number)) {            # "Fig" / "EDFig"
  is_ed   <- track %in% c("EDFig", "ED", "ed")
  folder_label <- if (is_ed) "Extended Data Fig" else "Figure"
  file_label   <- if (is_ed) "ED Fig" else "Fig"
  for (group in names(assignments$group_number[[track]])) {
    n        <- assignments$group_number[[track]][[group]]
    letters  <- assignments$panel_letter[[group]]            # token -> "a"/"b"/...
    fig_dir  <- file.path(OUT_ROOT, sprintf("%s %d", folder_label, n))
    single   <- length(letters) <= 1L

    # composite (skip for single-panel figures — the panel IS the figure)
    if (!single) {
      comp <- find_composite_hits(group, n, is_ed, all_files)
      if (length(comp)) copy_as(comp, fig_dir, sprintf("%s %d (composite)", file_label, n))
      else missing <- c(missing, sprintf("%s %d composite [%s]", file_label, n, group))
    }

    # panels
    for (token in names(letters)) {
      L    <- toupper(letters[[token]])
      hits <- find_panel_hits(token, all_files)
      stem <- if (single) sprintf("%s %d", file_label, n)
              else        sprintf("%s %d%s", file_label, n, L)
      if (length(hits)) copy_as(hits, fig_dir, stem)
      else missing <- c(missing, sprintf("%s [%s]", stem, token))
    }
  }
}

# ---- tables ----------------------------------------------------------------
tbl_root <- file.path(OUT_ROOT, "Tables")
t1 <- list.files(here::here("results", "Table1"), full.names = TRUE,
                 pattern = "\\.(xlsx|docx|pdf|png)$")
if (length(t1)) {
  dir.create(file.path(tbl_root, "Table 1"), recursive = TRUE, showWarnings = FALSE)
  file.copy(t1, file.path(tbl_root, "Table 1"), overwrite = TRUE)
}
st <- list.files(here::here("results", "SupplementaryTables"), full.names = TRUE,
                 pattern = "\\.xlsx$")
if (length(st)) {
  dir.create(file.path(tbl_root, "Supplementary Tables"), recursive = TRUE, showWarnings = FALSE)
  file.copy(st, file.path(tbl_root, "Supplementary Tables"), overwrite = TRUE)
}

# ---- report ----------------------------------------------------------------
n_files <- length(list.files(OUT_ROOT, recursive = TRUE))
message(sprintf("✓ organized %d files → results/figures/", n_files))
if (length(missing)) {
  writeLines(missing, file.path(OUT_ROOT, "MISSING.txt"))
  message(sprintf("  %d panel(s)/composite(s) had no source file — see results/figures/MISSING.txt",
                  length(missing)))
}
