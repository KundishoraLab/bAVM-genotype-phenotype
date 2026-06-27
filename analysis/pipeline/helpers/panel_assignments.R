# =============================================================================
# panel_assignments.R — registry-derived lookup API for producers + captions
# -----------------------------------------------------------------------------
# Phase 2 entry point. Exposes three lookup functions that producer scripts
# and caption functions call to get the currently-assigned figure/table
# numbers and panel letters:
#
#   figure_number(group)       -> integer  (the manuscript-facing Fig N)
#   panel_letter(token)        -> "a", "b", ... "" for single-panel group,
#                                  NA if the token exists in the registry
#                                  but is never cited in prose
#   table_number(token)        -> integer
#
# Values are driven by the panel_registry.R catalog and the narrative prose:
# first appearance in prose order assigns the number/letter. The computation
# is deterministic, pure, and cheap, but we still cache it to
# `results/stats/panel_assignments.rds` so producer scripts running in
# sequence do not re-parse the whole prose tree.
#
# Workflow:
#   1. run_all.R step [1.5] calls compute_panel_assignments() and writes the
#      cache (via resolve_sections_to_disk(), which invokes us internally).
#   2. Producer scripts (step [2]) call load_panel_assignments() once and
#      pass the result into their lookup calls:
#        assignments <- load_panel_assignments()
#        tag <- panel_letter("km_rupture", assignments)
#   3. Caption functions (loaded from _setup.qmd) call the lookup helpers
#      without arguments; the cache is read lazily on first call.
#
# Standalone: `Rscript analysis/pipeline/helpers/panel_assignments.R` regenerates the
# cache without touching the _resolved/ mirror or running analysis scripts.
# =============================================================================

suppressPackageStartupMessages(library(here))
source(here::here("analysis", "pipeline", "helpers", "resolve_panels.R"))

.default_assignments_path <- function() {
  file.path(here::here(), "results", "stats", "panel_assignments.rds")
}

# Walks the registry + narrative prose and returns the assignment maps.
# This is the single source of truth for manuscript-facing numbers/letters.
# Both the prose resolver and producer scripts consume its output.
compute_panel_assignments <- function(registry       = NULL,
                                      section_paths  = NULL,
                                      registry_path  = NULL) {
  if (is.null(registry)) registry <- load_panel_registry(registry_path)
  if (is.null(section_paths)) {
    sections_dir <- file.path(here::here(), "manuscript", "submission_prep", "sections")
    # Order matches main_text.qmd includes so first-citation order is stable.
    section_order <- c(
      "01_title_and_abstract",
      "02_introduction",
      "04_results",
      "05_discussion",
      "03_methods",
      "06_figure_legends"
    )
    section_paths <- file.path(sections_dir, paste0(section_order, ".qmd"))
    section_paths <- section_paths[file.exists(section_paths)]
  }

  token_index <- build_token_index(registry)

  all_citations <- list()
  for (i in seq_along(section_paths)) {
    txt <- paste(readLines(section_paths[i], warn = FALSE), collapse = "\n")
    # Strip comments / code fences so hidden tokens do not leak into
    # assignments. Mirrors resolve_panels.R's prose-scanning hygiene.
    txt <- gsub("(?s)<!--.*?-->", " ", txt, perl = TRUE)
    txt <- gsub("(?s)```.*?```",  " ", txt, perl = TRUE)
    txt <- gsub("`[^`\n]+`",      " ", txt, perl = TRUE)
    parsed <- parse_citations(txt, section_ord = i)
    all_citations <- c(all_citations, parsed$rows)
  }

  if (length(all_citations) > 0L) {
    ord <- order(
      vapply(all_citations, `[[`, integer(1L), "section_ord"),
      vapply(all_citations, `[[`, numeric(1L), "pos")
    )
    all_citations <- all_citations[ord]
  }

  maps <- assign_numbers(all_citations, token_index)
  maps$meta <- list(
    computed_at = Sys.time(),
    n_sections  = length(section_paths),
    n_cites     = length(all_citations),
    n_tokens    = length(token_index)
  )
  maps
}

write_panel_assignments <- function(assignments, path = NULL) {
  if (is.null(path)) path <- .default_assignments_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(assignments, path)
  invisible(path)
}

load_panel_assignments <- function(path = NULL) {
  if (is.null(path)) path <- .default_assignments_path()
  if (!file.exists(path)) {
    stop(
      "panel assignments cache not found at ", path, "\n",
      "  Regenerate with: Rscript analysis/pipeline/helpers/panel_assignments.R",
      call. = FALSE
    )
  }
  readRDS(path)
}

# ---- lookup helpers --------------------------------------------------------
# Each helper takes an optional `assignments` argument for producer scripts
# that have already loaded the cache (avoids re-reading rds on every call).
# Caption functions, called via Quarto, pass no argument and rely on a
# module-level cache that loads on first use.

.cached_assignments <- NULL

.get_assignments <- function(assignments = NULL) {
  if (!is.null(assignments)) return(assignments)
  if (is.null(.cached_assignments)) {
    .cached_assignments <<- load_panel_assignments()
  }
  .cached_assignments
}

figure_number <- function(group, assignments = NULL) {
  a <- .get_assignments(assignments)
  for (track in names(a$group_number)) {
    if (group %in% names(a$group_number[[track]])) {
      return(unname(a$group_number[[track]][[group]]))
    }
  }
  NA_integer_
}

panel_letter <- function(token, assignments = NULL) {
  a <- .get_assignments(assignments)
  for (grp in names(a$panel_letter)) {
    if (token %in% names(a$panel_letter[[grp]])) {
      return(unname(a$panel_letter[[grp]][[token]]))
    }
  }
  NA_character_
}

table_number <- function(token, assignments = NULL) {
  a <- .get_assignments(assignments)
  for (track in names(a$table_number)) {
    if (token %in% names(a$table_number[[track]])) {
      return(unname(a$table_number[[track]][[token]]))
    }
  }
  NA_integer_
}

# Convenience for caption authors: returns "a", "b", ... or "" for
# single-panel groups. Useful when splicing `**%s,**` into panel descriptions
# where single-panel figures should render as plain `**Fig. N** ...`
# rather than `**Fig. N** **,** ...`.
panel_tag <- function(token, assignments = NULL) {
  # Find the group for this token.
  a <- .get_assignments(assignments)
  for (grp in names(a$panel_letter)) {
    if (token %in% names(a$panel_letter[[grp]])) {
      # Count panels actually assigned letters in this group. If only one,
      # the caption should omit the letter entirely.
      if (length(a$panel_letter[[grp]]) == 1L) return("")
      return(unname(a$panel_letter[[grp]][[token]]))
    }
  }
  NA_character_
}

# ---- CLI entry point -------------------------------------------------------

.this_script <- "panel_assignments.R"
.invoked_as_script <- any(grepl(
  paste0("--file=.*/", .this_script, "$|--file=", .this_script, "$"),
  commandArgs(trailingOnly = FALSE)
))
if (.invoked_as_script) {
  info <- compute_panel_assignments()
  path <- write_panel_assignments(info)
  rel  <- sub(paste0("^", here::here(), "/?"), "", path)
  message(sprintf(
    "✓ panel assignments written to %s (%d tokens, %d first-citations)",
    rel, info$meta$n_tokens, info$meta$n_cites
  ))
}
