# =============================================================================
# check_panel_uniqueness.R — fail render on conflicting on-disk panels
# -----------------------------------------------------------------------------
# Walks every per-panel PNG/PDF under results/Figure*/ and results/ExtendedData/
# and validates four invariants:
#
#   (1) UNIQUE   — each registered panel token has files at exactly one
#                  on-disk path. If two producers wrote the same token to
#                  two different paths, the gate fails (this is the bug
#                  class behind the Figure3/panel_A duplicate issue:
#                  09_F1_km_age.R wrote Fig3A_km_presentation.png while
#                  14_ED7_vaf_phenotype.R wrote Fig3A_vaf_sm_total.png to
#                  the same panel_A directory; same dir, different
#                  filenames, so no fs collision but the slot held two
#                  panels at once).
#
#   (2) REGISTERED — every PNG/PDF under a panel dir must derive a token
#                    that is in panel_registry$figures. Stray files (e.g.
#                    Fig3D_vaf_eloquence.pdf left over after a rename
#                    where no producer writes it anymore) trip the gate.
#
#   (3) PLACED    — each file must live in its token's canonical group
#                   dir (Main: results/FigureN/, ED:
#                   results/ExtendedData/ed_<token>/). A file whose token
#                   says it belongs to ed_vaf_phenotype but sits in
#                   Figure3/panel_A/ trips the gate.
#
#   (4) NO_PREFIX — modern policy: panel filenames should be the bare
#                   registry token (e.g. parietal_kras_rupture.png), not
#                   the legacy "FigNX_<token>.png" form. The gate accepts
#                   the legacy form for now (just strips the prefix when
#                   parsing) but emits a soft warning so the disk
#                   gradually drifts toward the cleaner naming.
#
# Strictness: BAVM_PANEL_UNIQUE=off bypasses (debugging only). Default
# strict; (1)+(2)+(3) are hard errors, (4) warns.
# =============================================================================

suppressPackageStartupMessages({ library(here) })

# Group token -> on-disk dir name under results/. Main figures use the
# legacy "FigureN" naming (the producers write there); ED groups use the
# canonical ed_<token> naming.
#
# The mapping is derived dynamically from the panel-assignments cache
# (results/stats/panel_assignments.rds, populated by
# analysis/pipeline/helpers/panel_assignments.R from first-citation order in the
# resolved prose). This means the mapping always tracks the current
# manuscript figure-number assignment — if a re-citation moves the `age`
# group from Fig. 3 to Fig. 2, the disk dir Figure2/ that holds its
# panels stops looking misplaced without anyone hand-editing this file.
#
# Falls back to a static registry-declaration-order mapping only when the
# cache is missing (first run, or before run_all.R [1.5] has rebuilt it).
.main_fig_group_to_dir <- function() {
  cache_path <- file.path(here::here(), "results", "stats",
                          "panel_assignments.rds")
  if (file.exists(cache_path)) {
    a <- readRDS(cache_path)
    fig_nums <- a$group_number$Fig
    if (!is.null(fig_nums)) {
      out <- as.list(sprintf("Figure%d", as.integer(unname(fig_nums))))
      names(out) <- names(fig_nums)
      return(out)
    }
  }
  # Fallback: registry-declaration order (used pre-cache).
  list(
    cohort_variants     = "Figure1",
    null_phenotype      = "Figure2",
    age                 = "Figure3",
    venous_stenosis     = "Figure4",
    liquid_biopsy       = "Figure5",
    precision_framework = "Figure6"
  )
}

# Strip the resolver-prefix from filename basenames so the remainder
# matches a registry token. Two prefix grammars:
#   "FigNX_<token>"   main-figure human-navigation copies
#   "EDNNX_<token>"   extended-data human-navigation copies (sweeper-managed)
# Both are written by analysis/pipeline/helpers/sync_panel_prefixes.R alongside the
# canonical token-only file and are NOT legacy. The `had_prefix` flag is
# kept for the soft-warning path, but ED prefixes are not flagged because
# they are the current convention.
.parse_panel_token <- function(stem) {
  rx <- "^(Fig\\d+[A-Z]?_|ED\\d+[A-Z]?_)"
  m <- regmatches(stem, regexpr(rx, stem, perl = TRUE))
  if (length(m) == 1L && nzchar(m)) {
    list(
      token      = sub(rx, "", stem, perl = TRUE),
      had_prefix = startsWith(m, "Fig")
    )
  } else {
    list(token = stem, had_prefix = FALSE)
  }
}

.expected_dir_for_token <- function(tok, registry, group_to_dir = NULL) {
  if (is.null(group_to_dir)) group_to_dir <- .main_fig_group_to_dir()
  for (grp_name in names(registry$figures)) {
    grp <- registry$figures[[grp_name]]
    if (tok %in% grp$panels) {
      track <- grp$track
      if (track == "Fig") {
        d <- group_to_dir[[grp_name]]
        if (is.null(d)) {
          stop(sprintf(
            "[panel-uniqueness] no disk-dir mapping for main-fig group '%s'; ",
            "regenerate panel_assignments.rds via Rscript analysis/pipeline/helpers/panel_assignments.R"
          ), call. = FALSE)
        }
        return(list(group = grp_name, expected_dir = file.path("results", d)))
      } else if (track == "EDFig") {
        return(list(
          group = grp_name,
          expected_dir = file.path("results", "ExtendedData", grp_name)
        ))
      } else {
        stop(sprintf(
          "[panel-uniqueness] unknown track '%s' for group '%s'", track, grp_name
        ), call. = FALSE)
      }
    }
  }
  NULL
}

# Files we scan: any .png or .pdf under
#   results/Figure*/panel_*/
#   results/Figure*/                          (composite + flat layouts)
#   results/ExtendedData/ed_*/                (canonical ED dirs)
#   results/ExtendedData/ed_*/panel_*/        (ED panel subdirs)
# Skip _retired_*, _exploration/, the composite Fig{N}_composite files
# (they aren't per-panel), and Fig5/Fig6 dirs whose external assets carry
# a different naming convention (Figure5.png / Figure6.png).
.gather_panel_files <- function(root = NULL) {
  if (is.null(root)) root <- file.path(here::here(), "results")
  patterns <- c(
    "Figure1/panel_*/*.png", "Figure1/panel_*/*.pdf",
    "Figure2/panel_*/*.png", "Figure2/panel_*/*.pdf",
    "Figure3/panel_*/*.png", "Figure3/panel_*/*.pdf",
    "Figure4/panel_*/*.png", "Figure4/panel_*/*.pdf",
    "ExtendedData/ed_*/*.png", "ExtendedData/ed_*/*.pdf",
    "ExtendedData/ed_*/panel_*/*.png", "ExtendedData/ed_*/panel_*/*.pdf"
  )
  files <- unlist(lapply(patterns, function(p) {
    Sys.glob(file.path(root, p))
  }), use.names = FALSE)
  # Composite filenames live at the group-dir root and are not per-panel
  # artefacts. Filter them out:
  #   Main:  Figure3/Fig3_composite.{png,pdf}
  #   ED:    ed_X/anatomy_composite.{png,pdf} (named composites)
  #          ed_X/<anything>_combined.{png,pdf}
  #          ed_X/<anything>_composite.{png,pdf}
  files <- files[!grepl(
    "/(Fig\\d+_composite|[^/]+_combined|[^/]+_composite)\\.(png|pdf)$",
    files, perl = TRUE)]
  # Files inside per-panel-letter subdirs of ED groups (ed_X/panel_Y/) are
  # legitimate; nothing to filter further.
  files
}

check_panel_uniqueness <- function(registry = NULL, strict = NULL) {
  if (is.null(registry)) registry <- load_panel_registry()
  if (is.null(strict)) {
    strict <- !identical(tolower(Sys.getenv("BAVM_PANEL_UNIQUE")), "off")
  }

  files <- .gather_panel_files()
  root <- file.path(here::here(), "results")

  # Aggregate: token -> character vector of repo-relative paths.
  token_files  <- list()    # token -> char(N)
  unregistered <- character(0L)
  misplaced    <- character(0L)
  legacy_prefixed <- character(0L)

  for (f in files) {
    rel <- sub(paste0("^", here::here(), "/?"), "", f)
    stem <- tools::file_path_sans_ext(basename(f))
    parsed <- .parse_panel_token(stem)
    tok <- parsed$token
    if (parsed$had_prefix) legacy_prefixed <- c(legacy_prefixed, rel)

    look <- .expected_dir_for_token(tok, registry)
    if (is.null(look)) {
      unregistered <- c(unregistered, sprintf("%s  (token: %s)", rel, tok))
      next
    }

    parent_dir_rel <- sub(paste0("^", root, "/?"), "results/", dirname(f))
    # `panel_X` subdirs are allowed under both Main and ED groups; collapse
    # them so a path like results/Figure3/panel_A counts as inside Figure3.
    parent_collapsed <- sub("/panel_[A-Za-z]+$", "", parent_dir_rel)
    if (parent_collapsed != look$expected_dir) {
      misplaced <- c(misplaced, sprintf(
        "%s  (token '%s' belongs in %s)", rel, tok, look$expected_dir
      ))
      next
    }

    token_files[[tok]] <- c(token_files[[tok]], rel)
  }

  # (1) Token claimed by files at >1 distinct *parent dir* is the duplicate
  # case. Same token at the same path emitted as both .png and .pdf is OK.
  duplicate_tokens <- list()
  for (tok in names(token_files)) {
    parents <- unique(dirname(token_files[[tok]]))
    if (length(parents) > 1L) {
      duplicate_tokens[[tok]] <- token_files[[tok]]
    }
  }

  # ---- report --------------------------------------------------------------
  hard_errors <- character(0L)
  if (length(duplicate_tokens) > 0L) {
    msg <- c("PANEL DUPLICATES (same token, multiple paths):")
    for (tok in names(duplicate_tokens)) {
      msg <- c(msg, sprintf("  ✖ %s ->", tok),
               paste0("      ", duplicate_tokens[[tok]]))
    }
    hard_errors <- c(hard_errors, paste(msg, collapse = "\n"))
  }
  if (length(unregistered) > 0L) {
    hard_errors <- c(hard_errors, paste(c(
      "UNREGISTERED PANEL FILES (filename-stem not in panel_registry):",
      paste0("  ✖ ", unregistered)
    ), collapse = "\n"))
  }
  if (length(misplaced) > 0L) {
    hard_errors <- c(hard_errors, paste(c(
      "MISPLACED PANEL FILES (token registered in different group):",
      paste0("  ✖ ", misplaced)
    ), collapse = "\n"))
  }

  if (length(hard_errors) > 0L) {
    # `message()` (not stop()) for the body so very long violation lists
    # aren't truncated by R's error-string length cap. The eventual
    # stop()/warning() carries only a short summary.
    message("")
    for (block in hard_errors) message("[panel-uniqueness] ", block)
    message("")
    summary_line <- sprintf(
      "panel-uniqueness FAILED: %d duplicate token(s), %d unregistered file(s), %d misplaced file(s).  Bypass: BAVM_PANEL_UNIQUE=off",
      length(duplicate_tokens), length(unregistered), length(misplaced)
    )
    if (strict) stop(summary_line, call. = FALSE)
    else        warning(summary_line, call. = FALSE)
  } else {
    message(sprintf("[panel-uniqueness] %d file(s) -> %d token(s) -> all unique, registered, and placed%s",
                    length(files), length(token_files),
                    if (length(legacy_prefixed) > 0L)
                      sprintf(" (%d still using legacy FigNX_ prefix)",
                              length(legacy_prefixed))
                    else ""))
  }

  invisible(list(
    n_files          = length(files),
    n_tokens         = length(token_files),
    duplicate_tokens = duplicate_tokens,
    unregistered     = unregistered,
    misplaced        = misplaced,
    legacy_prefixed  = legacy_prefixed,
    strict           = strict
  ))
}

# CLI guard so sourcing from run_all.R does not auto-execute the check.
.this_script <- "check_panel_uniqueness.R"
.invoked_as_script <- any(grepl(
  paste0("--file=.*/", .this_script, "$|--file=", .this_script, "$"),
  commandArgs(trailingOnly = FALSE)
))
if (.invoked_as_script) {
  suppressPackageStartupMessages({
    source(here::here("analysis", "pipeline", "helpers", "resolve_panels.R"))
  })
  check_panel_uniqueness()
}
