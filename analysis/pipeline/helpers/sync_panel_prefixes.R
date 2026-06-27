# =============================================================================
# sync_panel_prefixes.R — keep human-navigation panel filenames live.
# -----------------------------------------------------------------------------
# Each panel saves both:
#   <token>.rds   composer-input artefact (read by 4N_compose_*.R scripts)
#   <token>.{pdf,png}   producer's standalone render
# Producers do NOT compute prefixes themselves; this sweeper takes care of
# the human-navigation layer:
#
#   <prefix>_<token>.{pdf,png}   prefixed copy that surfaces manuscript Fig N +
#                                 panel letter in the filename so you can
#                                 click around in Finder and find any panel
#                                 by its caption letter immediately.
#
# `prefix` comes from compute_panel_prefix() (see utils.R), which queries
# panel_assignments.rds and returns "Fig{N}" / "Fig{N}{X}" for main figures
# and "ED{NN}" / "ED{NN}{X}" for ED figures (single-panel groups drop the
# letter; multi-panel groups include it).
#
# What this script does (idempotent — safe to re-run):
#
#   1. Walks every <token>.{pdf,png} under results/Figure*/panel_*/ and
#      results/ExtendedData/ed_*/, results/ExtendedData/ed_*/panel_*/.
#      For each, ensures a <prefix>_<token>.{pdf,png} sibling exists with
#      the same content (file.copy if missing or older than canonical).
#
#   2. Sweeps any <prefix>_<token>.{pdf,png} whose prefix does NOT match
#      the current resolver. These are stale (e.g. left over from a Fig 2
#      <-> Fig 3 swap) and would otherwise mislead a human reviewer.
#
#   3. Reports a summary of what was added / removed.
#
# Skips: Fig*_composite.{pdf,png} at the figure root (those are composite
# outputs, not per-panel files), the legacy Fig1A_consort_flow.* (handled
# specially below), and anything under _retired_*.
#
# Run as a step in run_all.R after producers complete + the panel
# assignments cache is regenerated, or standalone:
#   Rscript analysis/pipeline/helpers/sync_panel_prefixes.R
# =============================================================================

suppressPackageStartupMessages({
  library(here)
})

source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))

# ── Discovery ────────────────────────────────────────────────────────────────

#' Return all panel-holding directories under results/.
.discover_panel_dirs <- function() {
  all <- list.dirs(here("results"), recursive = TRUE)
  # Main figure panel sub-dirs:        results/Figure{N}/panel_{X}
  # Main figure root (single-panel):   results/Figure{N}            (rare)
  # ED group dirs:                     results/ExtendedData/ed_*
  # ED panel sub-dirs:                 results/ExtendedData/ed_*/panel_{X}
  keep <- grepl(
    "/Figure[0-9]+(/panel_[A-Z]+)?$|/ExtendedData/ed_[a-z][a-z0-9_]*(/panel_[A-Z]+)?$",
    all
  )
  drop <- grepl("/_retired_|/_review/|/_exploration/", all)
  all[keep & !drop]
}

#' Parse a basename like "Fig3A_km_presentation" into prefix + token.
#' Returns NULL if the basename doesn't match the prefix grammar.
.parse_prefixed <- function(base) {
  m <- regmatches(base, regexec("^((Fig|ED)([0-9]+)([A-Z]?))_(.+)$", base))[[1]]
  if (length(m) < 6) return(NULL)
  list(
    full_prefix = m[2],
    track       = m[3],
    n_str       = m[4],
    letter      = m[5],
    token       = m[6]
  )
}

# ── Main sweep ───────────────────────────────────────────────────────────────

sync_panel_prefixes <- function(verbose = TRUE) {
  pa <- load_panel_assignments()
  panel_dirs <- .discover_panel_dirs()

  added   <- character(0)
  removed <- character(0)
  kept    <- character(0)

  for (d in panel_dirs) {
    # Zeroth pass: sweep files for tokens the registry has assigned to a
    # DIFFERENT panel dir (e.g. a C ↔ D registry swap that left the old
    # canonical sibling behind). compute_panel_prefix() alone can't catch
    # this — it would happily mint a new "ED14D_<token>" sibling next to
    # the stale plain file. Use panel_slot_dir() to resolve the canonical
    # dir for the token; if it doesn't match `d`, the file is misplaced
    # and must go (both the plain <token>.* and any prefixed siblings).
    all_files <- list.files(d, pattern = "\\.(pdf|png|rds)$", full.names = TRUE)
    for (f in all_files) {
      base   <- tools::file_path_sans_ext(basename(f))
      parsed <- .parse_prefixed(base)
      token  <- if (is.null(parsed)) base else parsed$token
      # Strip the KM-split __curve / __table suffixes before lookup.
      token  <- sub("__(curve|table)$", "", token)
      slot   <- tryCatch(panel_slot_dir(token, pa), error = function(e) NULL)
      if (is.null(slot)) next               # not in registry — leave alone
      if (normalizePath(slot, mustWork = FALSE) ==
          normalizePath(d,    mustWork = FALSE)) next
      file.remove(f)
      removed <- c(removed, f)
    }

    files <- list.files(d, pattern = "\\.(pdf|png)$", full.names = TRUE)

    # First pass: handle plain <token>.{pdf,png} — copy to prefixed name.
    for (f in files) {
      ext  <- tools::file_ext(f)
      base <- tools::file_path_sans_ext(basename(f))
      if (!is.null(.parse_prefixed(base))) next  # already prefixed
      if (base %in% c("Fig1_composite", "Fig2_composite", "Fig3_composite",
                      "Fig4_composite", "Fig5_composite", "Fig6_composite"))
        next                                     # composite output, skip
      token  <- base
      prefix <- compute_panel_prefix(token, pa)
      if (is.null(prefix)) next                  # not in registry; skip silently

      target <- file.path(d, paste0(prefix, "_", token, ".", ext))
      if (file.exists(target) &&
          file.info(target)$mtime >= file.info(f)$mtime) {
        kept <- c(kept, target)
        next
      }
      file.copy(f, target, overwrite = TRUE)
      added <- c(added, target)
    }

    # Second pass — pdf/png: walk prefixed files and sweep stale ones.
    # Only sweep when the parsed token IS in the registry. Files that LOOK
    # like a prefix-token name but whose token isn't in the panel registry
    # (e.g. Fig{N}_composite for the figure composite, or the legacy
    # Fig1A_consort_flow CONSORT panel) are left alone.
    files <- list.files(d, pattern = "\\.(pdf|png)$", full.names = TRUE)
    for (f in files) {
      base   <- tools::file_path_sans_ext(basename(f))
      parsed <- .parse_prefixed(base)
      if (is.null(parsed)) next
      expected <- compute_panel_prefix(parsed$token, pa)
      if (is.null(expected)) next                     # token not in registry — leave it
      if (expected == parsed$full_prefix) next        # already correct
      file.remove(f)
      removed <- c(removed, f)
    }

    # Third pass — rds: sweep ANY prefix-named .rds file. The canonical
    # composer-input artefact is always <token>.rds (token-only, no
    # prefix); a Fig{N}{X}_<token>.rds or ED{NN}{X}_<token>.rds file is
    # always a stale leftover from a renamed-on-disk producer call.
    rds_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
    for (f in rds_files) {
      base   <- tools::file_path_sans_ext(basename(f))
      # KM-split files have __curve / __table suffix; strip before parsing.
      base_root <- sub("__(curve|table)$", "", base)
      parsed <- .parse_prefixed(base_root)
      if (is.null(parsed)) next                     # un-prefixed, leave alone
      file.remove(f)
      removed <- c(removed, f)
    }

    # Fourth pass — composite files: ED groups whose figure number changed
    # leave a stale EDNN_composite.{pdf,png} on disk (composer writes the
    # new EDMM_composite alongside; sweeper drops the old one). Identify
    # the dir's group via path, look up the current ED number, and
    # remove any EDNN_composite where NN doesn't match. Main-figure
    # composites (Fig{N}_composite) are unique per dir so they don't
    # need sweeping.
    composite_files <- list.files(d, pattern = "^ED[0-9]+_composite\\.(pdf|png)$",
                                  full.names = TRUE)
    if (length(composite_files) > 0L) {
      group <- basename(d)
      grp_num <- pa$group_number$EDFig[[group]]
      if (!is.null(grp_num)) {
        expected_stem <- sprintf("ED%02d_composite", grp_num)
        for (f in composite_files) {
          stem <- tools::file_path_sans_ext(basename(f))
          if (stem != expected_stem) {
            file.remove(f)
            removed <- c(removed, f)
          }
        }
      }
    }
  }

  if (verbose) {
    cat(sprintf("[sync_panel_prefixes] added/refreshed %d, removed %d stale, kept %d in-sync\n",
                length(added), length(removed), length(kept)))
    if (length(removed) > 0L) {
      cat("  Stale (removed):\n")
      for (f in removed) cat("    -", sub(paste0("^", here(), "/?"), "", f), "\n")
    }
    if (length(added) > 0L && length(added) <= 20L) {
      cat("  Refreshed:\n")
      for (f in added) cat("    +", sub(paste0("^", here(), "/?"), "", f), "\n")
    }
  }

  # Maintain per-group numbered symlinks under ExtendedData/numbered/.
  # 2026-04-27: previously these symlinks lived alongside the canonical
  # ed_<token> dirs at ExtendedData/, which mixed two views in one
  # parent. Splitting them into a sibling `numbered/` subdir keeps the
  # token-named (registry-stable) view at the top level and the
  # manuscript-numbered (resolver-driven) view in its own dir, so a
  # fresh manuscript render produces an updated numbered/ view without
  # cluttering the canonical layout.
  sym_added   <- character(0)
  sym_removed <- character(0)
  ed_root      <- here("results", "ExtendedData")
  numbered_dir <- file.path(ed_root, "numbered")
  if (dir.exists(ed_root)) {
    dir.create(numbered_dir, recursive = TRUE, showWarnings = FALSE)

    entries <- list.files(ed_root, full.names = TRUE)
    # Canonical token-named directories: ed_<token> (excluding _retired_*
    # and the numbered/ sibling itself).
    canonical <- entries[file.info(entries)$isdir &
                         grepl("/ed_[a-z][a-z0-9_]*$", entries) &
                         !grepl("/_retired_", entries)]
    for (canon in canonical) {
      token <- sub("^ed_", "", basename(canon))
      # Use the group's own panel_letter list to find any token in it,
      # then ask compute_panel_prefix for the group's number. We pick
      # the first panel under the group token; for single-panel groups
      # that IS the group itself.
      group_key <- paste0("ed_", token)
      panel_tokens <- names(pa$panel_letter[[group_key]])
      if (length(panel_tokens) == 0L) next                     # group not in registry
      prefix <- compute_panel_prefix(panel_tokens[1L], pa)
      if (is.null(prefix) || !startsWith(prefix, "ED")) next
      # Strip any panel-letter suffix from the prefix to get just the
      # group-level "ED{NN}" stem.
      ed_num_str <- sub("([A-Z])$", "", prefix)                # "ED02A" -> "ED02"
      # Lowercase for the directory name to match the token-style sibling.
      target_link <- file.path(numbered_dir,
                               sprintf("%s_%s", tolower(ed_num_str), token))
      # Symlink target is relative to numbered/ (one level deeper than
      # the canonical sibling), so prefix `..` to climb out before
      # naming the sibling.
      rel_target <- file.path("..", basename(canon))            # e.g. "../ed_prisma"
      if (file.exists(target_link) || .Platform$OS.type == "unix") {
        info <- tryCatch(file.info(target_link), error = function(e) NULL)
        is_symlink <- !is.null(info) && !is.na(Sys.readlink(target_link)) &&
                      nzchar(Sys.readlink(target_link))
        if (is_symlink) {
          current <- Sys.readlink(target_link)
          if (identical(current, rel_target)) next             # already correct
          file.remove(target_link)
        } else if (file.exists(target_link)) {
          # Not a symlink — leave it alone (could be a real dir created
          # by another process). Skip rather than risk data loss.
          next
        }
      }
      ok <- file.symlink(rel_target, target_link)
      if (isTRUE(ok)) sym_added <- c(sym_added, target_link)
    }

    # Sweep stale numbered symlinks under numbered/ — ed{NN}_<token>
    # whose number doesn't match the current resolver-assigned ED
    # number for the same token. Also picks up legacy symlinks that may
    # still live at the ExtendedData/ root from before the numbered/
    # split (those get removed unconditionally because they no longer
    # belong at that level).
    legacy_at_root <- entries[grepl("/ed[0-9]+_[a-z][a-z0-9_]*$", entries)]
    for (sym in legacy_at_root) {
      readout <- Sys.readlink(sym)
      if (!nzchar(readout)) next                            # not a symlink
      file.remove(sym); sym_removed <- c(sym_removed, sym)
    }
    nested <- list.files(numbered_dir, full.names = TRUE)
    candidates <- nested[grepl("/ed[0-9]+_[a-z][a-z0-9_]*$", nested)]
    for (sym in candidates) {
      readout <- Sys.readlink(sym)
      if (!nzchar(readout)) next                            # not a symlink
      # Recover the token from "ed{NN}_<token>".
      token <- sub("^ed[0-9]+_", "", basename(sym))
      group_key <- paste0("ed_", token)
      panel_tokens <- names(pa$panel_letter[[group_key]])
      if (length(panel_tokens) == 0L) {
        # Token no longer in registry — remove the stale symlink.
        file.remove(sym); sym_removed <- c(sym_removed, sym); next
      }
      prefix <- compute_panel_prefix(panel_tokens[1L], pa)
      if (is.null(prefix)) {
        file.remove(sym); sym_removed <- c(sym_removed, sym); next
      }
      ed_num_str <- sub("([A-Z])$", "", prefix)
      expected <- sprintf("%s_%s", tolower(ed_num_str), token)
      if (basename(sym) != expected) {
        file.remove(sym); sym_removed <- c(sym_removed, sym)
      }
    }
  }

  if (verbose && (length(sym_added) > 0L || length(sym_removed) > 0L)) {
    cat(sprintf("[sync_panel_prefixes] ED group symlinks: +%d, -%d\n",
                length(sym_added), length(sym_removed)))
    for (f in sym_added)   cat("    +", sub(paste0("^", here(), "/?"), "", f), "\n")
    for (f in sym_removed) cat("    -", sub(paste0("^", here(), "/?"), "", f), "\n")
  }

  invisible(list(added = added, removed = removed, kept = kept,
                 sym_added = sym_added, sym_removed = sym_removed))
}

# ── CLI entry ────────────────────────────────────────────────────────────────

.this_script <- "sync_panel_prefixes.R"
.invoked_as_script <- any(grepl(
  paste0("--file=.*/", .this_script, "$|--file=", .this_script, "$"),
  commandArgs(trailingOnly = FALSE)
))
if (.invoked_as_script) {
  sync_panel_prefixes()
}
